#!/usr/bin/env python3
"""Neisseria meningitidis typing (MLST, BAST/MenDeVAR, finetyping) through the PubMLST REST API.

Called by the WDL task `neisseria_typing`.
"""
import argparse
import base64
import http.client
import json
import logging
import os
import random
import sys
import tarfile
import time
import urllib.error
import urllib.request
from datetime import date

import pandas as pd

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
EMPTY_VALUE = '?'      # PubMLST answered but found nothing
API_ERROR = 'api_error'   # PubMLST could not be queried (HTTP error, timeout...): see the log for the type
NADA_NOVEL = 'novel'   # nadA hit (closest non-exact match) but no exact allele
NADA_MIN_IDENTITY = 90.0   # % protein identity: a real nadA variant is expected well above 90% (K. Jolley)
API_ROOT = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef'
API_TIMEOUT = 300     # seconds without any response (a BLAST on a whole genome can be slow)
API_RETRIES = 3
API_KEY = None   # personal PubMLST API key, sent as the X-API-Key header (set in main)

SCHEME_ID = {'mlst': '1', 'finetyping': '2', 'bast': '53'}

# Loci of each typing scheme (order = order of the values in the "profile" columns)
MLST_LOCI = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']
BAST_LOCI = ['fHbp_peptide', 'NHBA_peptide', 'NadA_peptide', 'PorA_VR1', 'PorA_VR2']
FINETYPING_LOCI = ['PorA_VR1', 'PorA_VR2', 'FetA_VR']

# Summary columns: column name -> loci joined with commas
PROFILES = {
    'mlst_profile': MLST_LOCI,
    'bast_profile': BAST_LOCI,
    'finetyping_profile': FINETYPING_LOCI,
}

# All loci (no duplicates, order kept)
ALL_LOCI = list(dict.fromkeys(MLST_LOCI + BAST_LOCI + FINETYPING_LOCI))

# "fields" returned by PubMLST -> report columns
FIELD_MAP = {
    'mlst': {'ST': 'st', 'clonal_complex': 'clonal_complex'},
    'bast': {
        'BAST': 'bast_type',
        'MenDeVAR_Bexsero_reactivity': 'bexsero_cross_reactivity',
        'MenDeVAR_Trumenba_reactivity': 'trumenba_cross_reactivity',
    },
}
BAST_COLUMNS = list(FIELD_MAP['bast'].values())

# CSV column order (the summary columns come last so existing columns do not shift)
COLUMNS = (
    ['date', 'sample', 'st', 'clonal_complex']
    + MLST_LOCI + ['FetA_VR']
    + ['bast_type', 'fHbp_peptide', 'NHBA_peptide', 'NadA_peptide', 'PorA_VR1', 'PorA_VR2',
       'bexsero_cross_reactivity', 'trumenba_cross_reactivity']
    + list(PROFILES)
)

# Queries that can explain a missing value of a column (used to tell "api_error" from "?")
RELATED_QUERIES = {}
for _locus in MLST_LOCI:
    RELATED_QUERIES[_locus] = ['mlst', f'locus:{_locus}']
RELATED_QUERIES['st'] = RELATED_QUERIES['clonal_complex'] = ['mlst', 'mlst_designations']
RELATED_QUERIES['FetA_VR'] = ['finetyping', 'locus:FetA_VR']
for _locus in ('PorA_VR1', 'PorA_VR2'):
    RELATED_QUERIES[_locus] = ['bast', 'finetyping', f'locus:{_locus}']
for _locus in ('fHbp_peptide', 'NHBA_peptide', 'NadA_peptide'):
    RELATED_QUERIES[_locus] = ['bast', f'locus:{_locus}']
for _column in BAST_COLUMNS:
    RELATED_QUERIES[_column] = ['bast', 'bast_designations'] + [f'locus:{_l}' for _l in BAST_LOCI]

log = logging.getLogger('neisseria_typing')


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def is_missing(value):
    return value in (None, '', EMPTY_VALUE)


def save_json(path, content):
    with open(path, 'w') as f:
        json.dump(content, f, indent=2)


def read_assembly_b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


class ApiError(Exception):
    """A PubMLST query failed. `reason` is short and goes to the log/summary: 'HTTP 404', 'timeout'..."""

    def __init__(self, reason, detail='', body=None):
        super().__init__(f'{reason} {detail}'.strip())
        self.reason = reason
        self.detail = detail
        self.body = body      # parsed error body, kept in the JSON files for inspection


class ApiAuthError(ApiError):
    """401/403: every following query would fail the same way."""


def classify(err):
    """Short, readable type of a network-level error."""
    if isinstance(err, TimeoutError) or isinstance(getattr(err, 'reason', None), TimeoutError):
        return 'timeout'
    if isinstance(err, ValueError):
        return 'bad response'
    return 'network'


def wait_before_retry(attempt, url, reason):
    """Log the failed attempt and sleep (growing delay + jitter so parallel jobs do not retry in sync)."""
    delay = 10 * attempt + random.uniform(0, 5)
    log.warning(f"Attempt {attempt}/{API_RETRIES} failed on {url.replace(API_ROOT, '')} "
                f"({reason}); retrying in {delay:.0f}s")
    time.sleep(delay)


def request_headers():
    headers = {'Content-Type': 'application/json'}
    if API_KEY:
        headers['X-API-Key'] = API_KEY
    return headers


def post_json(url, payload):
    """POST JSON. Returns the parsed response, or raises ApiError (ApiAuthError for 401/403).

    Network errors, timeouts, 5xx and 429 are retried (growing delay + jitter); any other
    4xx (404, 451...) fails immediately. Nothing is ever turned into an empty result silently.
    """
    body = json.dumps(payload).encode()
    where = url.replace(API_ROOT, '')
    last = None

    for attempt in range(1, API_RETRIES + 1):
        request = urllib.request.Request(url, data=body, method='POST', headers=request_headers())
        try:
            with urllib.request.urlopen(request, timeout=API_TIMEOUT) as response:
                return json.load(response)
        except urllib.error.HTTPError as err:
            raw = err.read()
            try:
                parsed = json.loads(raw)
            except ValueError:
                parsed = raw[:300].decode('utf-8', 'replace')
            reason, detail = f'HTTP {err.code}', f'{where}: {raw[:200]!r}'
            if err.code in (401, 403):
                raise ApiAuthError(reason, detail, parsed)
            if err.code < 500 and err.code != 429:
                raise ApiError(reason, detail, parsed)
            last = ApiError(reason, detail, parsed)
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                http.client.HTTPException, ValueError) as err:
            last = ApiError(classify(err), f'{where}: {err}')

        if attempt < API_RETRIES:
            wait_before_retry(attempt, url, last.reason)

    raise last


def sequence_url(kind, locus=None):
    if kind == 'locus':
        return f'{API_ROOT}/loci/{locus}/sequence'
    return f'{API_ROOT}/schemes/{SCHEME_ID[kind]}/sequence'


def query_pubmlst(url, sequence_b64, out_file):
    """Send the assembly to PubMLST and keep the raw response (or the error) on disk."""
    try:
        result = post_json(url, {'base64': True, 'sequence': sequence_b64})
    except ApiError as err:
        save_json(out_file, {'error': err.reason, 'detail': err.detail, 'body': err.body})
        raise
    save_json(out_file, result)
    return result if isinstance(result, dict) else {}


# --------------------------------------------------------------------------- #
# Reading the responses
# --------------------------------------------------------------------------- #
def store_exact_matches(result, data, locus=None):
    """Copy the exact allele matches into `data`. Returns True if there were any.

    PubMLST returns two shapes:
      - scheme query: {"abcZ": [{"allele_id": "2"}], "adk": [...]}  (dict)
      - locus query:  [{"allele_id": "2", ...}]                      (list, `locus` required)
    """
    exact = result.get('exact_matches')
    if not exact:
        return False

    if isinstance(exact, dict):
        pairs = [(name, hits[0]) for name, hits in exact.items() if hits]
    elif locus:
        pairs = [(locus, exact[0])]
    else:
        return False

    found = False
    for name, hit in pairs:
        if hit.get('allele_id') is not None:
            data[name] = str(hit['allele_id'])
            found = True
    return found


def store_fields(result, kind, data):
    """Copy the scheme fields (ST, CC, BAST, MenDeVAR...) into `data`."""
    fields = result['fields']
    for db_key, column in FIELD_MAP[kind].items():
        data[column] = fields.get(db_key, EMPTY_VALUE)


def type_scheme(kind, sequence_b64, out_file, data):
    """Query a whole scheme (mlst / bast / finetyping). Returns True if identified."""
    result = query_pubmlst(sequence_url(kind), sequence_b64, out_file)

    if not store_exact_matches(result, data):
        return False
    if 'fields' in result and kind in FIELD_MAP:
        store_fields(result, kind, data)
    return True


def lookup_designations(kind, loci, data, out_dir):
    """Ask scheme `kind` for the fields (ST, BAST...) matching the alleles already known.

    Returns {column: value} (empty if PubMLST does not know this profile: an HTTP 404 here means
    "profile not found" and is not treated as an error).
    """
    designations = {
        locus: [{'allele': str(data[locus])}] for locus in loci if not is_missing(data.get(locus))
    }
    out_file = os.path.join(out_dir, f"{data['sample']}_{kind}_type.json")
    try:
        result = post_json(f"{API_ROOT}/schemes/{SCHEME_ID[kind]}/designations",
                           {'designations': designations})
    except ApiError as err:
        save_json(out_file, {'error': err.reason, 'detail': err.detail, 'body': err.body})
        if err.reason != 'HTTP 404':
            raise
        log.info(f'[{kind} designations] HTTP 404: profile not found ({err.detail})')
        result = {}
    else:
        save_json(out_file, result)

    found = {}
    if isinstance(result, dict) and 'fields' in result:
        store_fields(result, kind, found)
    return found


def nada_is_novel(has_best_match, identity):
    """A best match is "novel" only if its protein identity reaches the threshold (default 90%)."""
    if not has_best_match:
        return False
    if NADA_MIN_IDENTITY is None:
        return True
    try:
        return float(identity) >= NADA_MIN_IDENTITY
    except (TypeError, ValueError):
        return False


def build_profile(data, loci):
    """'2,3,4,3,8,4,6' ; a missing locus is '?' (or 'api_error') ; everything missing -> '?' / 'api_error'."""
    values = []
    for locus in loci:
        value = data.get(locus)
        values.append(API_ERROR if value == API_ERROR else EMPTY_VALUE if is_missing(value) else str(value))
    if all(value in (EMPTY_VALUE, API_ERROR) for value in values):
        return API_ERROR if API_ERROR in values else EMPTY_VALUE
    return ','.join(values)


# --------------------------------------------------------------------------- #
# Typing of one sample
# --------------------------------------------------------------------------- #
def type_sample(assembly, sample, output_dir, output_filename, split):
    """Type one sample. Returns {query: reason} for every PubMLST query that failed."""
    data = {'date': date.today().strftime('%d-%m-%Y'), 'sample': sample}
    best_matches = {}
    best_identity = {}
    errors = {}      # query name -> short reason ('HTTP 404', 'timeout'...)
    blocked = []     # filled by an authentication error: all remaining queries are skipped
    sequence_b64 = read_assembly_b64(assembly)

    def attempt(name, func, *args):
        """Run one PubMLST query. On failure: log it, remember it in `errors`, return None."""
        if blocked:
            errors[name] = blocked[0]
            log.error(f'✖️ [{name}] skipped ({blocked[0]} on an earlier query)')
            return None
        try:
            return func(*args)
        except ApiAuthError as err:
            errors[name] = err.reason
            blocked.append(err.reason)
            log.error(f'✖️ [{name}] {err}')
            log.error('PubMLST refused the request: a valid personal API key is required '
                      '(--api-key-file or PUBMLST_API_KEY). The remaining queries are skipped.')
        except ApiError as err:
            errors[name] = err.reason
            log.error(f'✖️ [{name}] {err}')
        except Exception as err:   # unexpected problem: the sample must still get its outputs and log
            errors[name] = f'error ({type(err).__name__})'
            log.error(f'✖️ [{name}] unexpected error: {err!r}')
        return None

    log.info(f'\n---------------------\n    {sample}\n---------------------')

    # 1) The three schemes
    for kind in ('mlst', 'bast', 'finetyping'):
        log.info(f'[Running] {kind}')
        out_file = os.path.join(output_dir, f'{sample}_{kind}.json')
        identified = attempt(kind, type_scheme, kind, sequence_b64, out_file, data)
        if kind not in errors:
            log.info(f'✔️ {kind}' if identified else f'✖️ {kind}: PubMLST returned no exact match')

    # 2) Missing loci: query them one by one
    #    (unless the 3 schemes failed: PubMLST cannot be queried, no point insisting locus by locus)
    missing_loci = [locus for locus in ALL_LOCI if is_missing(data.get(locus))]
    if missing_loci and all(kind in errors for kind in ('mlst', 'bast', 'finetyping')):
        log.error('PubMLST could not be queried: skipping the locus-by-locus queries')
        missing_loci = []
    if missing_loci:
        log.info(f'Missing {missing_loci}...')

    for locus in missing_loci:
        name = f'locus:{locus}'
        out_file = os.path.join(output_dir, f'{sample}_{locus}.json')
        result = attempt(name, query_pubmlst, sequence_url('locus', locus), sequence_b64, out_file)
        if name in errors:
            continue

        best = result.get('best_match')
        if store_exact_matches(result, data, locus):
            log.info(f'Fetched {locus}!')
        elif isinstance(best, dict):
            best_matches[locus] = best.get('allele_id')
            best_identity[locus] = best.get('identity')
            log.info(f"{locus}: no exact match, best match allele {best.get('allele_id')} "
                     f"(identity {best.get('identity')}%)")
        else:
            log.info(f'{locus}: no match returned by PubMLST')

    # 3) NadA without an exact allele (only if PubMLST actually answered). The BAST scheme response is trusted
    #    first (it already carries NadA "0" when PubMLST found nothing even at low stringency); this step only
    #    applies when NadA is missing from that response, using the NadA locus query (K. Jolley):
    #    - best match with protein identity >= 90% (any new variant differs by 1-2 aa from a known one) -> "novel"
    #    - nothing, or a weaker hit (e.g. 28% identity with an unrelated protein)                        -> "0"
    #    With "novel" the BAST is NOT looked up: the closest BAST is not meaningful (one differing protein
    #    variant can change the predicted reactivity completely).
    #    The threshold can be changed with --nada-novel-min-identity (0 = any best match is "novel").
    #    If PubMLST could not be queried for NadA, the value stays missing (-> "api_error", never "0").
    if (is_missing(data.get('NadA_peptide'))
            and 'bast' not in errors and 'locus:NadA_peptide' not in errors):
        novel = nada_is_novel('NadA_peptide' in best_matches, best_identity.get('NadA_peptide'))
        data['NadA_peptide'] = NADA_NOVEL if novel else '0'
        if not novel and 'NadA_peptide' in best_matches:
            log.info(f"NadA: best match identity {best_identity.get('NadA_peptide')}% is below "
                     f"{NADA_MIN_IDENTITY}%: treated as absent (0)")
        log.info(f"NadA: {data['NadA_peptide']}")

    # 4) BAST / MenDeVAR are sometimes missing from the first response -> confirmation request
    if any(is_missing(data.get(column)) for column in BAST_COLUMNS):
        if data.get('NadA_peptide') == NADA_NOVEL:
            log.info('[BAST] NadA is novel (no exact allele): BAST cannot be assigned, the closest BAST is not reported')
        elif any(is_missing(data.get(locus)) for locus in BAST_LOCI):
            log.info('[BAST] some BAST loci are missing: BAST cannot be assigned (no confirmation request sent)')
        else:
            log.info('[Confirm BAST type and MenDeVar]')
            found = attempt('bast_designations', lookup_designations, 'bast', BAST_LOCI, data, output_dir)
            if 'bast_designations' not in errors:
                # only fill what is missing: a value already found is never overwritten by '?'
                for column in BAST_COLUMNS:
                    if is_missing(data.get(column)):
                        data[column] = found.get(column, EMPTY_VALUE)
                status = {'Bast Type': 'bast_type', 'Bexsero': 'bexsero_cross_reactivity',
                          'Trumenba': 'trumenba_cross_reactivity'}
                log.info(', '.join(f"{'✔️' if not is_missing(data.get(col)) else '✖️'} {label}"
                                   for label, col in status.items()))

    # 5) Full MLST profile but no ST: ask PubMLST before concluding it is a new strain
    if is_missing(data.get('st')) and not any(is_missing(data.get(locus)) for locus in MLST_LOCI):
        found = attempt('mlst_designations', lookup_designations, 'mlst', MLST_LOCI, data, output_dir)
        if 'mlst_designations' not in errors:      # on error the ST stays unknown, it is not "new"
            if is_missing(found.get('st')):
                log.info('- [MLST] Full profile but no ST assignment; it must be new!')
                data['st'] = 'new'
            else:
                log.info(f"- [MLST] ST found from the full profile: {found['st']}")
                data['st'] = found['st']
                data['clonal_complex'] = found.get('clonal_complex', EMPTY_VALUE)

    # 6) Values still missing: "api_error" if a related query failed ("?" = PubMLST answered, found nothing)
    for column, queries in RELATED_QUERIES.items():
        if is_missing(data.get(column)) and any(query in errors for query in queries):
            data[column] = API_ERROR
    for column, loci in PROFILES.items():
        data[column] = build_profile(data, loci)
    for column in COLUMNS:
        if is_missing(data.get(column)):
            data[column] = EMPTY_VALUE

    log.info('\n===== Summary =====')
    if errors:
        log.error('PubMLST queries that failed: ' + ', '.join(f'{query} = {reason}' for query, reason in errors.items()))
        flagged = [column for column in COLUMNS if data[column] == API_ERROR]
        log.error(f"Columns reported as '{API_ERROR}': {', '.join(flagged) or 'none'}")
    else:
        log.info('All PubMLST queries were answered')
    log.info(f'{sample}: {data}')

    # 7) Outputs
    row = {column: data[column] for column in COLUMNS}

    if split:  # one .txt file per column, read by the WDL
        for column, value in row.items():
            with open(os.path.join(output_dir, f'{column}.txt'), 'w') as f:
                f.write(f'{value}\n')

    pd.DataFrame([row]).to_csv(os.path.join(output_dir, output_filename), index=False)
    pd.DataFrame([best_matches]).to_csv(os.path.join(output_dir, f'{sample}_best_matches.csv'), index=False)
    return errors


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args():
    parser = argparse.ArgumentParser(description='Neisseria typing (MLST, BAST, finetyping) via PubMLST')
    parser.add_argument('--input', required=True, help='Assembly FASTA file')
    parser.add_argument('--output', dest='output_dir', required=True, help='Destination directory')
    parser.add_argument('--output_csv_file', dest='output_file', required=True, help='Typing report CSV filename')
    parser.add_argument('--log-file', dest='log_file', default=None,
                        help='Log filename in the output directory (default: <sample>_logs.txt)')
    parser.add_argument('--strict', action='store_true',
                        help='Exit with an error if any PubMLST query failed (default: always exit 0 once the '
                             "outputs are written, failed values are reported as 'api_error')")
    parser.add_argument('--sample', default=None,
                        help='Sample name (default: input filename up to the first dot)')
    parser.add_argument('--api-key-file', default=None,
                        help='Text file containing your personal PubMLST API key '
                             '(else the PUBMLST_API_KEY environment variable is used)')
    parser.add_argument('--nada-novel-min-identity', type=float, default=NADA_MIN_IDENTITY,
                        help="A nadA best match counts as 'novel' (BAST then skipped) only when its %% identity is "
                             "at least this value, otherwise NadA is '0'. Default: 90; "
                             "0 = any best match is 'novel'")
    parser.add_argument('--json-archive', default=None,
                        help='Name of the tar.gz bundling all raw PubMLST JSON responses and the log '
                             '(default: <sample>_api_json.tar.gz, written in the output directory)')
    parser.add_argument('--split', action='store_true', help='Also write one .txt file per column')
    parser.add_argument('--files', help='Deprecated and ignored (kept so existing WDL calls still work)')
    return parser.parse_args()


def setup_logging(log_path):
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_path, mode='a', encoding='utf-8')],
    )


def load_api_key(path):
    """API key from --api-key-file, else from PUBMLST_API_KEY. Never logged."""
    key = os.environ.get('PUBMLST_API_KEY', '')
    if path:
        if not os.path.isfile(path):
            sys.exit(f'API key file not found: {path}')
        with open(path) as f:
            key = f.read()
    return key.strip() or None


def archive_json(output_dir, sample, log_file, name=None):
    """Bundle every raw PubMLST response (*.json) and the log in one tar.gz, to inspect what the API returned."""
    name = name or f'{sample}_api_json.tar.gz'
    files = sorted(f for f in os.listdir(output_dir) if f.endswith('.json'))
    if os.path.isfile(os.path.join(output_dir, log_file)):
        files.append(log_file)
    with tarfile.open(os.path.join(output_dir, name), 'w:gz') as tar:
        for filename in files:
            # every file inside the archive carries the sample name (e.g. logs.txt -> <sample>_logs.txt)
            prefixed = filename if filename.startswith(f'{sample}_') else f'{sample}_{filename}'
            tar.add(os.path.join(output_dir, filename), arcname=f'{sample}_api_json/{prefixed}')
    log.info(f'Raw API responses archived in {name} ({len(files)} files)')


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f'Input not valid (expected one fasta file): {args.input}')

    global API_KEY, NADA_MIN_IDENTITY
    API_KEY = load_api_key(args.api_key_file)
    NADA_MIN_IDENTITY = args.nada_novel_min_identity

    sample = args.sample or os.path.basename(args.input).split('.')[0]
    log_file = args.log_file or f'{sample}_logs.txt'

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(os.path.join(args.output_dir, log_file))
    log.info('PubMLST API key: ' + ('provided' if API_KEY else 'NOT provided (sequence queries will be refused)'))

    errors = {}
    try:
        errors = type_sample(args.input, sample, args.output_dir, args.output_file, args.split)
    finally:   # also written when the run stops on an error: the raw responses are what you want then
        archive_json(args.output_dir, sample, log_file, args.json_archive)

    if errors and args.strict:
        sys.exit('PubMLST queries failed (see the log): ' + ', '.join(f'{q} = {r}' for q, r in errors.items()))


if __name__ == '__main__':
    main()

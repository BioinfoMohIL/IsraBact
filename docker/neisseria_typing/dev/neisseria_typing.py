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
import time
import urllib.error
import urllib.request
from datetime import date

import pandas as pd

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
EMPTY_VALUE = '?'      # unknown / not retrieved
NADA_NOVEL = 'novel?'  # only with --nada-novel-min-identity: strong nadA hit without an exact allele
NADA_MIN_IDENTITY = None   # None (default) = no exact nadA allele -> '0', like PubMLST's own BAST profiles
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


def wait_before_retry(attempt, url, reason):
    """Log the failed attempt and sleep (growing delay + jitter so parallel jobs do not retry in sync)."""
    delay = 10 * attempt + random.uniform(0, 5)
    log.warning(f"Attempt {attempt}/{API_RETRIES} failed on {url.replace(API_ROOT, '')} "
                f"({reason}); retrying in {delay:.0f}s")
    time.sleep(delay)


def fail_auth(code, detail):
    """Authentication problems stop the whole run: every following query would fail the same way."""
    log.error(f'HTTP {code} from PubMLST: {detail}')
    sys.exit(f'PubMLST refused the request (HTTP {code}): a valid personal API key is required. '
             'Provide it with --api-key-file or the PUBMLST_API_KEY environment variable. Nothing was typed.')


def request_headers():
    headers = {'Content-Type': 'application/json'}
    if API_KEY:
        headers['X-API-Key'] = API_KEY
    return headers


def post_json(url, payload):
    """POST JSON, retried on network errors, 5xx and 429 (with a jittered, growing delay).

    401/403 (authentication required or refused) stop the whole run with a clear message:
    every following query would fail the same way, and "?" everywhere must never look like a result.
    Any other 4xx response (e.g. 404 = no match) is returned as is, but logged as a warning
    so that a silent API error is never mistaken for "no result".
    If its body is not JSON, {} is returned.
    """
    body = json.dumps(payload).encode()
    last_error = None

    for attempt in range(1, API_RETRIES + 1):
        request = urllib.request.Request(
            url, data=body, method='POST', headers=request_headers())
        try:
            with urllib.request.urlopen(request, timeout=API_TIMEOUT) as response:
                return json.load(response)
        except urllib.error.HTTPError as err:
            raw = err.read()
            if err.code in (401, 403):
                fail_auth(err.code, f"{url.replace(API_ROOT, '')} {raw[:200]!r}")
            if err.code < 500 and err.code != 429:
                log.warning(f"HTTP {err.code} from PubMLST ({url.replace(API_ROOT, '')}): {raw[:200]!r}")
                try:
                    return json.loads(raw)
                except ValueError:
                    return {}
            last_error, reason = err, f'HTTP {err.code}'
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                http.client.HTTPException, ValueError) as err:
            last_error, reason = err, f'{type(err).__name__}: {err}'

        if attempt < API_RETRIES:
            wait_before_retry(attempt, url, reason)

    raise last_error


def sequence_url(kind, locus=None):
    if kind == 'locus':
        return f'{API_ROOT}/loci/{locus}/sequence'
    return f'{API_ROOT}/schemes/{SCHEME_ID[kind]}/sequence'


def query_pubmlst(url, sequence_b64, out_file):
    """Send the assembly to PubMLST and keep the raw response on disk."""
    result = post_json(url, {'base64': True, 'sequence': sequence_b64})
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

    Returns {column: value} (empty if PubMLST does not know this profile).
    """
    designations = {
        locus: [{'allele': str(data[locus])}] for locus in loci if not is_missing(data.get(locus))
    }
    result = post_json(f"{API_ROOT}/schemes/{SCHEME_ID[kind]}/designations",
                       {'designations': designations})
    save_json(os.path.join(out_dir, f'{kind}_type.json'), result)

    found = {}
    if isinstance(result, dict) and 'fields' in result:
        store_fields(result, kind, found)
    return found


def nada_is_novel(identity):
    """True only when the opt-in threshold is set and the best match reaches it."""
    if NADA_MIN_IDENTITY is None:
        return False
    try:
        return float(identity) >= NADA_MIN_IDENTITY
    except (TypeError, ValueError):
        return False


def build_profile(data, loci):
    """'2,3,4,3,8,4,6' ; a missing locus is '?' ; everything missing -> '?'."""
    values = [EMPTY_VALUE if is_missing(data.get(locus)) else str(data[locus]) for locus in loci]
    if all(value == EMPTY_VALUE for value in values):
        return EMPTY_VALUE
    return ','.join(values)


# --------------------------------------------------------------------------- #
# Typing of one sample
# --------------------------------------------------------------------------- #
def type_sample(assembly, sample, output_dir, output_filename, split):
    data = {'date': date.today().strftime('%d-%m-%Y'), 'sample': sample}
    best_matches = {}
    failed = set()
    sequence_b64 = read_assembly_b64(assembly)

    log.info(f'\n---------------------\n    {sample}\n---------------------')

    # 1) The three schemes
    for kind in ('mlst', 'bast', 'finetyping'):
        log.info(f'[Running] {kind}')
        out_file = os.path.join(output_dir, f'{sample}_{kind}.json')
        try:
            identified = type_scheme(kind, sequence_b64, out_file, data)
            log.info(f'✔️ {kind}' if identified else f'✖️ Cannot fetch {kind}')
        except Exception as err:
            failed.add(kind)
            log.error(f'Error in fetching {kind}: {err}')

    # 2) Missing loci: query them one by one
    #    (unless the 3 schemes failed: PubMLST is unreachable, no point insisting locus by locus)
    missing_loci = [locus for locus in ALL_LOCI if is_missing(data.get(locus))]
    if missing_loci and len(failed) == 3:
        log.error('PubMLST unreachable: skipping locus-by-locus queries')
        missing_loci = []
    if missing_loci:
        log.info(f'Missing {missing_loci}...')

    failed_loci = set()
    best_identity = {}
    for locus in missing_loci:
        out_file = os.path.join(output_dir, f'{sample}_{locus}.json')
        try:
            result = query_pubmlst(sequence_url('locus', locus), sequence_b64, out_file)
        except Exception as err:
            failed_loci.add(locus)
            log.error(f'✖️ Cannot fetch the missing {locus}: {err}')
            continue

        best = result.get('best_match')
        if store_exact_matches(result, data, locus):
            log.info(f'Fetched {locus}!')
        elif isinstance(best, dict):
            best_matches[locus] = best.get('allele_id')
            best_identity[locus] = best.get('identity')
            log.info(f"{locus}: no exact match, best match allele {best.get('allele_id')} "
                     f"(identity {best.get('identity')}%)")

    # 3) NadA without an exact allele (only if PubMLST actually answered) -> "0" (needed to look up the BAST).
    #    The API also returns weak "best matches" (e.g. paralogous adhesins) that PubMLST does not count as
    #    nadA: isolates with the same fHbp/NHBA/PorA profile get NadA 0 and a BAST in PubMLST. So a best
    #    match alone never changes the result. Opt-in: --nada-novel-min-identity N marks "novel?" only when
    #    the best match has an identity >= N (BAST is then not looked up).
    if (is_missing(data.get('NadA_peptide'))
            and 'bast' not in failed and 'NadA_peptide' not in failed_loci):
        data['NadA_peptide'] = NADA_NOVEL if nada_is_novel(best_identity.get('NadA_peptide')) else '0'

    # 4) BAST / MenDeVAR are sometimes missing from the first response -> confirmation request
    if any(is_missing(data.get(column)) for column in BAST_COLUMNS):
        if data.get('NadA_peptide') == NADA_NOVEL:
            log.info('[BAST] NadA has no exact allele (novel?): BAST cannot be assigned')
        else:
            log.info('[Confirm BAST type and MenDeVar]')
            try:
                found = lookup_designations('bast', BAST_LOCI, data, output_dir)
                # only fill what is missing: a value already found is never overwritten by '?'
                for column in BAST_COLUMNS:
                    if is_missing(data.get(column)):
                        data[column] = found.get(column, EMPTY_VALUE)
            except Exception as err:
                log.error(f'✖️ Cannot confirm BAST: {err}')

            status = {'Bast Type': 'bast_type', 'Bexsero': 'bexsero_cross_reactivity',
                      'Trumenba': 'trumenba_cross_reactivity'}
            log.info(', '.join(f"{'✔️' if not is_missing(data.get(col)) else '✖️'} {label}"
                               for label, col in status.items()))

    # 5) Full MLST profile but no ST: ask PubMLST before concluding it is a new strain
    if is_missing(data.get('st')) and not any(is_missing(data.get(locus)) for locus in MLST_LOCI):
        try:
            found = lookup_designations('mlst', MLST_LOCI, data, output_dir)
        except Exception as err:
            log.error(f'✖️ Cannot confirm ST: {err}')   # unknown, not "new"
        else:
            if is_missing(found.get('st')):
                log.info('- [MLST] Full profile but no ST assignment; it must be new!')
                data['st'] = 'new'
            else:
                log.info(f"- [MLST] ST found from the full profile: {found['st']}")
                data['st'] = found['st']
                data['clonal_complex'] = found.get('clonal_complex', EMPTY_VALUE)

    # 6) Summary columns, then empty values -> '?'
    for column, loci in PROFILES.items():
        data[column] = build_profile(data, loci)

    for column in COLUMNS:
        if is_missing(data.get(column)):
            data[column] = EMPTY_VALUE

    log.info(f'{sample}: {data}')

    # 7) Outputs
    row = {column: data[column] for column in COLUMNS}

    if split:  # one .txt file per column, read by the WDL
        for column, value in row.items():
            with open(os.path.join(output_dir, f'{column}.txt'), 'w') as f:
                f.write(f'{value}\n')

    pd.DataFrame([row]).to_csv(os.path.join(output_dir, output_filename), index=False)
    pd.DataFrame([best_matches]).to_csv(os.path.join(output_dir, f'{sample}_best_matches.csv'), index=False)

    # PubMLST could not be queried at all: the files above are empty, so fail loudly
    # (a row full of "?" must never look like a normal result)
    if len(failed) == 3:
        sys.exit('PubMLST could not be queried (the 3 schemes failed): results are empty.')


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args():
    parser = argparse.ArgumentParser(description='Neisseria typing (MLST, BAST, finetyping) via PubMLST')
    parser.add_argument('--input', required=True, help='Assembly FASTA file')
    parser.add_argument('--output', dest='output_dir', required=True, help='Destination directory')
    parser.add_argument('--output_csv_file', dest='output_file', required=True, help='Typing report CSV filename')
    parser.add_argument('--log-file', dest='log_file', default='logs.txt', help='Log filename (in the output directory)')
    parser.add_argument('--sample', default=None,
                        help='Sample name (default: input filename up to the first dot)')
    parser.add_argument('--api-key-file', default=None,
                        help='Text file containing your personal PubMLST API key '
                             '(else the PUBMLST_API_KEY environment variable is used)')
    parser.add_argument('--nada-novel-min-identity', type=float, default=None,
                        help="Opt-in: report NadA as 'novel?' (and skip the BAST lookup) when the best nadA "
                             "match has at least this %% identity. Default: off (no exact allele -> '0')")
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


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f'Input not valid (expected one fasta file): {args.input}')

    global API_KEY, NADA_MIN_IDENTITY
    API_KEY = load_api_key(args.api_key_file)
    NADA_MIN_IDENTITY = args.nada_novel_min_identity

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(os.path.join(args.output_dir, args.log_file))
    log.info('PubMLST API key: ' + ('provided' if API_KEY else 'NOT provided (sequence queries will be refused)'))

    sample = args.sample or os.path.basename(args.input).split('.')[0]
    type_sample(args.input, sample, args.output_dir, args.output_file, args.split)


if __name__ == '__main__':
    main()
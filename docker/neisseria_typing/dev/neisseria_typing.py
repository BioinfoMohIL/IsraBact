#!/usr/bin/env python3
"""Typage de Neisseria meningitidis (MLST, BAST/MenDeVAR, finetyping) via l'API REST de PubMLST.

Appelé par la tâche WDL `neisseria_typing`.
"""
import argparse
import base64
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import date

import pandas as pd

# --------------------------------------------------------------------------- #
# Constantes
# --------------------------------------------------------------------------- #
EMPTY_VALUE = '?'
API_ROOT = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef'
API_TIMEOUT = 60      # secondes
API_RETRIES = 3

SCHEME_ID = {'mlst': '1', 'finetyping': '2', 'bast': '53'}

# Loci de chaque typage (ordre = ordre d'affichage dans les colonnes "profile")
MLST_LOCI = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']
BAST_LOCI = ['fHbp_peptide', 'NHBA_peptide', 'NadA_peptide', 'PorA_VR1', 'PorA_VR2']
FINETYPING_LOCI = ['PorA_VR1', 'PorA_VR2', 'FetA_VR']

# Colonnes récapitulatives : nom de colonne -> loci concaténés (séparés par une virgule)
PROFILES = {
    'mlst_profile': MLST_LOCI,
    'bast_profile': BAST_LOCI,
    'finetyping_profile': FINETYPING_LOCI,
}

# Tous les loci (sans doublon, ordre conservé)
ALL_LOCI = list(dict.fromkeys(MLST_LOCI + BAST_LOCI + FINETYPING_LOCI))

# Champs "fields" renvoyés par PubMLST -> colonnes du rapport
FIELD_MAP = {
    'mlst': {'ST': 'st', 'clonal_complex': 'clonal_complex'},
    'bast': {
        'BAST': 'bast_type',
        'MenDeVAR_Bexsero_reactivity': 'bexsero_cross_reactivity',
        'MenDeVAR_Trumenba_reactivity': 'trumenba_cross_reactivity',
    },
}
BAST_COLUMNS = list(FIELD_MAP['bast'].values())

# Ordre des colonnes du CSV (les nouvelles colonnes sont à la fin pour ne rien décaler)
COLUMNS = (
    ['date', 'sample', 'st', 'clonal_complex']
    + MLST_LOCI + ['FetA_VR']
    + ['bast_type', 'fHbp_peptide', 'NHBA_peptide', 'NadA_peptide', 'PorA_VR1', 'PorA_VR2',
       'bexsero_cross_reactivity', 'trumenba_cross_reactivity']
    + list(PROFILES)
)

log = logging.getLogger('neisseria_typing')


# --------------------------------------------------------------------------- #
# Utilitaires
# --------------------------------------------------------------------------- #
def is_missing(value):
    return value in (None, '', EMPTY_VALUE)


def save_json(path, content):
    with open(path, 'w') as f:
        json.dump(content, f, indent=2)


def read_assembly_b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def post_json(url, payload):
    """POST JSON avec quelques essais en cas d'erreur réseau / 5xx.

    Une réponse 4xx (ex. 404 = pas de correspondance) est renvoyée telle quelle :
    ce n'est pas un crash, le corps JSON est exploitable.
    """
    body = json.dumps(payload).encode()
    last_error = None

    for attempt in range(1, API_RETRIES + 1):
        request = urllib.request.Request(
            url, data=body, method='POST', headers={'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(request, timeout=API_TIMEOUT) as response:
                return json.load(response)
        except urllib.error.HTTPError as err:
            if err.code < 500:
                return json.loads(err.read() or '{}')
            last_error = err
        except (urllib.error.URLError, TimeoutError) as err:
            last_error = err

        if attempt < API_RETRIES:
            time.sleep(2 * attempt)

    raise last_error


def sequence_url(kind, locus=None):
    if kind == 'locus':
        return f'{API_ROOT}/loci/{locus}/sequence'
    return f'{API_ROOT}/schemes/{SCHEME_ID[kind]}/sequence'


def query_pubmlst(url, sequence_b64, out_file):
    """Envoie l'assemblage à PubMLST et garde la réponse brute sur disque."""
    result = post_json(url, {'base64': True, 'sequence': sequence_b64})
    save_json(out_file, result)
    return result


# --------------------------------------------------------------------------- #
# Lecture des réponses
# --------------------------------------------------------------------------- #
def store_exact_matches(result, data):
    """Copie les allèles exacts dans `data`. Renvoie True s'il y en avait."""
    exact = result.get('exact_matches') or {}
    for locus, hits in exact.items():
        data[locus] = hits[0]['allele_id']
    return bool(exact)


def store_fields(result, kind, data):
    """Copie les champs du schéma (ST, CC, BAST, MenDeVAR...) dans `data`."""
    fields = result['fields']
    for db_key, column in FIELD_MAP[kind].items():
        data[column] = fields.get(db_key, EMPTY_VALUE)


def type_scheme(kind, sequence_b64, out_file, data):
    """Requête un schéma complet (mlst / bast / finetyping). True si identifié."""
    result = query_pubmlst(sequence_url(kind), sequence_b64, out_file)

    if not store_exact_matches(result, data):
        return False
    if 'fields' in result and kind in FIELD_MAP:
        store_fields(result, kind, data)
    return True


def confirm_bast(data, out_dir):
    """Redemande BAST + MenDeVAR à partir des allèles déjà trouvés (NadA doit valoir au moins "0")."""
    designations = {
        locus: [{'allele': data[locus]}] for locus in BAST_LOCI if not is_missing(data.get(locus))
    }
    result = post_json(f"{API_ROOT}/schemes/{SCHEME_ID['bast']}/designations",
                       {'designations': designations})
    save_json(os.path.join(out_dir, 'bast_type.json'), result)

    if 'fields' in result:
        store_fields(result, 'bast', data)
    else:
        for column in BAST_COLUMNS:
            data[column] = EMPTY_VALUE


def build_profile(data, loci):
    """'2,3,4,3,8,4,6' ; un locus manquant vaut '?' ; tout manquant -> '?'."""
    values = [EMPTY_VALUE if is_missing(data.get(locus)) else str(data[locus]) for locus in loci]
    if all(value == EMPTY_VALUE for value in values):
        return EMPTY_VALUE
    return ','.join(values)


# --------------------------------------------------------------------------- #
# Typage d'un échantillon
# --------------------------------------------------------------------------- #
def typing(assembly, sample, output_dir, output_filename, split):
    data = {'date': date.today().strftime('%d-%m-%Y'), 'sample': sample}
    best_matches = {}
    failed = set()
    sequence_b64 = read_assembly_b64(assembly)

    log.info(f'\n---------------------\n    {sample}\n---------------------')

    # 1) Les trois schémas
    for kind in ('mlst', 'bast', 'finetyping'):
        log.info(f'[Running] {kind}')
        out_file = os.path.join(output_dir, f'{sample}_{kind}.json')
        try:
            identified = type_scheme(kind, sequence_b64, out_file, data)
            log.info(f'✔️ {kind}' if identified else f'✖️ Cannot fetch {kind}')
        except Exception as err:
            failed.add(kind)
            log.error(f'Error in fetching {kind}: {err}')

    # 2) Loci manquants : on les redemande un par un
    missing_loci = [locus for locus in ALL_LOCI if is_missing(data.get(locus))]
    if missing_loci:
        log.info(f'Missing {missing_loci}...')

    for locus in missing_loci:
        out_file = os.path.join(output_dir, f'{sample}_{locus}.json')
        try:
            result = query_pubmlst(sequence_url('locus', locus), sequence_b64, out_file)
        except Exception as err:
            log.error(f'✖️ Cannot fetch the missing {locus}: {err}')
            continue

        if store_exact_matches(result, data):
            log.info(f'Fetched {locus}!')
        elif 'best_match' in result:
            best_matches[locus] = result['best_match'].get('allele_id')

    # 3) NadA absent = "0" (nécessaire pour obtenir le BAST)
    if is_missing(data.get('NadA_peptide')):
        data['NadA_peptide'] = '0'

    # 4) BAST / MenDeVAR parfois absents de la première réponse -> requête de confirmation
    if any(is_missing(data.get(column)) for column in BAST_COLUMNS):
        log.info('[Confirm BAST type and MenDeVar]')
        try:
            confirm_bast(data, output_dir)
        except Exception as err:
            log.error(f'✖️ Cannot confirm BAST: {err}')

        status = {'Bast Type': 'bast_type', 'Bexsero': 'bexsero_cross_reactivity',
                  'Trumenba': 'trumenba_cross_reactivity'}
        log.info(', '.join(f"{'✔️' if not is_missing(data.get(col)) else '✖️'} {label}"
                           for label, col in status.items()))

    # 5) Profil MLST complet mais pas de ST -> probablement une nouvelle souche
    if ('mlst' not in failed and is_missing(data.get('st'))
            and not any(is_missing(data.get(locus)) for locus in MLST_LOCI)):
        log.info('- [MLST] Full profile but no ST assignment; it must be new!')
        data['st'] = 'new'

    # 6) Colonnes récapitulatives, puis valeurs vides -> '?'
    for column, loci in PROFILES.items():
        data[column] = build_profile(data, loci)

    for column in COLUMNS:
        if is_missing(data.get(column)):
            data[column] = EMPTY_VALUE

    log.info(f'{sample}: {data}')

    # 7) Sorties
    row = {column: data[column] for column in COLUMNS}

    if split:  # un fichier .txt par colonne, lu par le WDL
        for column, value in row.items():
            with open(os.path.join(output_dir, f'{column}.txt'), 'w') as f:
                f.write(f'{value}\n')

    pd.DataFrame([row]).to_csv(os.path.join(output_dir, output_filename), index=False)
    pd.DataFrame([best_matches]).to_csv(os.path.join(output_dir, f'{sample}_best_matches.csv'), index=False)


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
    parser.add_argument('--split', action='store_true', help='Also write one .txt file per column')
    parser.add_argument('--files', help='Deprecated and ignored (kept so existing WDL calls still work)')
    return parser.parse_args()


def setup_logging(log_path):
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_path, mode='a', encoding='utf-8')],
    )


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f'Input not valid (expected one fasta file): {args.input}')

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(os.path.join(args.output_dir, args.log_file))

    sample = args.sample or os.path.basename(args.input).split('.')[0]
    typing(args.input, sample, args.output_dir, args.output_file, args.split)


if __name__ == '__main__':
    main()
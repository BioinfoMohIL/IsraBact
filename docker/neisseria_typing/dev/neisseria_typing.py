#!/usr/bin/env python3
import os
import json
from datetime import date
import argparse
import pandas as pd

EMPTY_VALUE        = '?'

COL_DATE           = 'date'
COL_SAMPLE         = 'sample'

KEY_ST             = 'st'
KEY_CC             = 'clonal_complex'
KEY_EXACT_MATCHES  = 'exact_matches'
KEY_BEST_MATCHES   = 'best_match'
KEY_ALLELE_ID      = 'allele_id'
KEY_FIELDS         = 'fields'
KEY_BAST           = "bast"
KEY_BAST_TYPE      = 'bast_type'
KEY_NADA_PEPTIDE   = 'NadA_peptide'
KEY_BEXERO         = 'bexsero_cross_reactivity'
KEY_TRUMENBA       = 'trumenba_cross_reactivity'

# keys name from pubmlt db output json
DB_NADA_PEPTIDE    = 'NadA_peptide'
DB_KEY_BEXERO      = "MenDeVAR_Bexsero_reactivity"
DB_KEY_TRUMENBA    = "MenDeVAR_Trumenba_reactivity"
DB_ST              = 'ST'
DB_CC              = 'clonal_complex'
DB_BAST_TYPE       = 'BAST'

TYPE_NADA          = 'nadA'
TYPE_LOCUS         = 'locus'
TYPE_MLST          = 'mlst'
TYPE_BAST          = 'bast'
TYPE_FINETYPING    = 'finetyping'

URL_LOCUS  = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci'
URL_SCHEME = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes'
URL_NADA   ="https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NadA_peptide/sequence"

NADA_PEPTIDE_COL = 'NadA_peptide'

BEXERO = "bexsero_cross_reactivity"
TRUMENBA = "trumenba_cross_reactivity"

SAMPLE_FIELD = 'sample'
SEROGROUP_FIELD = 'serogroup'
GENOGROUP_FIELD = 'genogroup'
SEROGROUP_RESULT_FILE_PATH = os.path.join('serogroup','serogroup_results.json')

db_id  = {
    'mlst': "1", 'bast': "53", 'finetyping': "2", 'nadA': '46',
    'B': '36', 'A': '29', 'C': '37', 'E': '35', 'H': '32', 'L': '33', 'W': '38',
    'X': '30', 'Y': '39', 'Z': '31'
}

alleles = {
    'bast':["fHbp_peptide", "NHBA_peptide", "NadA_peptide", "PorA_VR1", "PorA_VR2"],
    'mlst':['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm'] ,
    'finetyping': ['FetA_VR']
}

alleles_bast = ["fHbp_peptide", "NHBA_peptide", "NadA_peptide", "PorA_VR1", "PorA_VR2"]
alleles_mlst = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']

columns_ordered = [
    COL_DATE, COL_SAMPLE, KEY_ST, KEY_CC,
    'abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm','FetA_VR', 
    KEY_BAST_TYPE, "fHbp_peptide", "NHBA_peptide", KEY_NADA_PEPTIDE, 
    "PorA_VR1", "PorA_VR2", KEY_BEXERO, KEY_TRUMENBA
] 

wd      = os.path.dirname(os.path.realpath(__file__))
dev_dir = os.path.join(wd, 'dev')

def fetch_value_if_key(data, key):
    return data[key] if key in data else ''

def fetch_fields_values(data, type, dict):
    fields = data[KEY_FIELDS]

    if type == TYPE_MLST:
        dict[KEY_ST] = fetch_value_if_key(fields, DB_ST)
        dict[KEY_CC] = fetch_value_if_key(fields, DB_CC)

    elif type == TYPE_BAST:
        dict[KEY_BAST_TYPE] = fetch_value_if_key(fields, DB_BAST_TYPE)
        dict[KEY_BEXERO] = fetch_value_if_key(fields, DB_KEY_BEXERO)
        dict[KEY_TRUMENBA] = fetch_value_if_key(fields, DB_KEY_TRUMENBA)
    
# Fetch types alleles values
def fetch_allele_id(file, dict, type=None, locus=None, dict_best_matches=None):
    if os.path.isfile(file):
        f = open(file)
        data = json.load(f)
        f.close()

        is_exact_matches = KEY_EXACT_MATCHES in data and len(data[KEY_EXACT_MATCHES]) > 0

        if is_exact_matches:
            for i in data[KEY_EXACT_MATCHES]:
                dict[i] = data[KEY_EXACT_MATCHES][i][0][KEY_ALLELE_ID]
            
            if KEY_FIELDS in data:
                fetch_fields_values(data, type, dict)
            return True
        # If locus, that means the script has to fetch the value for each loci separatly
        # We can get exact matches or best matches
        elif type == TYPE_LOCUS:
            if is_exact_matches and len(data[KEY_EXACT_MATCHES]) > 0:
                for i in data[KEY_EXACT_MATCHES]:
                    dict[locus] = data[KEY_EXACT_MATCHES][i][0][KEY_ALLELE_ID]
            elif KEY_BEST_MATCHES in data:
                dict_best_matches[locus] = data[KEY_BEST_MATCHES][KEY_ALLELE_ID]
            return False
        else:
            return False

    else:
        return False

def fetch_bast(bast_alleles, dict, out_dir): 
    loci = populate_url(bast_alleles, dict)
    
    cmd = f'''
        curl -s -H "Content-Type: application/json" \
        -X POST "{URL_SCHEME}/{str(db_id['bast'])}/designations" \
        '''
    cmd+= "-d '{\"designations\": { " + loci + " }}' > " + out_dir + "/bast_type.json"
   
    os.system(cmd)

    return get_bast_from_data(out_dir + "/bast_type.json", dict)

def get_bast_from_data(file, dict): 
    if os.path.isfile(file):
        f = open(file)
        data = json.load(f) 
        f.close()
        
        if KEY_FIELDS in data:
            fetch_fields_values(data, TYPE_BAST, dict)
        else:
            dict[KEY_BAST_TYPE] = EMPTY_VALUE
            dict[KEY_BEXERO]    = EMPTY_VALUE
            dict[KEY_TRUMENBA]  = EMPTY_VALUE

def get_db_url(type, locus=None):
    if type == TYPE_LOCUS:
        return f"{URL_LOCUS}/{locus}/sequence"
    elif type == TYPE_NADA:
        return URL_NADA
    else:
        return f"{URL_SCHEME}/{db_id[type]}/sequence"

def get_profile(type, sequence, out_file, locus=None):
   
    cmd="(echo -n '{\"base64\":true,\"sequence\": \"';base64 "
    cmd+=f"\"{sequence}\""
    cmd+="; echo '\"}') | "
    cmd+=f"""
        curl -s \
        -H "Content-Type: application/json" \
        -X POST {get_db_url(type, locus)} \
        -d @-  | jq . > {out_file} 
    """

    try:
        os.system(cmd)
    except Exception as e:
        exit(f'Cannot get profile {type}\n{e}')
     
def populate_url(data, al):
    str = ''
    for val in data:
        if val in al:
            str += '"'+ val +'":[{"allele":"' + al[val] + '"}],'
    return str[:-1]

def check_if_locus_missing(data, alleles):
    missing_values = []
   
    for gene, locus in alleles.items():
        for loci in locus:
            if not loci in data:
                missing_values.append(loci)
    return missing_values

def get_value_back(type, locus, sequence, out_file):
    get_profile(type, sequence, out_file, locus)

def populate_file(type, file, content):
    with open(file, type) as f:
        f.write(f'{content}\n')

def populate_empty_values(keys, values=''):
    k = ''
    for key in keys:
        k += f",{key}"
        values += ','
    return k, values

def combine_keys(*elements):
    combined_list = []
    for element in elements:
        if isinstance(element, list):
            combined_list.extend(element)
        else:
            combined_list.append(element)
            
    return combined_list

def combine_type_keys(type):
    if type == TYPE_MLST:
        return combine_keys(KEY_ST, KEY_CC, alleles[TYPE_MLST])
    elif type == TYPE_BAST:
        return alleles[TYPE_BAST]
    elif type == TYPE_FINETYPING:
        return alleles[TYPE_FINETYPING]

def populate_empty_values(data):
    for k in columns_ordered:
        if k not in data:
            data[k] = EMPTY_VALUE

def create_report(df, file_name):
    df.to_csv(file_name, index=False)

def typing(assembly, output_dir, output_filename, files_name=None, log_file=None):
    data = {}
    best_matches = {}
    current_date = date.today().strftime("%d-%m-%Y")
    sample = os.path.basename(assembly).split('.')[0]
    
    data[COL_DATE]   = current_date
    data[COL_SAMPLE] = sample

    msg = f'''
    ---------------------
    
        {sample}
    
    ---------------------
    '''
    log_msg = msg
    print(msg)

    for type in [TYPE_MLST, TYPE_BAST, TYPE_FINETYPING]:
        print(f'\n[Running] : {type}')
        out_file = os.path.join(output_dir, f"{sample}_{type}.json")
        
        try:
            get_profile(type, assembly, out_file)

            is_identified = fetch_allele_id(out_file, data, type)
            
            if is_identified:
                msg = f'✔️ {type}'
                log_msg += msg
            else:
                msg = f"✖️ Cannot fetch {type}"
                log_msg += msg
            print(msg)

        except Exception as e:
            print(f'Error in fetching {type}\n{e}')
            
            keys = combine_type_keys(type)
            
            for k in keys:
                data[k] = EMPTY_VALUE

    # If locus missing, we will request the db locus by locus (one the missing ones)
    # On only the loci of MST, because we already double check bast and nada pedtide
    missing_locus = check_if_locus_missing(data, alleles)

    if missing_locus:
        print(f'\nMissing {missing_locus}...')
        for locus in missing_locus:
            
            out_file=os.path.join(output_dir, f"{sample}_{locus}.json")

            try:
                get_value_back(
                    type=TYPE_LOCUS,
                    locus=locus,
                    sequence=assembly,
                    out_file=out_file
                )
                msg = f'Fetched {locus}!'
            except Exception as e:
                msg = f'✖️ Cannot fetch the missing {locus}\n--> {e}'
                print(msg)
            
            fetch_allele_id(
                file=out_file, 
                dict=data, 
                dict_best_matches=best_matches, 
                type=TYPE_LOCUS, 
                locus=locus)

            log_msg += msg

    if (KEY_NADA_PEPTIDE not in data) or (data[KEY_NADA_PEPTIDE] == EMPTY_VALUE):
        data[KEY_NADA_PEPTIDE] = "0"
  
    # Sometime you cannot get bast type and MendeVar from the bast request (the 'fields' key
    # in the json file)
    # So we need to call a separate request
    # !We have to get a NadA Pedtide value (0 or something) to get Bast type, and MendeVar
    if any((key not in data) or (data[key] == EMPTY_VALUE) for key in [KEY_BAST_TYPE, KEY_BEXERO, KEY_TRUMENBA]):
        print('\n[Confirm BAST type and MenDeVar]')
        
        try:
            fetch_bast(alleles_bast, data, output_dir)
            
            msg = ''
        
            if data[KEY_BAST_TYPE] != EMPTY_VALUE:
                msg += f"✔️ Bast Type, "
            else:
                msg += "✖️ Bast Type, "
            
            if data[KEY_BEXERO] != EMPTY_VALUE:
                msg += f"✔️ Bexero, "
            else:
                msg += "✖️ Bexsero, "
            
            if data[KEY_TRUMENBA] != EMPTY_VALUE:
                msg += f"✔️ Trumenba "
            else:
                msg += "✖️ Trumenba"

        except Exception as e:
            msg += f"\n✖️ No confirm Bast \n{e}"
        
        log_msg += msg
        print(msg)

    os.system(f'echo "{sample}: {data}" >> {log_file}')
    
    # Set empty all value we couln't fetch
    populate_empty_values(data)
    
    df = pd.DataFrame([data])[columns_ordered]

    # If using WDL, need create txt files for each alleles to fetch the value 
    # at WDL output. So we need to create a txt file for each value
    if files_name:
        # Assume single row (e.g., from Neisseria typing results)
        data = df.iloc[0].to_dict()

        if args.split:
            for col in columns_ordered:
                with open(f"{output_dir}/{col}.txt", "w") as f:
                    f.write(f"{data.get(col, 'NA')}\n")

        
    create_report(df, f'{output_dir}/{output_filename}')
   
    if best_matches is not None:
        df_bm = pd.DataFrame([best_matches])
        create_report(df_bm, f'{output_dir}/{sample}_best_matches')

    if log_file:
        populate_file('a', log_file, log_msg)

   
parser = argparse.ArgumentParser()
parser.add_argument("--input", dest="input", required=True, help="Here your sample assembly (one fasta file or a fastas directory path)")
parser.add_argument("--output", dest="output_dir", required=True, help="Here your destination directory")
parser.add_argument("--output_csv_file", dest="output_file", required=True, help="Here your typing report csv filename")
parser.add_argument("--log-file", dest="log_file", default='logs.txt', required=False, help="Here your log file")
parser.add_argument('--files', help="Comma-separated keys", required=True)
parser.add_argument('--split', action='store_true', help="Output one file per key")

args = parser.parse_args()
output_dir = args.output_dir
my_input = args.input


if not os.path.isfile(my_input):
    exit('Input not valid (path folder of fastas or one fasta file)')

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

log_file = os.path.join(output_dir, args.log_file) if args.log_file else None

try:
    os.makedirs(output_dir, exist_ok=True) 
except Exception as e:
    print(f'Error in creating sample folder\n{e}')

files = args.files.split(',')

typing(
    assembly=my_input,
    output_dir=output_dir,
    log_file=log_file,
    output_filename=args.output_file,
    files_name=files
    )


#!/usr/bin/env python3
import os
import json
import argparse
import subprocess
import base64
import time
import random


URL_LOCUS  = 'https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci'
URL_SCHEME = 'https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes'
SCHEME = 'scheme'
LOCUS = 'locus'

loci = ["cyaA", "dnt"]

shemes = {
    "Autotransporters"      : '11',
    "Bp_vaccine_antigens"   : '7',
    "Macrolide_resistance"  : '6',
    "Phase"                 : '8',
    "PRN-test-Bp"           : '5',
    "T3SS"                  : '10'
}

all_alleles = [
    "23S_rRNA", "fhaBz-7180_10773", "fhaBy-3190_7183", "fhaBx-1_3193", "fhaB-2400_5550", "fhaB",
    "bipA", "bvgS", "bvgA",
    "cyaA", "dnt",
    "T3SStype", "bopN", "bopB", "bteA", "bsp22", "bopD",
    "brkA", "prn", "vag8", "bapC", "tcfA",
    "ptxC", "ptxA", "ptxP", "fim3", "ptxE", "ptxB", "ptxD", "BPagST"
]

def get_db_url(url, db_id):
    return f'{url}/{db_id}/sequence'

def get_profile(sequence_file, out, type, db_id, max_attempts=3, base_wait=5):
    """Fetch PubMLST profile JSON with retries and exponential backoff"""
    
    if type == SCHEME:
        url = URL_SCHEME
    elif type == LOCUS:
        url = URL_LOCUS
    else:
        raise ValueError("type must be SCHEME or LOCUS")

    # read and encode sequence
    with open(sequence_file, "rb") as f:
        b64_seq = base64.b64encode(f.read()).decode("utf-8")
    
    payload = json.dumps({"base64": True, "sequence": b64_seq})

    wait_time = base_wait * random.choice(range(10, 121, 20))

    for attempt in range(1, max_attempts+1):
        # optional random jitter (0–3 sec)
        time.sleep(random.uniform(0, 3))
        
        cmd = [
            "curl", "-s", "-w", "\n%{http_code}",
            "-H", "Content-Type: application/json",
            "-X", "POST", get_db_url(url, db_id),
            "-d", "@-"
        ]

        result = subprocess.run(cmd, input=payload, capture_output=True, text=True)
        lines = result.stdout.strip().split("\n")
        http_code = int(lines[-1])
        body = "\n".join(lines[:-1])

        if http_code == 200:
            # write JSON to output file
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with open(out, "w") as f:
                f.write(body)
            print(f"[✅] Success ({http_code}) -> {out}")
            return body
        elif http_code in (408, 429):
            print(f"[⚠️] Attempt {attempt}: HTTP {http_code}, retrying in {wait_time}s...")
            time.sleep(wait_time)
            wait_time *= 2  # exponential backoff
        else:
            print(f"[❌] HTTP {http_code}, body: {body}")
            break  # stop retries for other errors
    else:
        raise RuntimeError(f"Failed after {max_attempts} attempts, last HTTP code: {http_code}")

def typing(assembly, output_dir):
   
    for name, scheme_id in shemes.items():
        output_file = os.path.join(output_dir, f"{name}.json")
        # print(f"[🔎] Typing scheme {name} (ID: {scheme_id})...")
        get_profile(assembly, output_file, SCHEME, scheme_id) 


    for locus in loci:
        output_file = os.path.join(output_dir, f"{locus}.json")
        print(f"[🔎] Typing locus {locus} ...")
        get_profile(assembly, output_file, LOCUS, locus) 

def extract_alleles_ids(output_file, input_dir):
    """
    Extracts 'allele_id' values and additional 'fields' from JSON files in a given directory.
    Supports:
    - exact_matches as dict: gene -> list of allele dicts
    - exact_matches as list: anonymous list of alleles
    - fields as dict: additional info to extract (e.g., T3SStype)
    """
    seen_alleles = set() 
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_dir):
            if filename.endswith(".json"):
                scheme_name = filename.replace(".json", "")
                filepath = os.path.join(input_dir, filename)

                try:
                    with open(filepath, 'r') as f:
                        data = json.load(f)

                    # ---- Handle 'fields' if present ----
                    fields = data.get("fields", {})
                    if isinstance(fields, dict):
                        for key, value in fields.items():
                            print(f"{key}: {value}")
                            outfile.write(f"{key}: {value}\n")

                            seen_alleles.add(key)


                    # ---- Handle 'exact_matches' ----
                    exact_matches = data.get("exact_matches", {})

                    # Case 1: exact_matches is a dictionary (gene → list of allele matches)
                    if isinstance(exact_matches, dict):
                        for gene, alleles in exact_matches.items():
                            if isinstance(alleles, list) and alleles:
                                allele_id = alleles[0].get("allele_id", "NA")
                                print(f"{gene}: {allele_id}")
                                outfile.write(f"{gene}: {allele_id}\n")

                                seen_alleles.add(gene)


                    # Case 2: exact_matches is a list (no gene names provided)
                    elif isinstance(exact_matches, list):
                        for idx, allele_info in enumerate(exact_matches):
                            if isinstance(allele_info, dict):
                                allele_id = allele_info.get("allele_id", "NA")
                                
                                href = allele_info.get("href", "NA")
                                name = '<locus_name>'
                                if '/cyaA/' in href:
                                    name = 'cyaA'
                                elif '/dnt/' in href:
                                    name = 'dnt'
                                print(f"{name}: {allele_id}")
                                outfile.write(f"{name}: {allele_id}\n")
                    
                                seen_alleles.add(name)

                    # --- after processing all JSON files ---
                    for allele in all_alleles:
                        if allele not in seen_alleles:
                            print(f"{allele}: \"\"")
                            outfile.write(f"{allele}: \"\"\n")

                except Exception as e:
                    print(f"[❌] Failed to parse {filename}: {e}")     

                    
parser = argparse.ArgumentParser()
parser.add_argument("--assembly", dest="assembly", required=True, help="Here your sample assembly (one fasta file or a fastas directory path)")
parser.add_argument("--output_dir", dest="output_dir", default=f'results', required=False, help="Here your destination directory")
parser.add_argument("--output_ids_file", dest="output_file", default="alleles_ids.txt", required=False, help="Here your destination alleles ids file")


args = parser.parse_args()
output_dir = args.output_dir
alleles_ids_file = args.output_file
assembly = args.assembly


os.makedirs(output_dir, exist_ok=True)

typing(assembly, output_dir)

extract_alleles_ids(alleles_ids_file, output_dir)





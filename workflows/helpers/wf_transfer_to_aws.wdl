version 1.0

## ============================================================
## Workflow : TransferSampleToAWS
## Platform  : Terra.bio
## Description:
##   Prend un sample_id en input, localise ses reads (R1 obligatoire,
##   R2/R3/R4 optionnels) dans un TSV, les copie vers un bucket S3.
##   Si is_metadata=true (défaut), copie aussi la ligne metadata du
##   sample vers un second bucket S3 dédié aux métadonnées.
##
## TSV attendu (metadata_tsv) — colonnes obligatoires :
##   sample_id | read1 | read2 (opt) | read3 (opt) | read4 (opt) | ...
##   Les autres colonnes sont les métadonnées libres du sample.
##
## Utilisation :
##   1. Importer le WDL dans un workspace Terra
##   2. Remplir le JSON d'inputs (voir template)
##   3. Lancer avec dry_run=true pour valider, puis dry_run=false
## ============================================================

workflow wf_transfer_to_aws {

    input {
        # ── Sample ciblé ──────────────────────────────────────
        String  sample_id

        # ── TSV centralisé (dans le bucket GCS du workspace) ──
        # Colonnes obligatoires : sample_id, read1
        # Colonnes optionnelles : read2, read3, read4
        # Toutes les autres colonnes = métadonnées du sample
        File    metadata_tsv

        # ── Destinations AWS ──────────────────────────────────
        String  reads_bucket      # ex: s3://my-lab-reads/
        String  metadata_bucket   # ex: s3://my-lab-metadata/  (utilisé si is_metadata=true)

        # ── Credentials AWS ───────────────────────────────────
        String  aws_access_key_id
        String  aws_secret_access_key
        String  aws_region

        # ── Options ───────────────────────────────────────────
        Boolean is_metadata   = true     # copier les métadonnées du sample ?
        Boolean dry_run       = false    # simulation sans écriture réelle
        Boolean verify_md5    = true     # vérification intégrité post-transfer
        Int     max_retries   = 3
        String  storage_class = "INTELLIGENT_TIERING"

        # ── Resources ─────────────────────────────────────────
        Int     cpu       = 2
        Int     memory_gb = 8
        Int     disk_gb   = 100
    }

    # ── 1. Résoudre le sample : extraire les paths des reads + metadata ──
    # call ResolveSample {
    #     input:
    #         sample_id    = sample_id,
    #         metadata_tsv = metadata_tsv,
    #         cpu          = 1,
    #         memory_gb    = 4,
    #         disk_gb      = 10
    # }

    # # ── 2. Transférer les reads (R1 obligatoire, R2/R3/R4 optionnels) ──
    # call transfer_reads {
    #     input:
    #         sample_id             = sample_id,
    #         read1                 = ResolveSample.read1,
    #         read2                 = ResolveSample.read2,
    #         read3                 = ResolveSample.read3,
    #         read4                 = ResolveSample.read4,
    #         reads_bucket          = reads_bucket,
    #         aws_access_key_id     = aws_access_key_id,
    #         aws_secret_access_key = aws_secret_access_key,
    #         aws_region            = aws_region,
    #         dry_run               = dry_run,
    #         verify_md5            = verify_md5,
    #         max_retries           = max_retries,
    #         storage_class         = storage_class,
    #         cpu                   = cpu,
    #         memory_gb             = memory_gb,
    #         disk_gb               = disk_gb
    # }

    # # ── 3. (Optionnel) Transférer les métadonnées ─────────────────────
    # if (is_metadata) {
    #     call transfer_metadata {
    #         input:
    #             sample_id             = sample_id,
    #             metadata_row_json     = ResolveSample.metadata_row_json,
    #             metadata_bucket       = metadata_bucket,
    #             aws_access_key_id     = aws_access_key_id,
    #             aws_secret_access_key = aws_secret_access_key,
    #             aws_region            = aws_region,
    #             dry_run               = dry_run,
    #             storage_class         = storage_class,
    #             cpu                   = 1,
    #             memory_gb             = 4,
    #             disk_gb               = 10
    #     }
    # }

    output {
        File   resolve_log      = ResolveSample.resolve_log
        File   reads_log        = TransferReads.reads_log
        File?  metadata_log     = TransferMetadata.metadata_log
        File   md5_report       = TransferReads.md5_report
    }
}


## ─────────────────────────────────────────────────────────────
## TASK 1 — ResolveSample
##   Lit le TSV, localise le sample_id, valide read1,
##   détecte read2/read3/read4, extrait les métadonnées.
## ─────────────────────────────────────────────────────────────
task resolve_sample {
    input {
        String sample_id
        File   metadata_tsv
        Int    cpu
        Int    memory_gb
        Int    disk_gb
    }

    command <<<
        set -euo pipefail

        python3 << 'PYEOF'
        import csv, json, sys, os

        sample_id   = "~{sample_id}"
        tsv_path    = "~{metadata_tsv}"
        log_lines   = []

        log_lines.append(f"=== ResolveSample : {sample_id} ===")
        log_lines.append(f"TSV source : {tsv_path}")

        # ── Lire le TSV ──────────────────────────────────────────────
        with open(tsv_path, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames or []
            log_lines.append(f"Colonnes TSV : {', '.join(headers)}")

            # Vérifier colonnes obligatoires
            for col in ["sample_id", "read1"]:
                if col not in headers:
                    log_lines.append(f"ERREUR : Colonne obligatoire manquante → '{col}'")
                    print('\n'.join(log_lines))
                    with open("resolve_log.txt", "w") as lf:
                        lf.write('\n'.join(log_lines))
                    sys.exit(1)

            # Chercher le sample
            row = None
            for r in reader:
                if r.get("sample_id", "").strip() == sample_id:
                    row = r
                    break

        if row is None:
            log_lines.append(f"ERREUR : sample_id '{sample_id}' introuvable dans le TSV.")
            print('\n'.join(log_lines))
            with open("resolve_log.txt", "w") as lf:
                lf.write('\n'.join(log_lines))
            sys.exit(1)

        # ── Valider read1 ─────────────────────────────────────────────
        read1 = row.get("read1", "").strip()
        if not read1:
            log_lines.append(f"ERREUR : 'read1' est vide pour sample_id '{sample_id}'. "
                            "Au moins un fichier read est obligatoire.")
            print('\n'.join(log_lines))
            with open("resolve_log.txt", "w") as lf:
                lf.write('\n'.join(log_lines))
            sys.exit(1)

        log_lines.append(f"read1 détecté : {read1}")

        # ── Reads optionnels ─────────────────────────────────────────
        reads_found = {"read1": read1}
        for r in ["read2", "read3", "read4"]:
            val = row.get(r, "").strip()
            if val:
                reads_found[r] = val
                log_lines.append(f"{r} détecté : {val}")
            else:
                log_lines.append(f"{r} : absent (skipped)")

        # ── Écrire les paths dans des fichiers pour WDL ──────────────
        with open("read1.txt", "w") as f: f.write(reads_found["read1"])
        with open("read2.txt", "w") as f: f.write(reads_found.get("read2", ""))
        with open("read3.txt", "w") as f: f.write(reads_found.get("read3", ""))
        with open("read4.txt", "w") as f: f.write(reads_found.get("read4", ""))

        # ── Extraire toutes les métadonnées (hors colonnes read*) ─────
        meta_keys_to_skip = {"read1", "read2", "read3", "read4"}
        metadata_row = {k: v for k, v in row.items()
                        if k not in meta_keys_to_skip and v.strip()}

        log_lines.append(f"Métadonnées extraites : {list(metadata_row.keys())}")

        with open("metadata_row.json", "w") as f:
            json.dump(metadata_row, f, indent=2)

        log_lines.append("=== ResolveSample terminé — OK ===")
        print('\n'.join(log_lines))
        with open("resolve_log.txt", "w") as lf:
            lf.write('\n'.join(log_lines))
        PYEOF

    >>>

    output {
        File   resolve_log       = "resolve_log.txt"
        String read1             = read_string("read1.txt")
        String read2             = read_string("read2.txt")   # "" si absent
        String read3             = read_string("read3.txt")   # "" si absent
        String read4             = read_string("read4.txt")   # "" si absent
        File   metadata_row_json = "metadata_row.json"
    }

    runtime {
        docker:     "python:3.11-slim"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}


## ─────────────────────────────────────────────────────────────
## TASK 2 — TransferReads
##   Copie R1 (obligatoire) + R2/R3/R4 (si non vides) vers S3.
##   Échoue explicitement si R1 est vide (garde-fou).
## ─────────────────────────────────────────────────────────────
task transfer_reads {
    input {
        String  sample_id
        String  read1           # obligatoire
        String  read2           # "" si absent
        String  read3           # "" si absent
        String  read4           # "" si absent
        String  reads_bucket
        String  aws_access_key_id
        String  aws_secret_access_key
        String  aws_region
        Boolean dry_run
        Boolean verify_md5
        Int     max_retries
        String  storage_class
        Int     cpu
        Int     memory_gb
        Int     disk_gb
    }

    command <<<
        set -euo pipefail

        export AWS_ACCESS_KEY_ID="~{aws_access_key_id}"
        export AWS_SECRET_ACCESS_KEY="~{aws_secret_access_key}"
        export AWS_DEFAULT_REGION="~{aws_region}"

        SAMPLE="~{sample_id}"
        BUCKET="~{reads_bucket}"
        LOG="reads_transfer_${SAMPLE}.log"
        MD5_REPORT="md5_${SAMPLE}.txt"

        echo "=== TransferReads : $SAMPLE ===" | tee "$LOG"
        echo "Destination : $BUCKET"           | tee -a "$LOG"
        echo "Dry-run     : ~{dry_run}"        | tee -a "$LOG"
        echo "Date        : $(date -u)"        | tee -a "$LOG"
        echo ""                                | tee -a "$LOG"

        # Garde-fou : read1 ne doit jamais être vide ici
        if [ -z "~{read1}" ]; then
            echo "ERREUR FATALE : read1 est vide pour sample '~{sample_id}'." | tee -a "$LOG"
            echo "Aucun transfert effectué. Vérifiez le TSV." | tee -a "$LOG"
            exit 1
        fi

        # ── Construire la liste des reads à transférer ────────
        READS=()
        READS+=("~{read1}")
        [ -n "~{read2}" ] && READS+=("~{read2}")
        [ -n "~{read3}" ] && READS+=("~{read3}")
        [ -n "~{read4}" ] && READS+=("~{read4}")

        echo "Reads à transférer (${#READS[@]}) :" | tee -a "$LOG"
        for R in "${READS[@]}"; do echo "  - $R" | tee -a "$LOG"; done
        echo "" | tee -a "$LOG"

        # Initialiser le rapport MD5
        echo "sample_id,file,md5_local,etag_s3,status" > "$MD5_REPORT"

        # ── Transférer chaque read ────────────────────────────
        transfer_file() {
            local SRC="$1"
            local BASENAME
            BASENAME=$(basename "$SRC")
            local DEST="${BUCKET}${SAMPLE}/${BASENAME}"

            echo "--- Transfert : $BASENAME ---" | tee -a "$LOG"
            echo "  Source      : $SRC"          | tee -a "$LOG"
            echo "  Destination : $DEST"         | tee -a "$LOG"

            # MD5 local
            local MD5_LOCAL="N/A"
            if [ "~{verify_md5}" = "true" ]; then
                echo "  Calcul MD5..." | tee -a "$LOG"
                MD5_LOCAL=$(md5sum "$SRC" | awk '{print $1}')
                echo "  MD5 local : $MD5_LOCAL" | tee -a "$LOG"
            fi

            # Dry-run ou vrai transfert
            local STATUS="SUCCESS"
            local ETAG="N/A"

            if [ "~{dry_run}" = "true" ]; then
                echo "  [DRY-RUN] aws s3 cp $SRC $DEST --storage-class ~{storage_class}" | tee -a "$LOG"
                STATUS="DRY-RUN"
            else
                local ATTEMPT=0
                local OK=false
                while [ $ATTEMPT -lt ~{max_retries} ]; do
                    ATTEMPT=$((ATTEMPT + 1))
                    echo "  Tentative $ATTEMPT..." | tee -a "$LOG"
                    if aws s3 cp "$SRC" "$DEST" \
                            --storage-class "~{storage_class}" \
                            --no-progress 2>>"$LOG"; then
                        OK=true
                        break
                    fi
                    sleep $((ATTEMPT * 15))
                done

                if [ "$OK" = "false" ]; then
                    echo "  ERREUR : échec après ~{max_retries} tentatives pour $BASENAME" | tee -a "$LOG"
                    echo "~{sample_id},$BASENAME,$MD5_LOCAL,N/A,ERROR" >> "$MD5_REPORT"
                    return 1
                fi

                # ETag S3 post-transfer
                if [ "~{verify_md5}" = "true" ]; then
                    BUCKET_NAME=$(echo "$BUCKET" | sed 's|s3://||' | cut -d'/' -f1)
                    KEY=$(echo "$DEST" | sed "s|s3://$BUCKET_NAME/||")
                    ETAG=$(aws s3api head-object \
                        --bucket "$BUCKET_NAME" --key "$KEY" \
                        --query 'ETag' --output text | tr -d '"')
                    echo "  ETag S3   : $ETAG" | tee -a "$LOG"
                fi
                echo "  ✅ OK" | tee -a "$LOG"
            fi

            echo "~{sample_id},$BASENAME,$MD5_LOCAL,$ETAG,$STATUS" >> "$MD5_REPORT"
        }

        # Lancer le transfert pour chaque read
        ALL_OK=true
        for READ in "${READS[@]}"; do
            if ! transfer_file "$READ"; then
                ALL_OK=false
            fi
        done

        echo "" | tee -a "$LOG"
        if [ "$ALL_OK" = "false" ]; then
            echo "=== TERMINÉ AVEC ERREURS pour $SAMPLE ===" | tee -a "$LOG"
            exit 1
        else
            echo "=== TransferReads terminé — OK ($SAMPLE) ===" | tee -a "$LOG"
        fi
    >>>

    output {
        File reads_log  = "reads_transfer_~{sample_id}.log"
        File md5_report = "md5_~{sample_id}.txt"
    }

    runtime {
        docker:     "amazon/aws-cli:2.15.0"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}


## ─────────────────────────────────────────────────────────────
## TASK 3 — TransferMetadata  (appelée seulement si is_metadata=true)
##   Écrit un fichier JSON avec toutes les métadonnées du sample
##   et le pousse dans le bucket S3 de métadonnées.
##   Chemin S3 : <metadata_bucket>/<sample_id>/<sample_id>_metadata.json
## ─────────────────────────────────────────────────────────────
task transfer_metadata {
    input {
        String  sample_id
        File    metadata_row_json
        String  metadata_bucket
        String  aws_access_key_id
        String  aws_secret_access_key
        String  aws_region
        Boolean dry_run
        String  storage_class
        Int     cpu
        Int     memory_gb
        Int     disk_gb
    }

    command <<<
        set -euo pipefail

        export AWS_ACCESS_KEY_ID="~{aws_access_key_id}"
        export AWS_SECRET_ACCESS_KEY="~{aws_secret_access_key}"
        export AWS_DEFAULT_REGION="~{aws_region}"

        SAMPLE="~{sample_id}"
        BUCKET="~{metadata_bucket}"
        LOG="metadata_transfer_${SAMPLE}.log"
        META_FILE="${SAMPLE}_metadata.json"
        DEST="${BUCKET}${SAMPLE}/${META_FILE}"

        echo "=== TransferMetadata : $SAMPLE ===" | tee "$LOG"
        echo "Destination : $DEST"               | tee -a "$LOG"
        echo "Date        : $(date -u)"          | tee -a "$LOG"
        echo ""                                  | tee -a "$LOG"

        # ── Enrichir le JSON avec timestamp de migration ──────
        python3 - << 'PYEOF'
        import json, datetime

        with open("~{metadata_row_json}") as f:
            meta = json.load(f)

        meta["_migration_timestamp"] = datetime.datetime.utcnow().isoformat() + "Z"
        meta["_migration_dry_run"]   = "~{dry_run}"

        with open("~{sample_id}_metadata.json", "w") as f:
            json.dump(meta, f, indent=2)

        print(json.dumps(meta, indent=2))
        PYEOF

        echo "Contenu du fichier de métadonnées :" | tee -a "$LOG"
        cat "$META_FILE" | tee -a "$LOG"
        echo "" | tee -a "$LOG"

        # ── Dry-run ou transfert réel ─────────────────────────
        if [ "~{dry_run}" = "true" ]; then
            echo "[DRY-RUN] aws s3 cp $META_FILE $DEST --storage-class ~{storage_class}" | tee -a "$LOG"
            echo "STATUT : DRY-RUN — aucun fichier transféré." | tee -a "$LOG"
        else
            if aws s3 cp "$META_FILE" "$DEST" \
                    --storage-class "~{storage_class}" \
                    --no-progress 2>>"$LOG"; then
                echo "✅ Métadonnées transférées : $DEST" | tee -a "$LOG"
            else
                echo "ERREUR : Échec du transfert de métadonnées pour $SAMPLE" | tee -a "$LOG"
                exit 1
            fi
        fi

        echo "=== TransferMetadata terminé — OK ($SAMPLE) ===" | tee -a "$LOG"
    >>>

    output {
        File metadata_log = "metadata_transfer_~{sample_id}.log"
    }

    runtime {
        docker:     "amazon/aws-cli:2.15.0"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 2
    }
}

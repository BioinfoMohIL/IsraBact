# species_detection: everything used by the species detection workflows in ONE image
#   - Kraken2 2.14 (+ specie_detection.sh, as in the former species_detection:1 image)
#   - ML species classifier (predict_one.py + kmer_utils.py + strain_model.pkl)
#   - plate_check.py (plate check xlsx)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates tar g++ make rsync perl build-essential libz-dev procps \
        python3 python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- Kraken2 (same steps as species_detection:1) ---
WORKDIR /opt
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.14.tar.gz -O kraken2.tar.gz && \
    tar -xzf kraken2.tar.gz && mv kraken2-2.14 kraken2 && rm kraken2.tar.gz
WORKDIR /opt/kraken2
RUN chmod +x install_kraken2.sh && ./install_kraken2.sh /usr/local/bin
ENV PATH="/usr/local/bin:${PATH}"

# --- Python deps: pinned, the model (a pickle) must be read with the same scikit-learn version ---
RUN pip3 install --no-cache-dir numpy==1.26.4 pandas==2.2.3 scikit-learn==1.5.2 joblib==1.4.2

# --- ML + plate check (used by task_ml_species.wdl and task_plate_check.wdl) ---
COPY kmer_utils.py predict_one.py plate_check.py /opt/species_prediction/
COPY strain_model.pkl /opt/species_prediction/strain_model.pkl
COPY --from=bioinfomoh/specie_detection:1 /app/db /app/db
# Build-time guard: the build fails if the model was trained with another scikit-learn version
RUN python3 -c "import pickle, sklearn; b = pickle.load(open('/opt/species_prediction/strain_model.pkl', 'rb')); \
assert b['sklearn_version'] == sklearn.__version__, 'model trained with sklearn ' + b['sklearn_version'] + ', image: ' + sklearn.__version__; \
print('model OK: k =', b['k'], '| classes =', b['classes'])"

# --- Former species_detection:1 content ---
WORKDIR /app
COPY specie_detection.sh /app/specie_detection.sh
RUN chmod +x /app/specie_detection.sh

CMD ["bash"]
# ─────────────────────────────────────────────────────────────────────────────
# Dockerfile — draw_mst_network
# Génère un MST style GrapeTree (PNG) depuis un .nwk + métadonnées optionnelles
#
# Build :
#   docker build -t draw-mst-network:1.0 .
#
# Test local :
#   docker run --rm -v $(pwd):/data draw-mst-network:1.0 \
#       --input /data/tree.nwk --output /data/mst.png
#
# Push (adapter selon votre registry) :
#   docker tag draw-mst-network:1.0 your-dockerhub/draw-mst-network:1.0
#   docker push your-dockerhub/draw-mst-network:1.0
# ─────────────────────────────────────────────────────────────────────────────

FROM python:3.11-slim

LABEL maintainer="your-team@your-org.com" \
      description="MST network visualization from chewBBACA / GrapeTree output" \
      version="1.0"

# ── Système — fonts pour matplotlib ──────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        fontconfig \
        fonts-dejavu-core \
        libfreetype6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── Dépendances Python ────────────────────────────────────────────────────────
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# ── Script ───────────────────────────────────────────────────────────────────
COPY draw_mst_network.py /app/draw_mst_network.py

# Matplotlib : backend non-interactif (pas d'écran sur Terra/WDL)
ENV MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1

WORKDIR /data

ENTRYPOINT ["python", "/app/draw_mst_network.py"]
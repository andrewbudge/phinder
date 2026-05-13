#!/usr/bin/env bash
# phinder setup — create conda envs and download databases.
# Safe to re-run: skips anything already done.
#
# Usage:
#   bash setup.sh [--db-dir DIR] [--skip-envs] [--skip-dbs] [--with-phabox2]
#
# Options:
#   --db-dir DIR     Directory to store databases (default: ~/.phinder_dbs)
#   --skip-envs      Skip conda environment creation
#   --skip-dbs       Skip database downloads
#   --with-phabox2   Also install the phabox2 conda env (optional; requires git)

set -euo pipefail

# --- Defaults -----------------------------------------------------------------
DB_DIR="$HOME/.phinder_dbs"
SKIP_ENVS=false
SKIP_DBS=false
WITH_PHABOX2=false

# --- Argument parsing ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --db-dir)      DB_DIR="$2"; shift 2 ;;
        --skip-envs)   SKIP_ENVS=true; shift ;;
        --skip-dbs)    SKIP_DBS=true; shift ;;
        --with-phabox2) WITH_PHABOX2=true; shift ;;
        -h|--help)
            sed -n '2,11p' "$0" | sed 's/^# //'
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# --- Helpers ------------------------------------------------------------------
env_exists() { conda env list | grep -q "^${1} "; }

db_exists()  { [ -d "$1" ] && [ "$(ls -A "$1")" ]; }

# Detect downloader
if command -v curl &>/dev/null; then
    download() { curl -L -O "$1"; }
elif command -v wget &>/dev/null; then
    download() { wget "$1"; }
else
    echo "Error: neither curl nor wget found. Please install one and retry." >&2
    exit 1
fi

echo ""
echo "=== phinder setup ==="
echo "Database directory: $DB_DIR"
echo ""

# --- Conda environments -------------------------------------------------------
if [ "$SKIP_ENVS" = false ]; then
    echo "--- conda environments ---"

    if env_exists genomad_phinder; then
        echo "[envs] genomad_phinder already exists — skipping"
    else
        echo "[envs] creating genomad_phinder (geNomad 1.11.2)"
        conda create -y -n genomad_phinder -c bioconda -c conda-forge genomad=1.11.2
    fi

    if env_exists checkv_phinder; then
        echo "[envs] checkv_phinder already exists — skipping"
    else
        echo "[envs] creating checkv_phinder (CheckV 1.0.3)"
        conda create -y -n checkv_phinder -c bioconda -c conda-forge checkv=1.0.3
    fi

    if env_exists pharokka_phinder; then
        echo "[envs] pharokka_phinder already exists — skipping"
    else
        echo "[envs] creating pharokka_phinder (Pharokka 1.8.2)"
        conda create -y -n pharokka_phinder -c bioconda -c conda-forge pharokka=1.8.2
    fi

    if [ "$WITH_PHABOX2" = true ]; then
        if env_exists phabox2_phinder; then
            echo "[envs] phabox2_phinder already exists — skipping"
        else
            echo "[envs] creating phabox2_phinder (PhaBOX2 latest)"
            conda create -y -n phabox2_phinder -c conda-forge -c bioconda phabox=2.1.13
            conda run -n phabox2_phinder --no-capture-output conda install -y git
            conda run -n phabox2_phinder --no-capture-output bash -c "
                git clone https://github.com/KennthShang/PhaBOX.git /tmp/PhaBOX
                pip install /tmp/PhaBOX
                rm -rf /tmp/PhaBOX
            "
            echo "[envs] phabox2_phinder done"
        fi
    fi

    echo ""
fi

# --- Databases ----------------------------------------------------------------
if [ "$SKIP_DBS" = false ]; then
    mkdir -p "$DB_DIR"
    echo "--- databases ---"

    # geNomad
    if db_exists "$DB_DIR/genomad_db"; then
        echo "[dbs] genomad_db already exists — skipping"
    else
        echo "[dbs] downloading geNomad database"
        conda run -n genomad_phinder --no-capture-output \
            genomad download-database "$DB_DIR/genomad_db"
    fi

    # CheckV
    if db_exists "$DB_DIR/checkv_db"; then
        echo "[dbs] checkv_db already exists — skipping"
    else
        echo "[dbs] downloading CheckV database"
        conda run -n checkv_phinder --no-capture-output \
            checkv download_database "$DB_DIR/checkv_db"
    fi

    # Pharokka
    if db_exists "$DB_DIR/pharokka_db"; then
        echo "[dbs] pharokka_db already exists — skipping"
    else
        echo "[dbs] downloading Pharokka database"
        conda run -n pharokka_phinder --no-capture-output \
            install_databases.py -o "$DB_DIR/pharokka_db"
    fi

    # PhaBOX
    if db_exists "$DB_DIR/phabox_db_v2_2"; then
        echo "[dbs] phabox_db_v2_2 already exists — skipping"
    else
        echo "[dbs] downloading PhaBOX database"
        cd "$DB_DIR"
        download https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_2.zip
        unzip -q phabox_db_v2_2.zip
        rm phabox_db_v2_2.zip
        cd - > /dev/null
    fi

    echo ""
fi

# --- Summary ------------------------------------------------------------------
CONDA_BASE=$(conda info --base)

echo "=== setup complete ==="
echo ""
echo "Run phinder with:"
echo ""
echo "  nextflow run andrewcbudge/phinder \\"
echo "      --input contigs.fasta \\"
echo "      --genomad_db  ${DB_DIR}/genomad_db \\"
echo "      --checkv_db   ${DB_DIR}/checkv_db \\"
echo "      --pharokka_db ${DB_DIR}/pharokka_db \\"
echo "      --genomad_env  ${CONDA_BASE}/envs/genomad_phinder \\"
echo "      --checkv_env   ${CONDA_BASE}/envs/checkv_phinder \\"
echo "      --pharokka_env ${CONDA_BASE}/envs/pharokka_phinder \\"
echo "      -profile conda"
echo ""
if [ "$WITH_PHABOX2" = true ]; then
    echo "To enable PhaBOX, add:"
    echo "      --phabox2_env ${CONDA_BASE}/envs/phabox2_phinder \\"
    echo "      --phabox_db   ${DB_DIR}/phabox_db_v2_2"
else
    echo "To add PhaBOX (optional), rerun with --with-phabox2 or install manually and add:"
    echo "      --phabox2_env /path/to/phabox2/env \\"
    echo "      --phabox_db   ${DB_DIR}/phabox_db_v2_2"
fi
echo ""

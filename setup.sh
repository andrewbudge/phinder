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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$SCRIPT_DIR/envs"
DB_DIR="$HOME/.phinder_dbs"
SKIP_ENVS=false
SKIP_DBS=false
WITH_PHABOX2=false

# PhaBOX DB is the one database we fetch by pinned URL (its tool has no downloader).
# Single source of truth for both the download and the provenance record.
PHABOX_DB_VER="2.2"
PHABOX_DB_DIR="phabox_db_v2_2"
PHABOX_DB_URL="https://github.com/KennthShang/PhaBOX/releases/download/v2/${PHABOX_DB_DIR}.zip"

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

MANIFEST="$DB_DIR/DB_MANIFEST.tsv"

# --- Helpers ------------------------------------------------------------------
env_exists() { conda env list | grep -q "^${1} "; }

db_exists()  { [ -d "$1" ] && [ "$(ls -A "$1")" ]; }

# --- Provenance ---------------------------------------------------------------
# phinder doesn't pin DB versions — it records what each tool's downloader gave
# us, so a run can be described after the fact. Versions are read off disk, so
# the manifest reflects ground truth (fresh download or pre-existing install).

detect_genomad_ver() {
    if [ -f "$1/version.txt" ]; then tr -d '[:space:]' < "$1/version.txt"; else printf -- '-'; fi
}
detect_checkv_ver() {  # version lives in the checkv-db-v* subdir name
    local d
    d="$(find "$1" -maxdepth 1 -type d -name 'checkv-db-v*' -print -quit 2>/dev/null)" || true
    if [ -n "$d" ]; then basename "$d" | sed 's/^checkv-db-v//'; else printf -- '-'; fi
}
detect_pharokka_ver() {  # version lives in a VERSION_x_y_z marker filename
    local f
    f="$(find "$1" -maxdepth 1 -name 'VERSION_*' -print -quit 2>/dev/null)" || true
    if [ -n "$f" ]; then basename "$f" | sed 's/^VERSION_//; s/_/./g'; else printf -- '-'; fi
}

# tool version of the conda env that fetched a DB ("-" if that env isn't present)
tool_version() {  # $1 env, $2 pkg
    local v
    v="$(conda list -n "$1" "$2" 2>/dev/null | awk -v p="$2" '$1==p{x=$2} END{print x}')" || true
    printf '%s' "${v:--}"
}

init_manifest() {
    [ -f "$MANIFEST" ] && return
    mkdir -p "$DB_DIR"
    {
        echo "# phinder database manifest — records which databases this install fetched."
        echo "# phinder does not pin DB versions: it orchestrates tools and records what"
        echo "# their downloaders provided. To use a specific version, point the matching"
        echo "# --<tool>_db parameter at it. See the README \"Database provenance\" section."
        printf 'database\tdb_version\ttool_version\tsource\trecorded_utc\n'
    } > "$MANIFEST"
}

manifest_upsert() {  # $1 name, $2 db_version, $3 tool_version, $4 source
    init_manifest
    local now tmp
    now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    tmp="$(mktemp)"
    # keep comments, the header, and every other DB's row; replace this DB's row
    awk -F'\t' -v n="$1" '/^#/ || $1=="database" || $1!=n' "$MANIFEST" > "$tmp"
    printf '%s\t%s\t%s\t%s\t%s\n' "$1" "${2:--}" "${3:--}" "$4" "$now" >> "$tmp"
    mv "$tmp" "$MANIFEST"
}

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

    # Build each env from its pinned yml — the same files the pipeline uses,
    # so setup and run can never drift.
    if env_exists genomad_phinder; then
        echo "[envs] genomad_phinder already exists — skipping"
    else
        echo "[envs] creating genomad_phinder (from envs/genomad.yml)"
        conda env create -n genomad_phinder -f "$ENV_DIR/genomad.yml"
    fi

    if env_exists checkv_phinder; then
        echo "[envs] checkv_phinder already exists — skipping"
    else
        echo "[envs] creating checkv_phinder (from envs/checkv.yml)"
        conda env create -n checkv_phinder -f "$ENV_DIR/checkv.yml"
    fi

    if env_exists pharokka_phinder; then
        echo "[envs] pharokka_phinder already exists — skipping"
    else
        echo "[envs] creating pharokka_phinder (from envs/pharokka.yml)"
        conda env create -n pharokka_phinder -f "$ENV_DIR/pharokka.yml"
    fi

    if [ "$WITH_PHABOX2" = true ]; then
        if env_exists phabox2_phinder; then
            echo "[envs] phabox2_phinder already exists — skipping"
        else
            echo "[envs] creating phabox2_phinder (from envs/phabox2.yml)"
            conda env create -n phabox2_phinder -f "$ENV_DIR/phabox2.yml"
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
    if db_exists "$DB_DIR/$PHABOX_DB_DIR"; then
        echo "[dbs] $PHABOX_DB_DIR already exists — skipping"
    else
        echo "[dbs] downloading PhaBOX database"
        cd "$DB_DIR"
        download "$PHABOX_DB_URL"
        unzip -q "${PHABOX_DB_DIR}.zip"
        rm "${PHABOX_DB_DIR}.zip"
        cd - > /dev/null
    fi

    # --- Provenance manifest --------------------------------------------------
    # Record whatever landed on disk (freshly downloaded or pre-existing).
    if db_exists "$DB_DIR/genomad_db"; then
        manifest_upsert genomad  "$(detect_genomad_ver "$DB_DIR/genomad_db")" \
            "$(tool_version genomad_phinder genomad)"  "genomad download-database"
    fi
    if db_exists "$DB_DIR/checkv_db"; then
        manifest_upsert checkv   "$(detect_checkv_ver "$DB_DIR/checkv_db")" \
            "$(tool_version checkv_phinder checkv)"    "checkv download_database"
    fi
    if db_exists "$DB_DIR/pharokka_db"; then
        manifest_upsert pharokka "$(detect_pharokka_ver "$DB_DIR/pharokka_db")" \
            "$(tool_version pharokka_phinder pharokka)" "install_databases.py"
    fi
    if db_exists "$DB_DIR/$PHABOX_DB_DIR"; then
        manifest_upsert phabox   "$PHABOX_DB_VER" \
            "$(tool_version phabox2_phinder phabox2)"  "$PHABOX_DB_URL"
    fi
    echo "[dbs] provenance written to $MANIFEST"

    echo ""
fi

# --- Summary ------------------------------------------------------------------
CONDA_BASE=$(conda info --base)

echo "=== setup complete ==="
echo ""
if [ -f "$MANIFEST" ]; then
    echo "Databases recorded in: $MANIFEST"
    echo ""
fi
echo "Run phinder with:"
echo ""
echo "  nextflow run andrewbudge/phinder \\"
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

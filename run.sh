
#!/usr/bin/env sh

# ------------------------------------------
# Configuration
# ------------------------------------------
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
R_SCRIPT_NAME="main.R"
MAIN_R="$SCRIPT_DIR/$R_SCRIPT_NAME"
RSCRIPT_EXE=""

# ------------------------------------------
# Check that the R script exists
# ------------------------------------------
if [ ! -f "$MAIN_R" ]; then
    echo "Error: Could not find $R_SCRIPT_NAME in:"
    echo "$SCRIPT_DIR"
    echo
    exit 1
fi

# ------------------------------------------
# Function: compare two version strings
# returns 0 if $1 >= $2
# ------------------------------------------
version_ge() {
    [ "$1" = "$2" ] && return 0
    [ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | tail -n 1)" = "$1" ]
}

# ------------------------------------------
# Function: search a root directory for the
# newest R installation
# ------------------------------------------
find_latest_rscript() {
    SEARCH_ROOT="$1"

    [ -d "$SEARCH_ROOT" ] || return 1

    BEST_VERSION=""
    BEST_RSCRIPT=""

    for d in "$SEARCH_ROOT"/R-*; do
        [ -d "$d" ] || continue

        BASE="$(basename "$d")"
        VERSION="${BASE#R-}"

        if [ -x "$d/bin/Rscript" ]; then
            CANDIDATE="$d/bin/Rscript"
        elif [ -x "$d/bin/x64/Rscript" ]; then
            CANDIDATE="$d/bin/x64/Rscript"
        else
            continue
        fi

        if [ -z "$BEST_VERSION" ] || version_ge "$VERSION" "$BEST_VERSION"; then
            BEST_VERSION="$VERSION"
            BEST_RSCRIPT="$CANDIDATE"
        fi
    done

    if [ -n "$BEST_RSCRIPT" ]; then
        RSCRIPT_EXE="$BEST_RSCRIPT"
        return 0
    fi

    return 1
}

# ------------------------------------------
# Search common install locations first
# ------------------------------------------
find_latest_rscript "/Library/Frameworks/R.framework/Versions"
[ -n "$RSCRIPT_EXE" ] || find_latest_rscript "/usr/local/lib/R"
[ -n "$RSCRIPT_EXE" ] || find_latest_rscript "/opt/R"
[ -n "$RSCRIPT_EXE" ] || find_latest_rscript "/usr/lib/R"

# ------------------------------------------
# Fallback to PATH
# ------------------------------------------
if [ -z "$RSCRIPT_EXE" ]; then
    if command -v Rscript >/dev/null 2>&1; then
        RSCRIPT_EXE="$(command -v Rscript)"
    fi
fi

# ------------------------------------------
# If still not found, fail clearly
# ------------------------------------------
if [ -z "$RSCRIPT_EXE" ]; then
    echo "Error: Rscript was not found on this system."
    echo "Please install R and make sure Rscript is available."
    echo "See the README.md for setup instructions."
    echo
    exit 1
fi

# ------------------------------------------
# Run the R script
# ------------------------------------------
echo "=========================================="
echo "MCD UI Launcher"
echo "=========================================="
echo "Using R at:"
echo "$RSCRIPT_EXE"
echo
echo "Running script:"
echo "$MAIN_R"
echo "=========================================="
echo

"$RSCRIPT_EXE" --vanilla "$MAIN_R" "$@"
EXITCODE=$?

echo
if [ "$EXITCODE" -ne 0 ]; then
    echo "The script ended with an error."
else
    echo "The script finished successfully."
fi

exit "$EXITCODE"
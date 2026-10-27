#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Density spatial convergence test
#
# Fixed:
#   M = 2000
#
# Spatial grids:
#   N = 20, 40, 80, 160, 320
#
# Reference:
#   Nref = 1280
#
# Error norms:
#   L1, L2, Linf
#
# Order:
#   p = log(E_N / E_2N) / log(2)
# ============================================================

M=2000
NREF=1280
NS=(20 40 80 160 320)

OUTDIR="build/output"
RESULT="${OUTDIR}/rho_space_convergence_M${M}_Nref${NREF}.csv"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# Solver executable
# ------------------------------------------------------------
if [[ -x ./build/staggerdeGrid ]]; then
    SOLVER=./build/staggerdeGrid
elif [[ -x ./build/staggeredGrid ]]; then
    SOLVER=./build/staggeredGrid
else
    echo "Error: cannot find ./build/staggerdeGrid or ./build/staggeredGrid" >&2
    exit 1
fi


# ------------------------------------------------------------
# Output filename
#
# According to your C code:
# build/output/staggeredGrid_M%zu_N%zu_Rho.csv
# ------------------------------------------------------------
rho_file()
{
    local n=$1

    printf '%s/staggeredGrid_M%s_N%s_Rho.csv\n' \
        "$OUTDIR" "$M" "$n"
}


# ------------------------------------------------------------
# Check that rho output contains exactly N physical cells
# ------------------------------------------------------------
check_file()
{
    local file=$1
    local expected=$2

    if [[ ! -f "$file" ]]; then
        echo "Error: missing file: $file" >&2
        exit 1
    fi

    local nlines
    nlines=$(awk -F',' 'NF >= 2 {n++} END {print n+0}' "$file")

    if [[ "$nlines" -ne "$expected" ]]; then
        echo "Error: $file has $nlines rows, expected $expected." >&2
        echo "Density output should contain exactly N physical cells:" >&2
        echo "    j = 1,...,N" >&2
        echo "Do not output the ghost cell." >&2
        exit 1
    fi
}


# ============================================================
# 1. Compute reference solution
# ============================================================
echo "============================================================"
echo "Computing reference solution"
echo "M    = $M"
echo "Nref = $NREF"
echo "============================================================"

"$SOLVER" "$M" "$NREF"

REF_FILE=$(rho_file "$NREF")
check_file "$REF_FILE" "$NREF"


# Temporary file:
#
# N L1 L2 Linf
#
RAW_FILE=$(mktemp)
trap 'rm -f "$RAW_FILE"' EXIT


# ============================================================
# 2. Compute coarse-grid errors
# ============================================================
for N in "${NS[@]}"; do

    echo
    echo "============================================================"
    echo "Computing M=$M, N=$N"
    echo "============================================================"

    if (( NREF % N != 0 )); then
        echo "Error: NREF=$NREF is not divisible by N=$N" >&2
        exit 1
    fi

    RATIO=$((NREF / N))

    "$SOLVER" "$M" "$N"

    COARSE_FILE=$(rho_file "$N")
    check_file "$COARSE_FILE" "$N"


    # --------------------------------------------------------
    # Restrict fine reference solution to coarse cells.
    #
    # Coarse cell j corresponds to:
    #
    #   fine cell (j-1)*RATIO+1
    #       ...
    #   fine cell j*RATIO
    #
    # Reference density at coarse cell:
    #
    #   rho_ref(j)
    #     = average of corresponding fine-cell densities
    #
    # CSV first column (Euler position) is deliberately ignored.
    # --------------------------------------------------------
    read -r L1 L2 LINF < <(

        awk -F',' \
            -v ratio="$RATIO" \
            -v nc_expected="$N" \
            -v nr_expected="$NREF" '

        # First input file: fine reference solution
        FNR == NR {
            if (NF >= 2) {
                nr++
                ref[nr] = $2 + 0.0
            }
            next
        }

        # Second input file: coarse solution
        NF >= 2 {

            nc++

            rho_c = $2 + 0.0

            first = (nc - 1) * ratio + 1
            last  = nc * ratio

            sum_ref = 0.0

            for (k = first; k <= last; k++) {
                sum_ref += ref[k]
            }

            rho_ref = sum_ref / ratio

            err = rho_c - rho_ref

            if (err >= 0.0) {
                aerr = err
            } else {
                aerr = -err
            }

            sum1 += aerr
            sum2 += err * err

            if (aerr > maxerr) {
                maxerr = aerr
            }
        }

        END {

            if (nr != nr_expected) {
                printf "Reference row count mismatch: %d != %d\n", nr, nr_expected > "/dev/stderr"
                exit 2
            }

            if (nc != nc_expected) {
                printf "Coarse row count mismatch: %d != %d\n", nc, nc_expected > "/dev/stderr"
                exit 2
            }

            # Domain length X = 1
            h = 1.0 / nc

            L1   = h * sum1
            L2   = sqrt(h * sum2)
            Linf = maxerr

            printf "%.17e %.17e %.17e\n", L1, L2, Linf
        }

        ' "$REF_FILE" "$COARSE_FILE"
    )


    printf '%d %.17e %.17e %.17e\n' \
        "$N" "$L1" "$L2" "$LINF" >> "$RAW_FILE"

done


# ============================================================
# 3. Compute convergence orders
# ============================================================
{
    echo "N,L1,order_L1,L2,order_L2,Linf,order_Linf"

    awk '

    {
        N    = $1
        L1   = $2
        L2   = $3
        Linf = $4

        if (NR == 1) {

            printf "%d,%.16e,,%.16e,,%.16e,\n", N, L1, L2, Linf

        } else {

            p1 = log(prev1 / L1)   / log(2.0)
            p2 = log(prev2 / L2)   / log(2.0)
            pi = log(previ / Linf) / log(2.0)

            printf "%d,%.16e,%.8f,%.16e,%.8f,%.16e,%.8f\n", N, L1, p1, L2, p2, Linf, pi
        }

        prev1 = L1
        prev2 = L2
        previ = Linf
    }

    ' "$RAW_FILE"

} > "$RESULT"


# ============================================================
# 4. Print table
# ============================================================
echo
echo "============================================================"
echo "Density spatial convergence"
echo "M = $M, reference N = $NREF"
echo "============================================================"

printf '%-8s %-16s %-10s %-16s %-10s %-16s %-10s\n' \
    "N" "L1" "order" "L2" "order" "Linf" "order"

printf '%-8s %-16s %-10s %-16s %-10s %-16s %-10s\n' \
    "--------" "----------------" "----------" \
    "----------------" "----------" \
    "----------------" "----------"

awk -F',' '
NR > 1 {
    printf "%-8s %-16s %-10s %-16s %-10s %-16s %-10s\n", \
           $1, $2, $3, $4, $5, $6, $7
}
' "$RESULT"


echo
echo "Result saved to:"
echo "    $RESULT"

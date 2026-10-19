#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")/.."
mkdir -p build/tests

# Test only functions changed for the indexing migration.
# GNU C retains the original M_PI definition without editing scalar routines.
flags=(-std=gnu11 -Wall -Wextra -Wpedantic -Werror -g -O1
       -fsanitize=undefined,bounds -fno-sanitize-recover=all
       -fno-omit-frame-pointer)

gcc "${flags[@]}" tests/mesh_calculation_test.c mesh_calculation.c \
    linear_algebra_2d.c -lm -o build/tests/mesh_calculation_test
gcc "${flags[@]}" tests/initial_boundary_test.c initial_bound_condition.c \
    linear_algebra_2d.c -lm \
    -o build/tests/initial_boundary_test
gcc "${flags[@]}" tests/data_output_test.c data_output.c \
    -lm -o build/tests/data_output_test

for test_name in mesh_calculation_test initial_boundary_test data_output_test; do
    printf 'Running %s\n' "$test_name"
    "./build/tests/$test_name"
done
printf 'All three focused test executables passed.\n'

#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")/.."
mkdir -p build/tests

# Optional first argument: a GCC sanitizer list, or "none".
# Optional second argument: a platform-specific sanitizer options source.
sanitizers="${1:-undefined,bounds}"
flags=(-std=c11 -Wall -Wextra -Wpedantic -Werror -g -O1
       -fno-omit-frame-pointer)
if [[ "$sanitizers" != "none" ]]; then
    flags+=("-fsanitize=$sanitizers" -fno-sanitize-recover=all)
fi
support_sources=()
if [[ -n "${2:-}" ]]; then
    support_sources+=("$2")
fi

gcc "${flags[@]}" tests/linearAlgebra_test.c linear_algebra_2d.c \
    "${support_sources[@]}" -lm -o build/tests/linearAlgebra_test
gcc "${flags[@]}" tests/mesh_calculation_test.c mesh_calculation.c \
    "${support_sources[@]}" -lm -o build/tests/mesh_calculation_test
gcc "${flags[@]}" tests/initial_boundary_test.c initial_bound_condition.c \
    linear_algebra_2d.c "${support_sources[@]}" -lm \
    -o build/tests/initial_boundary_test
gcc "${flags[@]}" tests/data_output_test.c data_output.c \
    "${support_sources[@]}" -lm -o build/tests/data_output_test
# solver_test includes the production main source with its entry point renamed.
gcc "${flags[@]}" tests/solver_test.c linear_algebra_2d.c mesh_calculation.c \
    physical_calculation.c initial_bound_condition.c data_output.c \
    "${support_sources[@]}" -lm -o build/tests/solver_test

for test_name in linearAlgebra_test mesh_calculation_test \
                 initial_boundary_test data_output_test solver_test; do
    printf 'Running %s\n' "$test_name"
    "./build/tests/$test_name"
done
printf 'All five test executables passed.\n'

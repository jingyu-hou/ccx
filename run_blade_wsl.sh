#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/ccx/test
mkdir -p /mnt/d/ccx/test_logs
log="/mnt/d/ccx/test_logs/blade_$(date +%Y%m%d_%H%M%S).log"
solver="/mnt/d/ccx/src/WeICME_MT"

{
  echo "== Blade test run =="
  date
  echo "Solver: $solver"
  echo "Input:  /mnt/d/ccx/test/blade.inp"

  if [ ! -x "$solver" ]; then
    echo "ERROR: solver is missing or not executable: $solver"
    exit 1
  fi
  if [ ! -f blade.inp ]; then
    echo "ERROR: input file not found: /mnt/d/ccx/test/blade.inp"
    exit 1
  fi

  echo
  echo "== Runtime environment =="
  : "${OMP_NUM_THREADS:=2}"
  : "${CCX_NPROC_EQUATION_SOLVER:=$OMP_NUM_THREADS}"
  export OMP_NUM_THREADS
  export CCX_NPROC_EQUATION_SOLVER
  echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
  echo "CCX_NPROC_EQUATION_SOLVER=$CCX_NPROC_EQUATION_SOLVER"
  uname -a

  echo
  echo "== Start solver =="
  "$solver" blade
  status=$?

  echo
  echo "== Solver exit status =="
  echo "$status"

  echo
  echo "== Output files =="
  ls -lh blade.* || true

  echo
  echo "== Tail blade.sta =="
  tail -n 80 blade.sta 2>/dev/null || true

  echo
  echo "== Tail blade.dat =="
  tail -n 80 blade.dat 2>/dev/null || true

  exit "$status"
} 2>&1 | tee "$log"

echo "Log written to: $log"
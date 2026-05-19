#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/ccx
mkdir -p build_logs
log="build_logs/build_$(date +%Y%m%d_%H%M%S).log"

{
  echo "== System =="
  if command -v lsb_release >/dev/null 2>&1; then
    lsb_release -a
  else
    cat /etc/os-release
  fi
  uname -a

  echo
  echo "== Toolchain check =="
  missing=0
  for tool in gcc gfortran make ar ranlib perl; do
    if command -v "$tool" >/dev/null 2>&1; then
      echo "-- $tool"
      "$tool" --version 2>&1 | head -n 1 || true
    else
      echo "MISSING: $tool"
      missing=1
    fi
  done

  if ! command -v cmake >/dev/null 2>&1; then
    echo "MISSING: cmake (not required by current Makefile, but recommended)"
  else
    echo "-- cmake"
    cmake --version 2>&1 | head -n 1 || true
  fi

  if [ "$missing" -ne 0 ]; then
    echo
    echo "Installing required build tools with apt..."
    sudo apt update
    sudo apt install -y build-essential gfortran make perl cmake libblas-dev liblapack-dev libarpack2-dev
  fi

  echo
  echo "== Project files =="
  ls -la
  test -f Makefile
  test -d src
  test -d ARPACK
  test -d SPOOLES.2.2

  echo
  echo "== Build =="
  make

  echo
  echo "== Build outputs =="
  find /mnt/d/ccx -maxdepth 3 -type f \( -name 'WeICME*' -o -name 'ccx*' -o -name '*.exe' \) -printf '%p\n' 2>/dev/null || true

  echo
  echo "BUILD_FINISHED"
} 2>&1 | tee "$log"

echo "Log written to: /mnt/d/ccx/$log"
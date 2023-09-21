#!/bin/bash
set -e

cleanup=0

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    -l|--local-install)
      install_prefix="--prefix $(pwd)/external"
      echo "Local install enabled"
      shift 1
      ;;
    -c|--cleanup)
      cleanup=1
      shift 1
      ;;
  esac
done

# In case you want CMake to use a different compiler than the system default, use
#
#    export CC=path/to/your-c-compiler
#    export FC=path/to/your-fortran-compiler

mkdir -p external
cd external

if [ -d "arpack-ng" ] && [ -n "$(ls -A arpack-ng)" ]; then
  cd arpack-ng
  git pull
else
  git clone https://github.com/opencollab/arpack-ng.git
  cd arpack-ng
fi

cmake -B build
cmake --build build --parallel
cmake --install build $install_prefix

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/arpack-ng/
fi

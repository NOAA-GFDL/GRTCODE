#!/bin/bash

# Get the directory where the environment script lives.
workflow_home="/ncrc/home1/Jing.Feng/scripts/grtcode/GRTworkflow" #"$( dirname -- "$( readlink -f -- "$0"; )"; )"

# Parse the command line arguments to see if the user wants to use a custom version of grtcode.
argument_list="$0 [-h|--help] [-p path to grtcode repository]"
default_grtcode_repository="$( readlink -f -- $workflow_home/../grtcode )"
grtcode_repository="$default_grtcode_repository"
while [[ $# -gt 0 ]]; do
  argument="$1"
  case $argument in
    -h|--help)
      printf "Usage: $argument_list\n"
      printf "\nOptional arguments:"
      printf "-h, --help:    Print this message.\n"
      printf "-p <path>:     Path to a grtcode repository.\n"
      exit 0
    ;;
    -p)
      shift
      if [ "$#" -gt 0 ]; then
        grtcode_repository="$1"
        shift
      else
        printf "Error: please specify the path to the grtcode repository.\n"
        printf "Usage: $argument_list\n"
        exit 1
      fi
    ;;
    *)
      printf "Error: Unknown argument.\n"
      printf "Usage: $argument_list\n"
      exit 1
    ;;
  esac
done
grtcode_repository="$( readlink -f -- $grtcode_repository )"
commit="$( git -C $grtcode_repository rev-parse HEAD )"
remote_url="$( git -C $grtcode_repository config --get remote.origin.url )"
local_changes="$( git -C $grtcode_repository status --porcelain=v1 2>/dev/null | wc -l )"
printf "Using grtcode repository: $grtcode_repository\n"
printf "  [$remote_url @ $commit]\n"
if [ "$local_changes" -gt 0 ]; then
  printf "Warning: $grtcode_repository has local changes.\n"
fi

# If you are not using the basic installation, then build your custom version.
  # Source the environment script.
  source /ncrc/home1/Jing.Feng/scripts/grtcode/GRTworkflow/environment.sh

  # Prepare to build the libraries and executables.
  current_location="$PWD"
  cd $grtcode_repository
  mkdir -p build

  # Clean and then re-build the executables.
  makefile="$grtcode_repository/Makefile"
  make -f $makefile clean
  CC="${CC}" \
    CPPFLAGS="${CPPFLAGS}" \
    CFLAGS="${CFLAGS} -diag-disable=10441" \
    CXX="${CXX}" \
    CXXFLAGS="${CXXFLAGS} -diag-disable=10441" \
    LDFLAGS="${LDFLAGS}" \
    MPICC="${MPICC}" \
    MPICXX="${MPICXX}" \
    make -f $makefile all rfmip-irf era5
  printf "Built executable in $grtcode_repository/build\n"

# Check which executables exist.
printf "Checking executables:\n"
for executable in circ rfmip-irf era5; do
  if [ -x "$grtcode_repository/build/$executable" ]; then
    printf "\t$grtcode_repository/build/$executable ... created successfully.\n"
  else
    printf "\t$grtcode_repository/build/$executable ... failed.\n"
  fi
done

if [ "$USER" == "Raymond.Menzel" ]; then
  # Build the combiner.
  $CC -o $workflow_home/grtcode-results-combiner -g -O3 -qopenmp -diag-disable=10441 \
    -I$grtcode_repository/utilities/src \
    $workflow_home/grtcode-results-combiner.c \
    $grtcode_repository/build/argparse.o -lnetcdf
fi

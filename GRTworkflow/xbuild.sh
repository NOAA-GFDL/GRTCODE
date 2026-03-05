# Get the directory where the environment script lives.
workflow_home="$( dirname -- "$( readlink -f -- "$0"; )"; )"

# Source the environment script.
source $workflow_home/environment.sh

# Build the combiner.
$CC -o combine -g -O0 -qopenmp -diag-disable=10441 \
    -I$workflow_home/../grtcode/utilities/src \
    $workflow_home/combine-v2.c \
    $workflow_home/../grtcode/build/argparse.o -lnetcdf

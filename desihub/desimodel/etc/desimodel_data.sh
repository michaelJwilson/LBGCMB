#!/bin/bash -x
set -e
#
# If this script is being run by Travis, the Travis install script will
# already be in the correct directory.
# If this script is being run by desiInstall, then we need to make sure
# we are running this in ${WORKING_DIR}.
#
[[ -n "${WORKING_DIR}" ]] && cd ${WORKING_DIR}
#
# Make sure DESIMODEL_VERSION is set.
#
if [[ -z "${DESIMODEL_VERSION}" ]]; then
    echo "DESIMODEL_VERSION is not set!"
    exit 1
fi
svn export https://desi.lbl.gov/svn/code/desimodel/${DESIMODEL_VERSION}/data
#
# Set this for subsequent Travis tests.  For desiInstall, this environment
# variable should already be set when the desimodel Module file is
# processed.
#
[[ -z "${DESIMODEL}" ]] && export DESIMODEL=$(pwd)

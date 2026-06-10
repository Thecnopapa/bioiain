#!/bin/bash

VIRTUAL_ENV="./venv"
PATH="./venv/bin:$PATH"
which python
pwd
ls


if [[ -z "$TARGET_VERSION" ]]; then
  echo "Updating bioiain to the newest version"
  pip3 install bioiain -U

else
  echo "Installing specific bioiain version: ${TARGET_VERSION}"
  pip3 install bioiain==${TARGET_VERSION}
fi

#echo "$(pip3 show bioiain)"
#echo "$(pip3 show bioiain | grep Version)"


BIOIAIN_VERSION=$(pip3 show bioiain | grep Version)
#echo $BIOIAIN_VERSION
BIOIAIN_VERSION=($BIOIAIN_VERSION)
#echo ${BIOIAIN_VERSION[1]}
BIOIAIN_VERSION=${BIOIAIN_VERSION[1]}
echo "BIOIAIN_VERSION=${BIOIAIN_VERSION}"

OUTPUT_FOLDER="/docs/${BIOIAIN_VERSION}"
mkdir -p OUTPUT_FOLDER

echo "OUTPUT_FOLDER=${OUTPUT_FOLDER}"

pdoc --html -o OUTPUT_FOLDER -f --skip-errors bioiain
#!/bin/sh

VIRTUAL_ENV="./venv"
PATH="./venv/bin:$PATH"
which python
pwd
ls

echo "$(pip3 show bioiain)"
echo "$(pip3 show bioiain | grep Version)"

BIOIAIN_VERSION=$(pip3 show bioiain | grep Version)
echo $BIOIAIN_VERSION
BIOIAIN_VERSION=($BIOIAIN_VERSION)
echo $BIOIAIN_VERSION
BIOIAIN_VERSION=${BIOIAIN_VERSION[1]}
echo "BIOIAIN_VERSION=${BIOIAIN_VERSION}"

OUTPUT_FOLDER="/docs/${BIOIAIN_VERSION}"
mkdirs OUTPUT_FOLDER

echo "OUTPUT_FOLDER=${OUTPUT_FOLDER}"

python pdoc --html -o OUTPUT_FOLDER --skip-errors bioiain
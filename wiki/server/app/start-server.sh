#!/bin/bash

VIRTUAL_ENV="./venv"
PATH="./venv/bin:$PATH"
which python
pwd
ls

echo "Starting server..."

python ./server.py
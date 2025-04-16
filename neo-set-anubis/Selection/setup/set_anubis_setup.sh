#!/bin/bash

# check conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda is not installed. Please install Anaconda or Miniconda first."
    exit
fi

# create env from the yaml
conda env create --file set_anubis_setup.yml

# activate the env
echo "To activate this environment, use:"
echo "conda activate <your_env_name>"

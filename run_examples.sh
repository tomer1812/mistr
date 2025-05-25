#!/bin/bash

echo "Starting MISTR and MISTR-IV examples execution..."

echo "Running src/example_data_generation.R..."
Rscript src/example_data_generation.R

echo "Running src/example_qrist.py..."
python src/example_qrist.py

echo "Running src/example_estimate_hte_mistr.R..."
Rscript src/example_estimate_hte_mistr.R

echo "Running src/example_estimate_hte_mistr-iv.R..."
Rscript src/example_estimate_hte_mistr-iv.R
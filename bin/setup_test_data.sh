#!/bin/bash

# Script to setup test data for the pipeline
# This ensures test data exists regardless of where the pipeline is run from

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PIPELINE_DIR/data"

echo "Setting up test data in: $DATA_DIR"

# Create data directory if it doesn't exist
mkdir -p "$DATA_DIR"

# Check if test data already exists
if [[ -f "$DATA_DIR/genome.fasta" && -f "$DATA_DIR/nanopore_reads.fastq.gz" ]]; then
    echo "Test data already exists."
    exit 0
fi

# Generate test data
echo "Generating test data..."
cd "$PIPELINE_DIR"

# Run the data generation scripts
if [[ -f "generate_example_data.py" ]]; then
    python generate_example_data.py
fi

if [[ -f "data/generate_nanopore_reads.py" ]]; then
    cd data && python generate_nanopore_reads.py
    cd "$PIPELINE_DIR"
fi

echo "Test data setup complete!"
echo "You can now run the pipeline with: nextflow run . -profile test"

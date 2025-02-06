#!/bin/bash

# Help message
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Process ND2 files for NFISH analysis"
    echo ""
    echo "Options:"
    echo "  -d, --directory DIR    Directory containing ND2 files (required)"
    echo "  -o, --output NAME      Output folder name (default: nfish_results)"
    echo "  -t, --threads NUM      Number of threads to use (default: 8)"
    echo "  -e, --env NAME         Conda environment name (optional)"
    echo "  -h, --help            Show this help message"
}

# Default values
THREADS=8
OUTPUT_DIR="nfish_results"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--directory)
            ND2_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -e|--env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check if directory is provided
if [ -z "$ND2_DIR" ]; then
    echo "Error: Directory not specified"
    show_help
    exit 1
fi

# Activate conda environment if specified
if [ ! -z "$CONDA_ENV" ]; then
    echo "Activating conda environment: $CONDA_ENV"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate "$CONDA_ENV"
fi

# Run the processor
echo "Processing ND2 files in: $ND2_DIR"
echo "Output folder: $OUTPUT_DIR"
echo "Using $THREADS threads"
python nfish_processor.py "$ND2_DIR" --output "$OUTPUT_DIR" --threads "$THREADS"
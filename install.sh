#!/bin/bash
# Installation script for genetique-analysis package

set -e

echo "Installing genetique-analysis package..."

# Check if uv is available
if command -v uv &> /dev/null; then
    echo "Using uv for installation..."
    uv pip install -e .
elif command -v pip &> /dev/null; then
    echo "Using pip for installation..."
    pip install -e .
else
    echo "Error: Neither uv nor pip found. Please install one of them."
    exit 1
fi

echo "Installation completed successfully!"
echo ""
echo "You can now use the package with:"
echo "  genetique-analysis your_project_name"
echo "  genetique-error-tests [project_name] [selection_name] [pop_name] [nb_inds] [nb_iterations]"
echo ""
echo "Or import it in Python:"
echo "  from genetique_analysis import FileConfiguration"
echo "  from genetique_analysis.cli import main_error_recaptures"
echo ""
echo "See README.md for more information."

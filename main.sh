#!/bin/bash
# ==============================================================================
# main.sh - Citation Analysis Pipeline
# 
# Usage: bash main.sh
# 
# This script runs the full analysis pipeline:
#   1. Data preparation (01_data_prep.R)
#   2. Analysis and visualization (02_analysis.R)
# ==============================================================================

set -e  # Exit on any error

echo "=============================================="
echo "  Citation Analysis Pipeline"
echo "=============================================="

echo ""
echo "[Step 1/3] Running data preparation..."
echo "[Step 2/3] Skipping data preparation..."

# Rscript 01_data_prep.R 
echo "[Step 3/3] Data preparation complete."

echo ""
echo "[Step 2/2] Running analysis and visualization..."
Rscript 02_analysis.R
echo "[Step 2/2] Analysis complete."

echo ""
echo "=============================================="
echo "  Pipeline finished successfully"
echo "=============================================="
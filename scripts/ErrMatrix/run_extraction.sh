#!/bin/bash

# Script to extract testing positions from both adaptive and fixed methods
# Run this on the server to extract the data for bug analysis

echo "Extracting testing positions for bug analysis..."
echo "Testing positions: 19780000, 19790000, 19800000, 19810000, 19820000, 19830000, 19840000"
echo ""

# Run the extraction script
Rscript extract_positions_server.R

echo ""
echo "Extraction complete!"
echo "Files created:"
echo "- testing_positions_comparison.rds"
echo ""
echo "You can now download this file to your local machine for analysis."

#!/bin/bash
set -e

echo "Installing dtcc-core with build dependencies..."

# Install dtcc-pyspade-native first (build dependency)
echo "Step 1/2: Installing dtcc-pyspade-native..."
pip install "dtcc-pyspade-native@git+https://github.com/dtcc-platform/dtcc-pyspade-native.git@main"

# Install dtcc-core
echo "Step 2/2: Installing dtcc-core..."
pip install -e .

echo "Installation complete!"

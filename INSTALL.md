# BEER Python Implementation - Installation Guide

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements_python.txt

# 2. Run example
python example_usage.py

# 3. Run tests
python test_beer.py
```

## Detailed Installation

### Option 1: Using pip (Recommended)

```bash
# Create virtual environment (recommended)
python -m venv beer_env
source beer_env/bin/activate  # On Windows: beer_env\Scripts\activate

# Install core dependencies
pip install numpy scipy pandas scikit-learn anndata scanpy

# Install optional dependencies
pip install combat-python bbknn matplotlib seaborn

# Install development dependencies (optional)
pip install pytest black flake8 mypy
```

### Option 2: Using conda

```bash
# Create conda environment
conda create -n beer python=3.9
conda activate beer

# Install from conda-forge
conda install -c conda-forge numpy scipy pandas scikit-learn anndata scanpy

# Install pip-only packages
pip install combat-python bbknn

# Install visualization
conda install -c conda-forge matplotlib seaborn
```

### Option 3: Using requirements file

```bash
# Install all dependencies at once
pip install -r requirements_python.txt
```

## Verifying Installation

### Check Python Version

```bash
python --version
# Should be Python 3.8 or higher
```

### Test Imports

```python
# test_imports.py
import numpy as np
import pandas as pd
import scipy
import sklearn
import anndata as ad
import scanpy as sc

print("✓ Core dependencies installed successfully")

# Test optional dependencies
try:
    from combat.pycombat import pycombat
    print("✓ combat-python installed")
except ImportError:
    print("✗ combat-python not installed (optional)")

try:
    import bbknn
    print("✓ bbknn installed")
except ImportError:
    print("✗ bbknn not installed (optional)")

try:
    import matplotlib.pyplot as plt
    print("✓ matplotlib installed")
except ImportError:
    print("✗ matplotlib not installed (optional)")

# Test BEER import
try:
    from beer import BEER
    print("✓ BEER module can be imported")
except ImportError as e:
    print(f"✗ BEER module import failed: {e}")
```

Run the test:
```bash
python test_imports.py
```

### Run Unit Tests

```bash
# Run all tests
python test_beer.py

# Expected output:
# ================================================================================
# BEER Python Implementation - Test Suite
# ================================================================================
#
# test_aggregate_by_group_mean (test_beer.TestUtilityFunctions) ... ok
# test_aggregate_by_group_sum (test_beer.TestUtilityFunctions) ... ok
# ...
# ================================================================================
# Test Summary
# ================================================================================
# Tests run: 30+
# Successes: 30+
# Failures: 0
# Errors: 0
# ================================================================================
```

## Troubleshooting

### Issue: ImportError for combat-python

```
ImportError: No module named 'combat'
```

**Solution:**
```bash
pip install combat-python
```

### Issue: ImportError for bbknn

```
ImportError: No module named 'bbknn'
```

**Solution:**
```bash
pip install bbknn
```

Note: BBKNN is optional. BEER will work without it, but BBKNN enhancement will be unavailable.

### Issue: NumPy/SciPy installation fails

```
ERROR: Failed building wheel for numpy
```

**Solution (Linux/Mac):**
```bash
# Install system dependencies first
sudo apt-get install python3-dev  # Ubuntu/Debian
# or
brew install python  # macOS

# Then install NumPy/SciPy
pip install numpy scipy
```

**Solution (Windows):**
```bash
# Use pre-built wheels
pip install --upgrade pip
pip install numpy scipy --only-binary :all:
```

### Issue: Scanpy/AnnData installation fails

**Solution:**
```bash
# Try installing from conda-forge
conda install -c conda-forge scanpy

# Or use pip with pre-built wheels
pip install scanpy --only-binary :all:
```

### Issue: Memory errors during testing

```
MemoryError: Unable to allocate array
```

**Solution:**
- Close other applications
- Reduce test data size
- Use a machine with more RAM (recommended: 8 GB+)

### Issue: sklearn import fails

```
ImportError: No module named 'sklearn'
```

**Solution:**
```bash
pip install scikit-learn
```

Note: Package name is `scikit-learn` but import name is `sklearn`.

## Platform-Specific Notes

### Linux (Ubuntu/Debian)

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install python3-pip python3-dev python3-venv
sudo apt-get install build-essential

# Create virtual environment
python3 -m venv beer_env
source beer_env/bin/activate

# Install BEER
pip install -r requirements_python.txt
```

### macOS

```bash
# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python
brew install python

# Create virtual environment
python3 -m venv beer_env
source beer_env/bin/activate

# Install BEER
pip install -r requirements_python.txt
```

### Windows

```bash
# Install Python from python.org (if not already installed)
# https://www.python.org/downloads/

# Create virtual environment
python -m venv beer_env
beer_env\Scripts\activate

# Install BEER
pip install -r requirements_python.txt
```

## Docker Installation (Advanced)

### Dockerfile

```dockerfile
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements
COPY requirements_python.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements_python.txt

# Copy BEER files
COPY beer.py .
COPY example_usage.py .
COPY test_beer.py .

# Default command
CMD ["python", "test_beer.py"]
```

### Build and Run

```bash
# Build Docker image
docker build -t beer-python .

# Run tests
docker run beer-python

# Run interactive shell
docker run -it beer-python bash

# Run with mounted data directory
docker run -v /path/to/data:/data beer-python python your_script.py
```

## Upgrading

### Upgrade BEER

```bash
# Pull latest code from GitHub
git pull origin master

# No additional installation needed (pure Python)
```

### Upgrade Dependencies

```bash
# Upgrade all dependencies
pip install --upgrade -r requirements_python.txt

# Upgrade specific package
pip install --upgrade scanpy
```

## Uninstallation

### Remove Virtual Environment

```bash
# Deactivate environment
deactivate

# Remove environment directory
rm -rf beer_env
```

### Remove System-wide Installation

```bash
# Uninstall dependencies
pip uninstall numpy scipy pandas scikit-learn anndata scanpy combat-python bbknn matplotlib seaborn
```

## Getting Help

### Check Version Information

```python
import numpy as np
import scipy
import pandas as pd
import scanpy as sc

print(f"NumPy: {np.__version__}")
print(f"SciPy: {scipy.__version__}")
print(f"pandas: {pd.__version__}")
print(f"scanpy: {sc.__version__}")
```

### Report Issues

If you encounter installation issues:

1. Check this guide for common solutions
2. Verify Python version (should be 3.8+)
3. Check dependency versions match requirements
4. Open an issue on GitHub with:
   - Operating system and version
   - Python version
   - Error message (full traceback)
   - Output of `pip list`

### Community Support

- **GitHub Issues**: https://github.com/jumphone/BEER/issues
- **Original BEER**: https://github.com/jumphone/BEER
- **Scanpy Help**: https://discourse.scverse.org/

## Next Steps

After successful installation:

1. **Run Examples**: `python example_usage.py`
2. **Read Documentation**: See [README_PYTHON.md](README_PYTHON.md)
3. **Explore API**: Read docstrings in [beer.py](beer.py)
4. **Try Your Data**: Follow examples to analyze your own data

---

**Last Updated**: October 2025
**Python Version**: 3.8+
**Platform**: Linux, macOS, Windows

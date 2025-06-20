# Contributing to SHALSTAB

Thank you for your interest in contributing to the SHALSTAB (Shallow Landsliding STABility) project! This document provides guidelines and information for contributors.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Code Standards](#code-standards)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)
- [Documentation](#documentation)

## Getting Started

### Prerequisites

- Python 3.8 or higher
- Git
- Basic understanding of slope stability analysis and geospatial data processing

### Required Dependencies

```bash
pip install -r requirements.txt
```

## Development Setup

1. **Fork and Clone the Repository**

   ```bash
   git clone https://github.com/federicogmz/shalstab.git
   cd shalstab
   ```

2. **Create a Virtual Environment**

   ```bash
   python -m venv shalstab_env
   source shalstab_env/bin/activate
   # On Windows:
   shalstab_env\Scripts\activate
   ```

3. **Install Dependencies**

   ```bash
   pip install -r requirements.txt
   pip install -e . # Install in development mode
   ```

4. **Verify Installation**
   ```bash
   python -c "from shalstab import ShalstabAnalyzer; print('Installation successful!')"
   ```

## Code Standards

### Python Style Guide

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use meaningful variable and function names
- Maximum line length: 88 characters (Black formatter standard)
- Use type hints for function parameters and return values

### Docstring Format

We use NumPy-style docstrings. All public methods and classes must include:

```python
def example_method(self, param1: float, param2: str) -> xr.DataArray:
    """
    Brief description of the method.

    Longer description if needed, explaining the method's purpose,
    algorithm, or important implementation details.

    Parameters
    ----------
    param1 : float
        Description of param1.
    param2 : str
        Description of param2.

    Returns
    -------
    xarray.DataArray
        Description of the returned data.

    Raises
    ------
    ValueError
        When invalid parameters are provided.
    """
```

### Code Organization

- Keep methods focused and single-purpose
- Use private methods (prefixed with `_`) for internal functionality
- Group related methods with clear section headers
- Maintain consistent error handling patterns

## Submitting Changes

### Workflow

1. **Create a Feature Branch**

   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Your Changes**

   - Write clean, well-documented code
   - Add tests for new functionality
   - Update documentation as needed

3. **Run Quality Checks**

   ```bash
   # Format code
   black shalstab/ tests/

   # Check style
   flake8 shalstab/ tests/

   # Run tests
   pytest tests/
   ```

4. **Commit Your Changes**

   ```bash
   git add .
   git commit -m "Add feature: brief description of changes"
   ```

5. **Push and Create Pull Request**
   ```bash
   git push origin feature/your-feature-name
   ```

### Pull Request Guidelines

- Provide a clear title and description
- Reference any related issues: `Fixes #123`
- Include screenshots for UI changes
- Ensure all tests pass
- Request review from maintainers

### Commit Message Format

Use clear, descriptive commit messages:

```
Add feature: slope stability visualization options

- Add color scheme customization for stability plots
- Include legend positioning options
- Update documentation with visualization examples

Fixes #45
```

## Reporting Issues

### Bug Reports

When reporting bugs, please include:

- Python version and operating system
- SHALSTAB version
- Minimal code example that reproduces the issue
- Expected vs. actual behavior
- Error messages and stack traces
- Sample data files (if applicable)

### Feature Requests

For feature requests, please describe:

- The use case and motivation
- Proposed implementation approach
- Potential impact on existing functionality
- Alternative solutions considered

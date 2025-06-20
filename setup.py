"""Setup configuration for SHALSTAB package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

# Read requirements
requirements = [
    "rasterio>=1.3.0",
    "rioxarray>=0.13.0",
    "numpy>=1.21.0",
    "xarray>=2022.6.0",
    "geopandas>=0.12.0",
    "scipy>=1.9.0",
    "pysheds>=0.3.0",
    "matplotlib>=3.5.0",
    "geocube>=0.3.0",
]

# Optional dependencies for development
dev_requirements = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=22.0.0",
    "flake8>=5.0.0",
    "sphinx>=5.0.0",
    "sphinx-rtd-theme>=1.0.0",
    "jupyter>=1.0.0",
    "pre-commit>=2.20.0",
]

setup(
    name="shalstab",
    version="1.0.0",
    author="Federico Gómez",
    author_email="fjgomezc@eafit.edu.co",
    description="Shallow Landsliding STABility (SHALSTAB) model for slope stability analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/federicogmz/shalstab",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": dev_requirements,
        "all": dev_requirements,
    },
    keywords=[
        "shalstab",
        "slope stability",
        "landslide",
        "geohazard",
        "hydrology",
        "geotechnical",
        "gis",
        "raster",
    ],
    project_urls={
        "Bug Reports": "https://github.com/federicogmz/shalstab/issues",
        "Source": "https://github.com/federicogmz/shalstab",
        "Documentation": "https://github.com/federicogmz/shalstab/blob/main/README.md",
    },
    include_package_data=True,
    zip_safe=False,
)

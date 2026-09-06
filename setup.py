from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="chemcompute",
    version="0.1.1",
    author="Mohammad Keifari",
    author_email="mohammadkeifari2007@gmail.com",
    description="A Python library for chemical reaction simulation, including kinetic modeling and thermodynamic equilibrium calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MohammadKeifari/ChemCompute",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.19.0",
        "matplotlib>=3.3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.0.0",
        ],
    },
    keywords="chemistry, chemical reactions, kinetics, thermodynamics, equilibrium, simulation",
    project_urls={
        "Bug Reports": "https://github.com/MohammadKeifari/ChemCompute/issues",
        "Source": "https://github.com/MohammadKeifari/ChemCompute",
        "Documentation": "https://github.com/MohammadKeifari/ChemCompute#readme",
    },
)


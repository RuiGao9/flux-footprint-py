from setuptools import setup, find_packages

setup(
    name="flux-footprint-py",
    version="0.1.0",
    packages=find_packages(),
    # automatically include required libraries when installing the package
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "openpyxl",   
        "rasterio",   
    ],
    author="Rui Gao, Mohammad Safeeq, and Joshua H. Viers",
    description="A Python tool for Kljun flux footprint modeling for research purposes.",
)
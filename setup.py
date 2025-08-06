#!/usr/bin/env python3
"""
Setup script for Monte Carlo Radiation Transport Simulator
Build system using pybind11 for C++/Python integration
"""

from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools.command.build_ext import build_ext as _build_ext
import pybind11
import sys
import os

# Project metadata
__version__ = "0.1.0"
__author__ = "Lucas Holik"
__email__ = "lucas.holik88@gmail.com"

# Compiler arguments for different platforms
def get_compiler_flags():
    """Get platform-specific compiler flags"""
    flags = []
    
    if sys.platform.startswith('win'):
        # Windows-specific flags
        flags.extend(['/O2', '/std:c++17'])
    else:
        # Unix-like systems (Linux, macOS)
        flags.extend(['-O3', '-std=c++17', '-ffast-math'])
        
        # Add OpenMP support if available
        if sys.platform.startswith('linux'):
            flags.extend(['-fopenmp'])
    
    return flags

def get_linker_flags():
    """Get platform-specific linker flags"""
    flags = []
    
    if sys.platform.startswith('linux'):
        flags.extend(['-fopenmp'])
    elif sys.platform.startswith('darwin'):  # macOS
        # Note: May need special handling for OpenMP on macOS
        pass
        
    return flags

# Define the pybind11 extension module
ext_modules = [
    Pybind11Extension(
        "mcshield",  # Module name that will be importable in Python
        [
            # C++ source files - start with a simple hello world
             "cpp/tests/test_python_binding.cpp",
        ],
        include_dirs=[
            # Include directories
            "cpp/include",
            pybind11.get_cmake_dir() + "/../../../include",
        ],
        cxx_std=17,  # C++17 standard
        extra_compile_args=get_compiler_flags(),
        extra_link_args=get_linker_flags(),
        define_macros=[
            ('VERSION_INFO', f'"{__version__}"'),
        ],
    ),
]

# Custom build_ext class for better error handling
class CustomBuildExt(build_ext):
    """Custom build extension with better error messages"""
    
    def build_extensions(self):
        try:
            super().build_extensions()
        except Exception as e:
            print(f"Build failed: {e}")
            print("\nTroubleshooting tips:")
            print("1. Ensure you have a C++ compiler installed")
            print("2. Check that pybind11 is properly installed: pip install pybind11")
            print("3. Verify your C++ source files exist in cpp/src/")
            raise

# Setup configuration
setup(
    name="monte-carlo-shielding",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="Monte Carlo Radiation Transport Simulator for Nuclear Shielding Applications",
    long_description="""
A professional-grade Monte Carlo radiation transport simulator built for nuclear 
shielding calculations. Features multi-particle transport (photons, neutrons, 
electrons) with realistic physics and high-performance C++ computational engine.
    """.strip(),
    long_description_content_type="text/plain",
    url="https://github.com/your-username/monte-carlo-shielding",  # Update when you create repo
    
    # Python package configuration
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    
    # C++ extension modules
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt},
    
    # Dependencies
    install_requires=[
        "numpy>=1.19.0",
        "matplotlib>=3.3.0",
        "streamlit>=1.0.0",
        "plotly>=5.0.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
    ],
    
    # Development dependencies
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "test": [
            "pytest>=6.0",
            "pytest-cov",
        ],
    },
    
    # Python version requirement
    python_requires=">=3.8",
    
    # Classification
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: C++",
    ],
    
    # Include additional files
    include_package_data=True,
    zip_safe=False,
    
    # Entry points for command-line tools (optional)
    entry_points={
        "console_scripts": [
            "mc-shield=mcshield.cli:main",  # We can add this later
        ],
    },
)

# Post-installation message
if __name__ == "__main__":
    print("\n" + "="*60)
    print("Monte Carlo Radiation Transport Simulator")
    print("="*60)
    print(f"Version: {__version__}")
    print(f"Author: {__author__}")
    print("\nBuild configuration:")
    print(f"- C++ Standard: C++17")
    print(f"- Platform: {sys.platform}")
    print(f"- Python: {sys.version.split()[0]}")
    print(f"- Compiler flags: {' '.join(get_compiler_flags())}")
    
    if get_linker_flags():
        print(f"- Linker flags: {' '.join(get_linker_flags())}")
    
    print("\nNext steps:")
    print("1. Create C++ source files in cpp/src/")
    print("2. Run: python setup.py build")
    print("3. Run: pip install -e .")
    print("="*60)
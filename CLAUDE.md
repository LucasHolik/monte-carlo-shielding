# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a professional Monte Carlo radiation transport simulator built for nuclear shielding applications. The project features a high-performance C++ computational engine with Python bindings and a Streamlit web interface for interactive simulations. It models multi-particle transport (photons, neutrons, electrons) through complex geometries with realistic physics interactions.

## Build System and Commands

### Core Build Commands

```bash
# Build C++ extensions and install package
python setup.py build
pip install -e .

# Alternative modern build using pyproject.toml
pip install -e .

# Install dependencies
pip install -r requirements.txt

# Install development dependencies
pip install -e .[dev]
```

### Testing Commands

```bash
# Run Python tests
python -m pytest tests/
pytest --cov=mcshield tests/  # with coverage

# Compile and test C++ components
g++ -std=c++17 -I cpp/include cpp/tests/test_compile.cpp -o test_compile
./test_compile

# Note: User will run .exe files manually and provide output to Claude
# Claude should create .exe files but not execute them directly
```

### Code Quality Commands

```bash
# Format code
black python/
flake8 python/

# Type checking
mypy python/
```

### Running the Application

```bash
# Launch Streamlit web interface
streamlit run python/web-app/app.py

# Use Python API directly
python -c "import mcshield; print('Import successful')"
```

## Architecture Overview

### Core C++ Engine (`cpp/`)

The simulation engine is built in modern C++17 with these key components:

- **Transport System** (`Transport.hpp/cpp`): Main particle transport engine that handles particle movement through geometries, boundary crossings, and interaction sampling
- **Monte Carlo Sampling** (`MonteCarloSampling.hpp/cpp`): Comprehensive sampling utilities for interaction distances, directions, energies, and physics-specific sampling (Compton scattering, photoelectric effect)
- **Geometry System** (`Geometry.hpp`, `Box.hpp`): Spatial geometry handling with ray-tracing capabilities for complex 3D structures
- **Particle System** (`Particle.hpp/cpp`): Particle data structures and state management
- **Materials** (`Material.hpp/cpp`): Material properties and composition handling
- **Random Number Generation** (`RandomNumberGenerator.hpp/cpp`): High-quality random number generation with validation (`RandomNumberQuality.hpp`)
- **Debugging Tools** (`SamplingDebugger.hpp`): Statistical validation and debugging utilities for Monte Carlo sampling

### Python Integration (`python/`)

- **mcshield/**: Main Python package with pybind11 bindings to C++ engine
- **visualization/**: Plotting and analysis tools using matplotlib/plotly
- **web-app/**: Streamlit interface for interactive simulations

### Key Design Patterns

- **Modern C++17**: Extensive use of smart pointers, STL containers, and RAII
- **Functional Programming**: Cross-section functions as std::function objects for flexible physics models
- **Template-based Design**: Generic algorithms for different particle types
- **Comprehensive Error Handling**: Validation throughout the transport pipeline

## Development Guidelines

### C++ Code Structure

- Header files in `cpp/include/` with corresponding implementations in `cpp/src/`
- Comprehensive test files in `cpp/tests/` following naming convention `test_*.cpp`
- Use modern C++17 features: smart pointers, auto, range-based loops
- Physics calculations follow established nuclear data libraries (NIST, ENDF)
- **Include Path Convention**: Use simple header file names in #include statements (e.g., `#include "Material.hpp"` not `#include "../include/Material.hpp"` or `#include "cpp/include/Material.hpp"`)

### Python Integration

- All C++ classes exposed via pybind11 in the `mcshield` module
- NumPy array integration for efficient data transfer
- Follow Python naming conventions in the Python layer while maintaining C++ naming in the core

### Physics Implementation

- **Photon Transport**: Photoelectric absorption, Compton scattering, pair production, Rayleigh scattering
- **Neutron Transport**: Elastic scattering, absorption reactions, thermal neutron physics
- **Cross-sections**: Energy-dependent with interpolation capabilities
- **Geometry**: Ray-tracing algorithms for complex 3D geometries

### Testing Strategy

- C++ unit tests for individual components
- Python integration tests via pytest
- Statistical validation of Monte Carlo sampling
- Physics validation against analytical solutions

## Important Implementation Notes

- The project follows a 5-week development timeline (see plan.md for detailed schedule)
- All physics implementations should be validated against published nuclear data
- Performance is critical - use OpenMP for parallelization when possible
- The Transport class is the main entry point for particle simulation
- Cross-section functions are pluggable via std::function interfaces
- Extensive use of validation and error checking throughout the transport pipeline

## File Organization

```
cpp/
├── include/           # Header files (.hpp)
├── src/              # Implementation files (.cpp)  
└── tests/            # Unit tests (test_*.cpp)

python/
├── mcshield/         # Main Python package
├── visualization/    # Plotting tools
└── web-app/          # Streamlit interface

Root files:
├── setup.py          # Build configuration
├── pyproject.toml    # Modern Python packaging
├── requirements.txt  # Dependencies
└── plan.md          # Development timeline
```

This is an active development project focused on demonstrating advanced computational physics and software engineering skills for nuclear industry applications.

## Git Commit Guidelines

When creating git commits, do not mention Claude or AI assistance in commit messages. Focus on the technical changes made.

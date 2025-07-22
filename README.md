# Monte Carlo Radiation Transport Simulator

Work in Progress - Under Active Development

A professional-grade Monte Carlo radiation transport simulator built in C++ with Python bindings and a Streamlit web interface. This project demonstrates advanced computational physics and software engineering skills for nuclear industry applications.

## Project Overview

This simulator models the transport of radiation through matter using Monte Carlo methods, supporting multiple particle types (photons, neutrons, electrons) with realistic physics interactions. The project showcases:

- **High-performance C++ computational engine** with OpenMP parallelisation
- **Realistic nuclear physics** with energy-dependent cross-sections
- **Flexible geometry handling** (boxes, cylinders, spheres, composite materials)
- **Python integration** via pybind11 for data analysis and visualisation
- **Professional web interface** using Streamlit for interactive simulations

## Applications

Radiation shielding calculations, dose assessment, and nuclear safety analysis - directly relevant to the nuclear industry and facilities like Sellafield.

## Current Status

- See plan.md for progress

## Planned Features

- Multi-particle Monte Carlo transport (photons, neutrons, electrons, protons)
- Energy-dependent interaction cross-sections
- Advanced geometry support with material composition
- Parallel processing capabilities
- Comprehensive Python interface with scientific computing integration
- Interactive web application for simulation and visualisation
- Validation against analytical solutions and published data

## Technology Stack

- **C++14/17** - Core simulation engine
- **Python 3.8+** - Data analysis and interface
- **pybind11** - C++/Python integration
- **NumPy/Matplotlib** - Scientific computing and plotting
- **Streamlit** - Web application framework
- **OpenMP** - Parallel processing

## Quick Start

```bash
# Clone and setup
git clone <repository-url>
cd monte_carlo_shielding

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Mac/Linux
# or venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Build and test (coming soon)
python setup.py build
python -m pytest tests/
```

## Development Timeline

This is a 5-week intensive development project (January-February 2025) designed to demonstrate advanced programming skills and nuclear physics knowledge for graduate applications in the nuclear industry.

## Author

Lucas Holik - 4th Year MPhys Physics, University of Manchester  
Specialising in Nuclear and Medical Physics

## Licence

This project is developed for educational and portfolio purposes.

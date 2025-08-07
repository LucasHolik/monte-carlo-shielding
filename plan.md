# Monte Carlo Radiation Transport Simulator

## 5-Week Development Plan for Sellafield Graduate Application

---

## **Project Overview**

### **Objective**

Build a professional-grade Monte Carlo radiation transport simulator with C++ computational engine and Python/Streamlit web interface to demonstrate advanced programming skills and nuclear physics knowledge for Sellafield Technical Graduate Programme application.

### **Key Features**

- **Multi-particle transport**: Photons, neutrons, electrons, protons
- **Realistic physics**: Energy-dependent cross-sections, multiple interaction types
- **Flexible geometry**: Boxes, cylinders, spheres, composite materials
- **High performance**: C++ core with OpenMP parallelization
- **User-friendly interface**: Python bindings with Streamlit web application
- **Professional visualisation**: 3D geometry, particle tracks, dose distributions

### **Technical Architecture**

```txt
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   C++ Engine    │    │  Python Bindings │    │ Streamlit Web   │
│                 │    │                  │    │   Interface     │
│ • Particle      │◄──►│ • pybind11       │◄──►│ • User controls │
│   transport     │    │ • NumPy arrays   │    │ • Visualisation │
│ • Physics       │    │ • Data analysis  │    │ • Results       │
│ • Geometry      │    │ • Plotting       │    │ • Comparison    │
│ • Performance   │    │ • Validation     │    │ • Export        │
└─────────────────┘    └──────────────────┘    └─────────────────┘
```

---

## **5-Week Detailed Schedule**

### **Work Pattern**: 5 days/week, weekends free for problem-solving and catch-up

---

## **WEEK 1: Foundation & Basic Framework**

### **Monday - Day 1: Project Setup & Environment**

- [x] Install development environment:
  - Visual Studio/CLion for C++
  - Python 3.8+ with pip
  - Git for version control
- [x] Install dependencies:

  ```bash
  pip install pybind11 numpy matplotlib streamlit plotly
  ```

- [x] Create project directory structure:

  ```txt
  monte_carlo_shielding/
  ├── cpp/
  │   ├── src/
  │   ├── include/
  │   └── tests/
  ├── python/
  │   ├── mcshield/
  │   ├── visualization/
  │   └── web_app/
  ├── data/
  ├── docs/
  ├── tests/
  ├── pyproject.toml     # Modern Python packaging configuration
  ├── setup.py           # Build system for C++ extensions
  ├── .gitignore
  └── README.md
  ```

- [x] Set up basic `setup.py` and build system
- [x] Create simple "Hello World" pybind11 example to test setup
- [x] Initialize Git repository with proper .gitignore
- **Success Metric**: Working build system that compiles C++ and imports in Python

### **Tuesday - Day 2: Core Data Structures**

- [x] Design and implement `Vector3D` class:

  ```cpp
  class Vector3D {
      double x, y, z;
      Vector3D operator+(const Vector3D& other);
      Vector3D operator*(double scalar);
      double dot(const Vector3D& other);
      Vector3D normalize();
  };
  ```

- [x] Create comprehensive `Particle` class:

  ```cpp
  class Particle {
      Vector3D position, direction;
      double energy, weight;
      ParticleType type;
      bool alive;
  };
  ```

- [x] Implement basic `Material` class with atomic properties
- [x] Create `RandomNumberGenerator` class with proper seeding
- [x] Write unit tests for all basic functionality
- [x] Document all classes with clear comments
- **Success Metric**: All core classes working with passing unit tests

### **Wednesday - Day 3: Basic Geometry Engine**

- [x] Implement `Box` geometry class (rectangular shielding):

  ```cpp
  class Box {
      Vector3D min_corner, max_corner;
      Material material;
      bool isInside(const Vector3D& point);
      double distanceToSurface(const Vector3D& pos, const Vector3D& dir);
  };
  ```

- [x] Add boundary detection methods and ray-box intersection algorithms
- [x] Create simple `Geometry` container class for multiple regions
- [x] Test particle-boundary interactions thoroughly
- [x] Add material assignment to geometry regions
- [x] Implement boundary crossing logic
- **Success Metric**: Particles can navigate through box geometry with proper material transitions

### **Thursday - Day 4: Monte Carlo Sampling Foundation**

- [x] Research and implement Monte Carlo sampling techniques:
  - Exponential sampling for interaction distances
  - Uniform spherical direction sampling
  - Energy distribution sampling
- [x] Create sampling utility functions:

  ```cpp
  double sampleExponential(double lambda);
  Vector3D sampleIsotropicDirection();
  double sampleEnergy(const EnergySpectrum& spectrum);
  ```

- [x] Implement and test random number quality
- [x] Add statistical validation for sampling distributions
- [x] Create debugging tools for sampling verification
- **Success Metric**: Robust, validated sampling foundation with statistical tests

### **Friday - Day 5: Basic Particle Transport**

- [x] Implement core particle transport loop:

  ```cpp
  void Transport::trackParticle(Particle& particle) {
      while (particle.alive) {
          double distance = sampleDistanceToInteraction();
          moveParticle(particle, distance);
          if (hitBoundary()) handleBoundary();
          else performInteraction();
      }
  }
  ```

- [x] Add distance-to-interaction sampling
- [x] Create particle tracking through geometry
- [x] Implement boundary crossing and material changes
- [x] Test particle movement and geometry navigation
- [x] Add basic particle history logging
- **Success Metric**: Particles can move through materials, hit boundaries, and be tracked

---

## **WEEK 2: Monte Carlo Physics Implementation**

### **Monday - Day 8: Photon Physics Foundation**

- [x] Research photon interaction cross-sections from NIST data
- [x] Implement photoelectric absorption:

  ```cpp
  bool performPhotoelectricAbsorption(Particle& photon, const Material& mat);
  ```

- [x] Add Compton scattering physics (Klein-Nishina formula):

  ```cpp
  void performComptonScattering(Particle& photon, const Material& mat);
  ```

- [x] Create cross-section calculation functions
- [x] Test individual interaction physics against known values
- [x] Add energy and angle sampling for scattered photons
- **Success Metric**: Individual photon interactions working with correct physics

### **Tuesday - Day 9: Complete Photon Physics**

- [x] Add pair production for high-energy photons (>1.022 MeV)
- [x] Implement Rayleigh (coherent) scattering
- [x] Create total cross-section calculations:

  ```cpp
  double getTotalCrossSection(double energy, const Material& mat);
  InteractionType sampleInteractionType(double energy, const Material& mat);
  ```

- [x] Add interaction type sampling based on relative probabilities
- [x] Test complete photon physics package
- [x] Validate against published cross-section data
- **Success Metric**: Complete, validated photon interaction physics

### **Wednesday - Day 10: Energy-Dependent Cross-Sections**

- [ ] Implement cross-section database structure:

  ```cpp
  class CrossSectionDatabase {
      std::map<int, std::vector<EnergyPoint>> data;
      double interpolate(int Z, double energy, InteractionType type);
  };
  ```

- [ ] Add energy interpolation (log-log interpolation)
- [ ] Create material composition handling for mixtures
- [ ] Load real nuclear data from NIST or ENDF sources
- [ ] Test cross-section calculations vs. published data
- [ ] Add density and temperature corrections
- **Success Metric**: Realistic, energy-dependent physics matching published data

### **Thursday - Day 11: Integrated Photon Monte Carlo**

- [ ] Integrate all physics with transport engine
- [ ] Implement complete photon history tracking:

  ```cpp
  SimulationResults runPhotonSimulation(const Source& source, const Geometry& geom, int n_histories);
  ```

- [ ] Add energy deposition tallies and flux calculations
- [ ] Create particle termination handling (absorption, escape)
- [ ] Test full photon transport simulation
- [ ] Add statistical uncertainty calculations
- **Success Metric**: Working end-to-end photon Monte Carlo simulator

### **Friday - Day 12: Validation & Benchmarking**

- [ ] Create validation test cases using analytical solutions
- [ ] Compare results with Beer-Lambert law for simple cases
- [ ] Test different materials (lead, concrete, water) and energies
- [ ] Benchmark against published shielding data
- [ ] Debug and fix any physics discrepancies
- [ ] Document photon transport module comprehensively
- [ ] Create performance benchmarks
- **Success Metric**: Validated photon transport matching analytical/published results

---

## **WEEK 3: Multi-Particle Support & Advanced Features**

### **Monday - Day 15: Neutron Physics Implementation**

- [ ] Research neutron interaction cross-sections and libraries
- [ ] Implement elastic scattering in center-of-mass system:

  ```cpp
  void performElasticScattering(Particle& neutron, const Material& mat);
  ```

- [ ] Add neutron absorption reactions (n,γ), (n,p), (n,α)
- [ ] Create neutron energy loss calculations
- [ ] Implement thermal neutron physics
- [ ] Test neutron physics components against known data
- **Success Metric**: Basic neutron transport physics working correctly

### **Tuesday - Day 16: Electron Transport Physics**

- [ ] Implement electron energy loss (Bethe-Bloch formula):

  ```cpp
  double calculateEnergyLoss(const Particle& electron, const Material& mat, double distance);
  ```

- [ ] Add multiple scattering (Molière theory or simplified)
- [ ] Create bremsstrahlung physics (secondary photon production)
- [ ] Add electron range calculations
- [ ] Implement electron-electron interactions
- [ ] Test electron transport components
- **Success Metric**: Realistic electron transport with energy loss

### **Wednesday - Day 17: Proton Transport & Advanced Geometry**

- [ ] Implement basic proton energy loss and range calculations
- [ ] Add Coulomb scattering for charged particles
- [ ] Create cylindrical geometry class:

  ```cpp
  class Cylinder {
      Vector3D center, axis;
      double radius, height;
      bool isInside(const Vector3D& point);
  };
  ```

- [ ] Add spherical geometry support
- [ ] Create composite geometry handling (multiple regions)
- [ ] Test complex geometries with particle transport
- **Success Metric**: Multi-particle transport through complex geometries

### **Thursday - Day 18: Performance Optimization**

- [ ] Add OpenMP parallelization for particle histories:

  ```cpp
  #pragma omp parallel for
  for(int i = 0; i < n_histories; i++) {
      transportParticle(source.sample());
  }
  ```

- [ ] Optimize memory usage and data structures
- [ ] Implement basic variance reduction techniques:
  - Russian roulette for low-weight particles
  - Particle splitting in important regions
- [ ] Profile code and identify bottlenecks
- [ ] Test parallel performance and scaling
- **Success Metric**: Parallelized simulation with significant speedup

### **Friday - Day 19: Advanced Features & Testing**

- [ ] Add geometry visualization helpers
- [ ] Implement particle track recording for visualization
- [ ] Create comprehensive test suite for all physics
- [ ] Add simulation result analysis tools
- [ ] Implement error checking and validation
- [ ] Document all advanced features
- [ ] Performance testing and optimization
- **Success Metric**: Robust, well-tested simulation engine

---

## **WEEK 4: Python Integration & Core Interface**

### **Monday - Day 22: Python Bindings Foundation**

- [ ] Create comprehensive pybind11 interface:

  ```cpp
  PYBIND11_MODULE(mcshield, m) {
      py::class_<Particle>(m, "Particle")
          .def(py::init<double, ParticleType>())
          .def_readwrite("energy", &Particle::energy);
      
      py::class_<Simulation>(m, "Simulation")
          .def("run", &Simulation::run);
  }
  ```

- [ ] Expose all major C++ classes to Python
- [ ] Test basic C++ function calls from Python
- [ ] Add Python-friendly constructor overloads
- [ ] Debug any compilation or linking issues
- **Success Metric**: C++ simulation engine accessible from Python

### **Tuesday - Day 23: Python Data Interface**

- [ ] Complete all Python bindings for simulation results
- [ ] Add numpy array integration for large datasets:

  ```python
  results = sim.run_simulation(n_particles=100000)
  dose_array = np.array(results.dose_distribution)
  ```

- [ ] Create Python-friendly result structures
- [ ] Test large data transfer between C++ and Python
- [ ] Add comprehensive error handling for Python interface
- [ ] Create Python wrapper classes for ease of use
- **Success Metric**: Seamless, efficient data flow C++ ↔ Python

### **Wednesday - Day 24: Python Visualization Tools**

- [ ] Create dose distribution plotting functions:

  ```python
  def plot_dose_distribution(results):
      fig, ax = plt.subplots()
      ax.plot(results.thickness, results.dose_rate)
      ax.set_xlabel('Thickness (cm)')
      ax.set_ylabel('Dose Rate')
  ```

- [ ] Add 3D geometry visualization (matplotlib/plotly)
- [ ] Implement particle track plotting
- [ ] Create energy spectrum analysis tools
- [ ] Add statistical analysis functions
- [ ] Test all visualization functions with real data
- **Success Metric**: Professional visualization capabilities

### **Thursday - Day 25: Advanced Python Tools**

- [ ] Create simulation comparison tools
- [ ] Add optimization functions (finding optimal shielding thickness)
- [ ] Implement data export capabilities (CSV, JSON)
- [ ] Create validation and benchmarking utilities
- [ ] Add batch simulation capabilities
- [ ] Document all Python modules comprehensively
- **Success Metric**: Complete Python toolkit for analysis

### **Friday - Day 26: Python Testing & Integration**

- [ ] Create comprehensive Python test suite
- [ ] Test all Python-C++ integration
- [ ] Add example notebooks demonstrating usage
- [ ] Validate Python results against C++ direct calls
- [ ] Performance testing of Python interface
- [ ] Create Python package structure for easy installation
- **Success Metric**: Robust, tested Python interface ready for web app

---

## **WEEK 5: Streamlit Web Application & Professional Polish**

### **Monday - Day 29: Core Streamlit Application**

- [ ] Set up basic Streamlit application framework:

  ```python
  import streamlit as st
  import mcshield
  
  st.title("Monte Carlo Radiation Shielding Calculator")
  
  # Sidebar controls
  st.sidebar.header("Simulation Parameters")
  material = st.sidebar.selectbox("Material", ["Lead", "Concrete", "Steel"])
  ```

- [ ] Create sidebar for all simulation parameters
- [ ] Add particle type selection (photon, neutron, electron)
- [ ] Implement geometry configuration controls
- [ ] Add input validation and error handling
- [ ] Test parameter input and basic interface
- **Success Metric**: Working parameter input interface

### **Tuesday - Day 30: Simulation Integration & Basic Results**

- [ ] Integrate C++ simulation engine with Streamlit:

  ```python
  if st.button("Run Simulation"):
      with st.spinner("Running Monte Carlo..."):
          results = mcshield.run_simulation(params)
      
      col1, col2 = st.columns(2)
      col1.metric("Transmission", f"{results.transmission:.1%}")
      col2.metric("Dose Reduction", f"{results.attenuation:.1f}x")
  ```

- [ ] Add progress bars and real-time status updates
- [ ] Create results display section
- [ ] Implement basic plotting integration
- [ ] Add simulation history and logging
- [ ] Test end-to-end simulation workflow
- **Success Metric**: Working simulation workflow in web app

### **Wednesday - Day 31: Advanced Visualization**

- [ ] Add 3D geometry visualization in Streamlit:

  ```python
  fig = create_3d_geometry_plot(geometry)
  st.plotly_chart(fig, use_container_width=True)
  ```

- [ ] Create interactive particle track displays
- [ ] Implement energy spectrum plots
- [ ] Add dose rate vs. thickness curves
- [ ] Create material comparison visualizations
- [ ] Test all interactive features
- **Success Metric**: Rich, interactive visualizations

### **Thursday - Day 32: Professional Features & UX**

- [ ] Add comparison tools (multiple materials/configurations):

  ```python
  comparison_mode = st.checkbox("Compare Materials")
  if comparison_mode:
      materials = st.multiselect("Select materials", material_list)
  ```

- [ ] Implement result export functionality (PDF, CSV)
- [ ] Add simulation presets for common scenarios
- [ ] Create help documentation and tooltips
- [ ] Polish UI/UX with professional styling
- [ ] Add loading animations and user feedback
- **Success Metric**: Professional, user-friendly interface

### **Friday - Day 33: Final Polish & Documentation**

- [ ] Comprehensive testing of entire application
- [ ] Create user guide and documentation
- [ ] Add example scenarios and validation cases
- [ ] Performance optimization for web deployment
- [ ] Create demo video/screenshots
- [ ] Prepare GitHub repository for public viewing
- [ ] Final code cleanup and commenting
- **Success Metric**: Production-ready, portfolio-quality application

---

## **Tools & Technologies**

### **Development Environment**

- **C++ Compiler**: GCC 9+ or MSVC 2019+
- **IDE**: Visual Studio, CLion, or VS Code
- **Python**: 3.8+ with pip package manager
- **Git**: Version control

### **C++ Libraries**

```cpp
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <memory>
#include <omp.h>  // OpenMP for parallelization
```

### **Python Libraries**

```python
# Core
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st

# Building
import pybind11
from setuptools import setup, Extension

# Testing
import pytest
import unittest
```

### **Build System**

```python
# setup.py
from pybind11.setup_helpers import Pybind11Extension, build_ext
ext_modules = [
    Pybind11Extension(
        "mcshield",
        ["cpp/src/monte_carlo.cpp"],
        cxx_std=14,
    ),
]
```

---

## **Success Metrics & Milestones**

### **Weekly Milestones**

- **Week 1**: ✅ Basic C++ framework with particle transport
- **Week 2**: ✅ Complete photon Monte Carlo simulation
- **Week 3**: ✅ Multi-particle support with advanced geometry  
- **Week 4**: ✅ Complete Python interface and visualization
- **Week 5**: ✅ Professional web application ready for demo

### **Daily Success Criteria**

Each day should achieve:

- ✅ Specific functionality working (as detailed above)
- ✅ Code committed to Git with clear messages
- ✅ Basic tests passing for new functionality
- ✅ Documentation updated for new features
- ✅ No breaking changes to existing code

### **Final Project Requirements**

- ✅ Multi-particle Monte Carlo transport (photons + neutrons minimum)
- ✅ Realistic physics with energy-dependent cross-sections
- ✅ Parallel processing capability
- ✅ Python interface with comprehensive bindings
- ✅ Professional web application
- ✅ Validation against analytical solutions
- ✅ Comprehensive documentation
- ✅ GitHub repository ready for portfolio

---

## **Risk Management & Backup Plans**

### **If Behind Schedule**

**Priority 1 (Must Have)**:

- Photon transport with basic physics
- Simple geometry (boxes only)
- Basic Python bindings
- Simple Streamlit interface

**Priority 2 (Should Have)**:

- Neutron transport
- Advanced geometry (cylinders, spheres)
- Advanced visualizations
- Performance optimization

**Priority 3 (Nice to Have)**:

- Electron and proton transport
- Variance reduction techniques
- Advanced web features
- Comprehensive validation

### **Common Issues & Solutions**

**Build System Problems**:

- Use conda environment for consistent dependencies
- Test on multiple platforms if possible
- Keep backup simple Makefile approach

**Physics Implementation**:

- Start with simplified models, add complexity gradually
- Validate each physics component individually
- Use well-documented references (Turner, Evans)

**Performance Issues**:

- Profile early and often
- Focus on algorithmic improvements before micro-optimization
- Use OpenMP for embarrassingly parallel particle histories

### **Weekend Catch-up Strategy**

- **Saturday**: Focus on critical path items that are behind
- **Sunday**: Testing, debugging, and documentation
- **Emergency**: Reduce scope rather than compromise quality

---

## **Expected Outcomes**

### **Technical Skills Demonstrated**

- **Advanced C++ programming** (OOP, memory management, performance)
- **Monte Carlo methods** (sampling, statistics, variance reduction)
- **Nuclear physics** (radiation transport, cross-sections, interactions)
- **Software architecture** (modular design, interfaces, testing)
- **Python integration** (pybind11, numpy, scientific computing)
- **Web development** (Streamlit, visualization, UX design)
- **Project management** (planning, version control, documentation)

### **Portfolio Value**

- **Computational physics capability**
- **Multi-language proficiency**
- **Industry-relevant knowledge**
- **Professional software development**
- **Problem-solving skills**
- **Communication through visualization**

### **Sellafield Relevance**

This project directly demonstrates skills needed for Sellafield Technical graduates:

- **Radiation transport calculations** (core to nuclear safety)
- **Computer simulation** (essential for modern nuclear engineering)
- **Data analysis and visualization** (key to operational support)
- **Software development** (increasingly important in technical roles)
- **Physics understanding** (fundamental to nuclear applications)

---

## **Final Deliverables**

1. **GitHub Repository** with professional README
2. **Working C++ Monte Carlo engine** with comprehensive physics
3. **Python package** with full bindings and utilities
4. **Streamlit web application** for interactive simulation
5. **Technical documentation** explaining physics and implementation
6. **Validation report** comparing results to analytical solutions
7. **Demo video** showing key features and capabilities
8. **User guide** for installation and usage

This project will serve as a showcase piece for technical interviews and demonstrate the advanced programming and physics skills valued by Sellafield's Technical Graduate Programme.

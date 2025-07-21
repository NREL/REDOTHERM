# REDOTHERM

**REDOTHERM** is an open-source, MATLAB-based thermodynamic modeling framework developed to evaluate the performance of redox-active materials for water (H2O) and carbon dioxide (CO2) splitting. It includes models of all major unit operations and supports comparative analysis of different redox-active material candidates. The model is tailored for systems of moving oxide under a parallel/cocurrent flow (PF) and countercurrent flow (CF) configurations. Unvalidated mixed flow reactor (MFR, also known as CSTR) model is also included as an optional addition.

The development of REDOTHERM was done by Alon Lidor and Janna Martinek at the National Renewable Energy Laboratory.

---

## LICENSE

[License](https://github.com/NREL/REDOTHERM/blob/main/LICENSE)

---

## 🔧 Features

- Thermodynamic performance evaluation of redox-active materials
- Modeling of full system including key unit operations
- Modular structure for easy extension and customization
- MATLAB-based: leverage familiar scripting and visualization tools
- Open-source and community-driven
- Multiple technology options for different unit operations

---

## 📥 Installation

### Requirements
- MATLAB R2023a or later (recommended)
- For optmization capabilities: Global Optimization Toolbox
- [Cantera](https://cantera.org/) - chemical kinetics and thermodynamics toolkit
- [CoolProp](https://coolprop.org/) - thermophysical properties package

### Clone the repository
```bash
git clone https://github.com/NREL/REDOTHERM.git
cd REDOTHERM
```

---

### 🧱 Repository Structure

```bash
REDOTHERM/
│
├── core/                                  # Core thermodynamic functions
├── functions/                             # Functions for optimization and cycle performance
├── materials/                             # Thermodynamic properties of redox-active materials
├── examples/                              # Usage examples and test cases
├── utils/                                 # Utility functions (plotting, integration, etc.)
├── Redox_Countercurrent_Thermo_Main.m     # Running REDOTHERM for a single or parametric run
├── Redox_Countercurrent_Thermo_Opt.m      # Running REDOTHERM for a single optimization analysis
├── Redox_Countercurrent_Thermo_Opt_Para.m # Running REDOTHERM for multiple optimization runs
├── LICENSE            # License file
└── README.md          # Readme file
```

---

## 📚 Documentation

The model structure, functionality of each script, and the derviation of all the equations are provided in a manuscript current under review. A URL will be provided here once the paper is accepted and online. In the meantime, for any questions about the usage or methodology, please email [Dr. Alon Lidor](mailto:alon.lidor@nrel.gov)

# UDF Code for Heat Transfer Enhancement in Ferrofluids under External Magnetic Field

This repository contains a **User Defined Function (UDF)** written in **C++** for simulations in **ANSYS Fluent**. The project focuses on enhancing heat transfer using ferrofluids under the influence of an external magnetic field. This UDF was developed as part of my final year **Major Project**.

## Project Overview

The aim of this project is to explore how ferrofluids (colloidal suspensions of magnetic nanoparticles) can enhance heat transfer when subjected to an external magnetic field. The UDF code modifies the heat transfer calculations in **ANSYS Fluent**, incorporating the effects of both **Brownian motion** and **magnetic field** interactions to predict the fluid behavior more accurately.

### Key Features:

- **Temperature-dependent Thermal Conductivity**: The code accounts for temperature variations in ferrofluid properties, including the Brownian thermal conductivity.
- **External Magnetic Field Effects**: Incorporates the influence of a magnetic field on ferrofluid movement and heat transfer.
- **Dynamic Property Updates**: Ensures that fluid properties are dynamically updated based on local temperature and magnetic field strength.

## Requirements

- **ANSYS Fluent** (with UDF support)
- **C++ compiler** (such as GCC or Visual Studio)
- Basic knowledge of CFD (Computational Fluid Dynamics) and UDF implementation in ANSYS.

   
2. **Open the project in ANSYS Fluent** and load the UDF.
   
3. **Compile the UDF** using the in-built ANSYS UDF compiler or an external C++ compiler:
    - Navigate to `Define -> User-Defined -> Functions -> Interpreted/Compiled`.
   
4. **Assign Boundary Conditions**: Set up your simulation with appropriate boundary conditions, ensuring that the external magnetic field and fluid properties are defined correctly.

5. **Run the Simulation**: Once the UDF is compiled and applied, run the simulation and analyze the enhanced heat transfer effects.

## Simulation Workflow

1. **Pre-processing**:
    - Set up geometry and mesh in ANSYS Fluent.
    - Define ferrofluid properties and boundary conditions.
  
2. **UDF Integration**:
    - Import and compile the UDF.
    - Set the UDF for relevant boundary conditions (e.g., wall heat flux, magnetic field influence).

3. **Post-processing**:
    - Visualize the temperature distribution and velocity fields to observe the effects of ferrofluid and magnetic fields on heat transfer.

## Contact

For any inquiries or further collaboration, feel free to reach out:

- **Asmanya**  
- Email: asmanya.sh@gmail.com


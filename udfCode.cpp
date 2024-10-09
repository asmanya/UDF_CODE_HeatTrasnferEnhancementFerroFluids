#include <bits/stdc++.h>
#include "udf.h"
#include "math.h"


// Constants
#define MU_0 1.256637061e-6 // Permeability of free space
#define RHO_F 2446.0 // Density of fluid
#define CHI_0 1.0 // Initial magnetic susceptibility
#define BETA 0.1 // Coefficient for temperature dependence
#define T_0 298.15 // Reference temperature in Kelvin
#define K_STATIC 0.6 // Static thermal conductivity (example value)
#define K_BROWN 5.0e4 // Coefficient for Brownian thermal conductivity (example value)
#define PHI 0.1 // Volume fraction (example value)
#define RHO_NP 1000.0 // Density of nanoparticle
#define C_P_F 4200.0 // Specific heat of fluid
#define C_P_NP 800.0 // Specific heat of nanoparticle


// __________________________________________________________________________________________


// Function to calculate the magnetic field B
DEFINE_PROPERTY(magnetic_field, cell, thread)
{
    real H = C_H[cell]; // Magnetic field intensity (from Fluent)
    real chi_m = C_CHI_M[cell]; // Magnetic susceptibility
    real M = chi_m * H; // Magnetization
    return MU_0 * (M + H); // B = μ_0 (M + H)
}


// __________________________________________________________________________________________


// Function to calculate the temperature dependence of magnetic susceptibility
DEFINE_PROPERTY(susceptibility, cell, thread)
{
    real T = C_T[cell]; // Temperature
    return CHI_0 / (1 + BETA * (T - T_0)); // χ_m(T) = χ_0 / (1 + β(T - T_0))
}


// __________________________________________________________________________________________


// Function for calculating the magnetic force F_k
DEFINE_SOURCE(magnetic_force, cell, thread, dS, eqn)
{
    real T = C_T(cell);            // Temperature at the cell
    real T_0 = 298.15;             // Reference temperature (K)
    real H = H_FIELD;              // Magnetic field intensity (replace with actual value or calculation)
    real chi_m = CHI_M;            // Magnetic susceptibility (replace with actual value or calculation)
    real chi_0 = CHI_0;            // Chi_0 (replace with actual value or calculation)
    real beta = BETA;              // Beta (replace with actual value)
    real alpha = ALPHA;            // Alpha (replace with actual value)
    real mu_0 = 4 * M_PI * 1e-7;   // Magnetic permeability in vacuum (H/m or T·m/A)

    // Temperature gradient (∂T/∂y) in the y-direction
    real dT_dy = C_T_G(cell)[1];   // Gradient of temperature in the y-direction

    // Updated equation for F_k
    real F_k = mu_0 * alpha * chi_m * H * H * (-chi_0 * beta) / 
               pow(1 + beta * (T - T_0), 2) * dT_dy;

    dS[eqn] = 0.0;  // No contribution to the source derivative
    
    return F_k;
}


// __________________________________________________________________________________________


// Function to calculate the momentum equations for the fluid
DEFINE_MASS_TRANSFER(mass_transfer, cell, thread)
{
    real rho = RHO_F; // Density of fluid
    real C_p = C_P_F; // Specific heat of fluid
    real T = C_T[cell]; // Temperature
    real P = C_P[cell]; // Pressure

    // Continuity equation: ∇·(ρV) = 0
    real continuity = rho * C_V[cell]; 

    // Momentum equation: ∇·(ρVV) = -∇P + ∇·(τ_ij) + F_k
    real tau = MU_0 * C_CHI_M[cell]; // Simplified shear stress (placeholder)
    real momentum = -P + tau + force_k(cell, thread, NULL, NULL, 0);

    // Energy equation: ∇·(ρC_pT) = ∇·(k∇T)
    real energy = rho * C_p * T;

    return continuity + momentum + energy;
}


// __________________________________________________________________________________________


// Function to calculate the thermal conductivity
DEFINE_PROPERTY(thermal_conductivity, cell, thread)
{
    real k_static = K_STATIC; // Static thermal conductivity
    real k_brownian = K_BROWN * PHI * RHO_NP * C_P_F; // Brownian thermal conductivity

    // Total thermal conductivity k_nf = k_static + k_brownian
    return k_static + k_brownian;
}


// __________________________________________________________________________________________


// Function to calculate stress tensor
DEFINE_PROPERTY(stress_tensor, cell, thread)
{
    real mu = 0.001; // Dynamic viscosity (example value)
    real du_dx = (C_V[cell + 1] - C_V[cell - 1]) / (2 * DX); // Velocity gradient
    real du_dy = (C_V[cell + 1] - C_V[cell - 1]) / (2 * DY); // Velocity gradient in y direction

    // τ_ij = μ (∂v_i/∂x_j + ∂v_j/∂x_i) - (2/3)μδ_ij(∂v_k/∂x_k)
    real tau_ij = mu * (du_dx + du_dy) - (2.0 / 3.0) * mu * (du_dx + du_dy);
    
    return tau_ij; // Return the stress tensor
}


// __________________________________________________________________________________________


// Boundary condition function for velocity and temperature
DEFINE_BC(boundary_conditions, cell, thread)
{
    // Conditions for fluid
    real u = 0.0; // Velocity in x direction
    real v = 0.0; // Velocity in y direction
    real T_f = C_T[cell]; // Fluid temperature
    real T_s = C_T[cell]; // Solid temperature

    // Assuming Dirichlet conditions for velocity and temperature
    if (THREAD_ID(thread) == 1) { // Assuming thread ID 1 is a wall
        C_V[cell] = u; // Set velocity
        C_T[cell] = T_f; // Set fluid temperature
    }

    // Heat flux condition
    real k_f = thermal_conductivity(cell, thread);
    real k_s = 0.5; // Solid thermal conductivity (example value)
    
    // -k_f * (∂T_f/∂n) = -k_s * (∂T_s/∂n)
    // This will need to be implemented as a heat transfer boundary condition in Fluent.

    return 0.0; // No additional source term
}


// __________________________________________________________________________________________


// Function for effective density
DEFINE_PROPERTY(effective_density, cell, thread)
{
    return (1 - PHI) * RHO_F + PHI * RHO_NP; // ρ_{nf} = (1 - φ) ρ_f + φ ρ_{np}
}


// __________________________________________________________________________________________


// Function for effective heat capacity
DEFINE_PROPERTY(effective_heat_capacity, cell, thread)
{
    return (1 - PHI) * C_P_F + PHI * C_P_NP; // (ρC_p)_{nf} = (1 - φ) (ρC_p)_f + φ (ρC_p)_{np}
}


// __________________________________________________________________________________________


// Function for effective viscosity
DEFINE_PROPERTY(effective_viscosity, cell, thread)
{
    return RHO_F * (1 + 2.5 * PHI); // μ_{nf} = μ_f (1 + 2.5 φ)
}


// __________________________________________________________________________________________


// Function for effective thermal conductivity of nanoparticle
DEFINE_PROPERTY(nanoparticle_thermal_conductivity, cell, thread)
{
    real k_np = 0.5; // Thermal conductivity of nanoparticle (example value)
    
    return k_static + ((k_np + 2 * k_static) - 2 * PHI * (k_static - k_np)) / ((k_np + 2 * k_static)
        + PHI * (k_static - k_np)); 
}


// __________________________________________________________________________________________


// Function for calculating Brownian thermal conductivity
DEFINE_PROPERTY(brownian_conductivity, cell, thread)
{
    real T = C_T[cell];  // Temperature at the cell
    real k = 0.026;      // Example thermal conductivity of fluid (W/m·K) (replace with actual value)
    real rho_f = RHO_F;  // Fluid density
    real C_p_f = C_P_F;  // Specific heat capacity of fluid
    real D_np = 1e-9;    // Example nanoparticle diameter (m) (replace with actual value)
    real rho_np = 5000;  // Example nanoparticle density (kg/m³) (replace with actual value)

    // Calculate the temperature-dependent function g(Phi, T)
    real g_Phi_T = (-6.04 * PHI + 0.4705) * T + 1722.3 * PHI - 134.63;

    // Updated equation for k_Brownian = 5 * 10^4 * beta * Phi * rho_f * C_p_f * sqrt(k * T / (rho_np * D_np)) * g(Phi, T)
    return 5.0e4 * BETA * PHI * rho_f * C_p_f * sqrt((k * T) / (rho_np * D_np)) * g_Phi_T;
}


// __________________________________________________________________________________________


// Function for temperature-dependent function g(Phi, T)
DEFINE_PROPERTY(temp_dependent_function, cell, thread)
{
    real T = C_T[cell];  // Current temperature at the cell
    
    // Updated equation for g(Phi, T) = (-6.04 * Phi + 0.4705) * T + 1722.3 * Phi - 134.63
    return (-6.04 * PHI + 0.4705) * T + 1722.3 * PHI - 134.63; 
}


// __________________________________________________________________________________________


// Function for density of fluid
DEFINE_PROPERTY(density_fluid, cell, thread)
{
    return 2446 - 20.6747 + 0.115767 * C_T[cell] * C_T[cell] - 3.12895e-4 * C_T[cell]
        * C_T[cell] * C_T[cell] + 4.0505e-7 * C_T[cell] * C_T[cell] * C_T[cell]
        * C_T[cell] - 2.0546e-10 * C_T[cell] * C_T[cell] * C_T[cell] * C_T[cell] * C_T[cell]; 
}


// __________________________________________________________________________________________


// Function for viscosity of fluid
DEFINE_PROPERTY(viscosity_fluid, cell, thread)
{
    return 2.414e-5 * pow(10, (247.8 / (C_T[cell] - 140))); // μ_f = 2.414 × 10^{-5} × 10^{[247.8 / (T - 140)]}
}


// __________________________________________________________________________________________


// Function for Nusselt number
DEFINE_PROPERTY(nusselt_number, cell, thread)
{
    real h = (Q / (T_W - T_B)); // h = q'' / (T_w - T_b)
    real D = 0.1; // Characteristic length (example value)
    
    return h * D / (K_STATIC * 1.0); // Nu = hD / k
}


// __________________________________________________________________________________________


// Function for Reynolds number
DEFINE_PROPERTY(reynolds_number, cell, thread)
{
    real u = C_V[cell]; // Velocity
    real L = 0.1; // Characteristic length (example value)
    real mu = viscosity_fluid(cell, thread); // Viscosity
    
    return (RHO_F * u * L) / mu; // Re = (ρuL)/μ
}


// __________________________________________________________________________________________


// Function for Prandtl number
DEFINE_PROPERTY(prandtl_number, cell, thread)
{
    real mu = viscosity_fluid(cell, thread); // Dynamic viscosity
    real alpha = K_STATIC / (RHO_F * C_P_F); // Thermal diffusivity
    
    return mu / alpha; // Pr = μ / α
}


// __________________________________________________________________________________________


// Function for Schmidt number
DEFINE_PROPERTY(schmidt_number, cell, thread)
{
    real mu = viscosity_fluid(cell, thread); // Dynamic viscosity
    real D = 1.0; // Diffusion coefficient (example value)
    
    return mu / (RHO_F * D); // Sc = μ / (ρD)
}


// __________________________________________________________________________________________


// Function for additional thermal properties
DEFINE_PROPERTY(additional_thermal_properties, cell, thread)
{
    return C_P_F * (T_W - T_B) / RHO_F; // Function for additional thermal properties
}


// __________________________________________________________________________________________


// Function for calculating heat capacity of the mixture
DEFINE_PROPERTY(heat_capacity_mixture, cell, thread)
{
    return (1 - PHI) * C_P_F + PHI * C_P_NP; // (ρC_p)_{nf} = (1 - φ) (ρC_p)_f + φ (ρC_p)_{np}
}
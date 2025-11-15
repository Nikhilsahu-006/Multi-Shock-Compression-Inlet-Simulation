import math as m
from Theta_beta_mach import beta, theta_max
from Normal_shock import normal_shock, stagnation_temp, stagnation_pressure

def oblique_shock(M1, P1=1, T1=1, theta1=0, theta2=0, gamma=1.4):
    """
    Calculate flow properties across an oblique shock wave.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number (must be > 1)
    P1 : float, optional
        Upstream static pressure (default = 1, normalized units)
    T1 : float, optional
        Upstream static temperature (default = 1, normalized units)
    theta1 : float, optional
        First deflection angle in degrees (default = 0)
    theta2 : float, optional
        Second deflection angle in degrees (default = 0)
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    list or str
        If valid: [M2, P2, T2, T0, P01, P02, beta, Cp] where:
        - M2: Downstream Mach number
        - P2: Downstream static pressure
        - T2: Downstream static temperature
        - T0: Stagnation temperature (constant across shock)
        - P01: Upstream stagnation pressure
        - P02: Downstream stagnation pressure
        - beta: Shock wave angle in degrees
        - Cp: Pressure coefficient
        If invalid: Error message string
    
    Theory:
    -------
    Oblique shocks are analyzed by decomposing the flow into:
    1. Normal component: undergoes normal shock relations
    2. Tangential component: remains unchanged across shock
    
    The normal Mach component Mn1 = M1 * sin(β) determines shock strength.
    """
    # Input validation
    if M1 < 1:
        return "Mach number is less than 1, no oblique shock formed"
    if theta1 <= 0 or theta1 + theta2 <= 0:
        return "Deflection angle theta1 must be greater than 0 for oblique shock"
    if theta1 + theta2 >= theta_max(M1):
        return f"Deflection angle {theta1+theta2} is greater than max deflection {theta_max(M1)}"

    # Solve for shock angle beta using theta-beta-Mach relation
    beta_angle = beta(M1, (theta1 + theta2))  # Method to solve for beta (shock angle)
    beta_rad = m.radians(beta_angle)

    # Calculate normal Mach number component before the shock
    Mn1 = M1 * m.sin(beta_rad)

    # Solve the normal shock equations using the normal Mach component
    sol = normal_shock(Mn1, P1, T1, gamma)

    # Check if the normal shock solution is valid
    if isinstance(sol, str):
        return sol  # Return the error message if the solution is not valid

    # Calculate downstream Mach number using shock geometry
    M2 = sol[0] / m.sin(m.radians(beta_angle - theta1 - theta2))

    # Extract pressure, temperature, and stagnation properties
    P2 = sol[1]
    T2 = sol[2]
    T0 = stagnation_temp(M1, T1, gamma)
    P01 = stagnation_pressure(M1, P1, gamma)
    P02 = stagnation_pressure(M2, P2, gamma)
    
    # Calculate pressure coefficient
    Cp = round((2/(gamma*(M1**2)))*((P2/P1)-1), 3)

    return [round(M2, 5), P2, T2, T0, P01, P02, beta_angle, Cp]

def show_oblique(results):
    """
    Display oblique shock results in a formatted way.
    
    Parameters:
    -----------
    results : list or str
        Output from oblique_shock function
    """
    # Check if the output is a list or error message
    if isinstance(results, list):
        print(f"M2 = {results[0]}")      # Downstream Mach number
        print(f"P2 = {results[1]}")      # Downstream static pressure
        print(f"T2 = {results[2]}")      # Downstream static temperature
        print(f"T0 = {results[3]}")      # Stagnation temperature
        print(f"P01 = {results[4]}")     # Upstream stagnation pressure
        print(f"P02 = {results[5]}")     # Downstream stagnation pressure
        print(f"beta = {results[6]}")    # Shock wave angle (degrees)
        print(f"Cp = {results[7]}")      # Pressure coefficient
    else:
        print(results)

# Key Physical Insights for Oblique Shocks:
# ----------------------------------------
# 1. Decomposition Principle:
#    - Oblique shock = Normal shock (normal component) + Tangential flow
#    - Mn1 = M1 * sin(β) governs shock strength
#    - Tangential velocity component is conserved

# 2. Shock Angle Behavior:
#    - Weak shock: β ≈ μ (Mach angle) for small θ
#    - Strong shock: β → 90° (approaches normal shock)
#    - Maximum θ: Boundary between attached/detached shocks

# 3. Property Changes:
#    - Pressure: Increases (P2 > P1)
#    - Temperature: Increases (T2 > T1)
#    - Density: Increases (ρ2 > ρ1)
#    - Stagnation pressure: Decreases (P02 < P01)
#    - Stagnation temperature: Constant (T02 = T01)

# 4. Pressure Coefficient (Cp):
#    - Dimensionless pressure rise
#    - Cp = (P2 - P1) / (0.5 * γ * P1 * M1²)
#    - Useful for aerodynamic force calculations

# Applications:
# ------------
# - Supersonic airfoils and wings
# - Engine intakes and compression ramps
# - Supersonic inlets
# - Shock wave/boundary layer interactions
# - Supersonic compressor blades


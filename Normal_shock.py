import math as m

def Mach2(M1, gamma=1.4):
    """
    Calculate downstream Mach number after a normal shock.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number (must be > 1 for shock to exist)
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Downstream Mach number M2 (subsonic)
    
    Equation:
    ---------
    M2² = [1 + (γ-1)/2 * M1²] / [γ * M1² - (γ-1)/2]
    """
    term1 = 1 + ((gamma - 1) / 2) * (M1 ** 2)
    term2 = (gamma * (M1 ** 2)) - ((gamma - 1) / 2)
    return round(m.sqrt(term1 / term2), 5)

def static_pressure_ratio(M1, M2, gamma=1.4):
    """
    Calculate static pressure ratio across normal shock.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number
    M2 : float
        Downstream Mach number
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Static pressure ratio P2/P1
    
    Equation:
    ---------
    P2/P1 = (1 + γ * M1²) / (1 + γ * M2²)
    """
    term1 = 1 + gamma * (M1 ** 2)
    term2 = 1 + gamma * (M2 ** 2)
    return round(term1 / term2, 5)

def static_temp_ratio(M1, M2, gamma=1.4):
    """
    Calculate static temperature ratio across normal shock.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number
    M2 : float
        Downstream Mach number
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Static temperature ratio T2/T1
    
    Equation:
    ---------
    T2/T1 = [1 + (γ-1)/2 * M1²] / [1 + (γ-1)/2 * M2²]
    """
    term1 = 1 + ((gamma - 1) / 2) * (M1 ** 2)
    term2 = 1 + ((gamma - 1) / 2) * (M2 ** 2)
    return round(term1 / term2, 5)

def stagnation_temp(M1, T1, gamma=1.4):
    """
    Calculate stagnation temperature.
    
    Parameters:
    -----------
    M1 : float
        Mach number
    T1 : float
        Static temperature
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Stagnation temperature T0
    
    Equation:
    ---------
    T0/T1 = 1 + (γ-1)/2 * M1²
    """
    term1 = 1 + ((gamma - 1) / 2) * (M1 ** 2)
    return round(T1 * term1, 5)

def stagnation_pressure(M1, P1, gamma=1.4):
    """
    Calculate stagnation pressure.
    
    Parameters:
    -----------
    M1 : float
        Mach number
    P1 : float
        Static pressure
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Stagnation pressure P0
    
    Equation:
    ---------
    P0/P1 = [1 + (γ-1)/2 * M1²]^(γ/(γ-1))
    """
    term1 = (1 + ((gamma - 1) / 2) * (M1 ** 2)) ** (gamma / (gamma - 1))
    return round(P1 * term1, 5)

def normal_shock(M1, P1=1, T1=1, gamma=1.4):
    """
    Calculate all flow properties across a normal shock wave.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number (must be > 1 for shock to exist)
    P1 : float, optional
        Upstream static pressure (default = 1, normalized units)
    T1 : float, optional
        Upstream static temperature (default = 1, normalized units)
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    list or str
        If M1 >= 1: [M2, P2, T2, T0, P01, P02] where:
        - M2: Downstream Mach number
        - P2: Downstream static pressure
        - T2: Downstream static temperature
        - T0: Stagnation temperature (constant across shock)
        - P01: Upstream stagnation pressure
        - P02: Downstream stagnation pressure
        If M1 < 1: Error message string
    """
    if M1 < 1:
        return "Mach number is less than 1, no shock formed"
    else:
        M2 = Mach2(M1, gamma)
        P2 = P1 * static_pressure_ratio(M1, M2, gamma)
        T2 = T1 * static_temp_ratio(M1, M2, gamma)
        T0 = stagnation_temp(M1, T1, gamma)
        P01 = stagnation_pressure(M1, P1, gamma)
        P02 = stagnation_pressure(M2, P2, gamma)
        
        return [M2, P2, T2, T0, P01, P02]

def show_normal(results):
    """
    Display normal shock results in a formatted way.
    
    Parameters:
    -----------
    results : list or str
        Output from normal_shock function
    """
    if isinstance(results, list):
        print(f"M2 = {results[0]}")      # Downstream Mach number
        print(f"P2 = {results[1]}")      # Downstream static pressure
        print(f"T2 = {results[2]}")      # Downstream static temperature
        print(f"T0 = {results[3]}")      # Stagnation temperature
        print(f"P01 = {results[4]}")     # Upstream stagnation pressure
        print(f"P02 = {results[5]}")     # Downstream stagnation pressure
    else:
        print(results)

# Key Physical Insights:
# ----------------------
# 1. Normal shocks always decelerate flow from supersonic to subsonic
# 2. Stagnation temperature remains constant across shock (T0)
# 3. Stagnation pressure decreases across shock (P02 < P01) - entropy increase
# 4. Static pressure and temperature increase across shock
# 5. Shock strength increases with upstream Mach number M1

# Typical Applications:
# --------------------
# - Supersonic inlets and diffusers
# - Shock tubes
# - Supersonic wind tunnels
# - Rocket engine nozzles
# - High-speed aircraft intakes
import numpy as np

def isentropic_diffuser(M1, P1=1, T1=1, M2=0.1, gamma=1.4):
    """
    Calculate isentropic flow properties through a diffuser.
    
    This function computes the exit Mach number, pressure, temperature, 
    and area ratio for isentropic flow in a diffuser using gas dynamics relations.
    
    Parameters:
    -----------
    M1 : float
        Inlet Mach number
    P1 : float, optional
        Inlet pressure (default = 1, normalized units)
    T1 : float, optional  
        Inlet temperature (default = 1, normalized units)
    M2 : float, optional
        Exit Mach number (default = 0.1, subsonic exit)
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    list
        [M2, P2, T2, A_At] where:
        - M2: Exit Mach number
        - P2: Exit pressure
        - T2: Exit temperature  
        - A_At: Area ratio (A/A*) where A* is throat area
    
    Equations:
    ----------
    Temperature ratio: T2/T1 = [1 + 0.5*(γ-1)*M1²] / [1 + 0.5*(γ-1)*M2²]
    Pressure ratio: P2/P1 = (T2/T1)^(γ/(γ-1))
    Area ratio: A/A* = (1/M) * [2/(γ+1) * (1 + (γ-1)/2 * M²)]^((γ+1)/(2(γ-1)))
    """
    
    # Calculate temperature ratio using isentropic relation
    T_ratio = ((1 + 0.5 * (gamma - 1) * M1**2) / 
               (1 + 0.5 * (gamma - 1) * M2**2))

    # Calculate pressure ratio using isentropic relation
    P2_P1 = T_ratio**(gamma / (gamma - 1))

    # Calculate exit temperature and pressure
    T2 = round(T_ratio * T1, 6)  # Exit temperature
    P2 = round(P2_P1 * P1, 6)    # Exit pressure

    # Calculate area ratio (A/A*) where A* is throat area
    term = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M1**2)
    A_At = round((1 / M1) * term ** ((gamma + 1) / (2 * (gamma - 1))), 4)

    return [M2, P2, T2, A_At]
import numpy as np
from scipy.optimize import minimize_scalar, fsolve
import math as m

def theta_beta_mach_relation(beta_rad, M1, theta_deg, gamma=1.4):
    """
    Calculate the theta-beta-Mach relation for oblique shocks.
    
    Parameters:
    -----------
    beta_rad : float
        Shock wave angle in radians
    M1 : float
        Upstream Mach number
    theta_deg : float
        Flow deflection angle in degrees
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Residual value (should be zero for valid solution)
    
    Equation:
    ---------
    tan(θ) = 2cot(β) * (M1²sin²(β)-1) / [M1²(γ+cos(2β)) + 2]
    
    This is the fundamental relation linking flow deflection angle (θ),
    shock angle (β), and upstream Mach number (M1).
    """
    theta = np.radians(theta_deg)
    term1 = 2 / np.tan(beta_rad)
    term2 = (M1**2 * np.sin(beta_rad)**2 - 1)
    term3 = (M1**2 * (gamma + np.cos(2 * beta_rad)) + 2)
    
    # Relation between theta, beta, and Mach number
    return np.tan(theta) - (term1 * term2 / term3)

def theta_max(M1, gamma=1.4):
    """
    Calculate the maximum flow deflection angle for attached oblique shock.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Maximum flow deflection angle θ_max in degrees
    
    Physical Significance:
    ----------------------
    - For θ > θ_max: No attached oblique shock solution exists
    - For θ < θ_max: Two possible solutions (weak and strong shock)
    - For θ = θ_max: Single solution (detachment condition)
    
    The function finds the shock angle β that maximizes θ for given M1.
    """
    # Define the theta-beta-Mach relation
    def f(beta_rad):
        term1 = 2 / np.tan(beta_rad)
        term2 = (M1**2 * np.sin(beta_rad)**2 - 1)
        term3 = (M1**2 * (gamma + np.cos(2 * beta_rad)) + 2)
        theta_rad = np.arctan((term1 * term2) / term3)
        return theta_rad

    # Minimize the negative of the function to find the maximum theta
    result = minimize_scalar(lambda beta: -f(beta), bounds=(0.01, np.pi/2), 
                           method='bounded')

    # The maximum theta in radians
    max_theta_rad = f(result.x)

    # Convert from radians to degrees
    max_theta_deg = m.degrees(max_theta_rad)

    return round(max_theta_deg, 10)

def beta(M1, theta_deg, gamma=1.4):
    """
    Calculate shock wave angle for given Mach number and deflection angle.
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number
    theta_deg : float
        Flow deflection angle in degrees
    gamma : float, optional
        Specific heat ratio (default = 1.4 for air)
    
    Returns:
    --------
    float
        Shock wave angle β in degrees
    
    Solution Types:
    ---------------
    - Returns the WEAK shock solution by default (smaller β)
    - For θ < θ_max, two solutions exist: weak and strong shock
    - Strong shock solution has larger β and subsonic downstream flow
    
    Note:
    -----
    Uses numerical root-finding to solve the transcendental theta-beta-Mach equation.
    """
    beta_guess = np.radians(45)  # Guess of 45 degrees in radians
    
    # Use fsolve to find the root (shock angle) beta
    beta_rad = fsolve(theta_beta_mach_relation, beta_guess, 
                     args=(M1, theta_deg, gamma))[0]
    
    # Convert beta from radians to degrees and round it
    return round(np.degrees(beta_rad), 4)

# Key Physical Concepts:
# ----------------------
# 1. Oblique Shock Waves:
#    - Form when supersonic flow encounters a compressive turn
#    - Flow deflects by angle θ through shock at angle β
#    - Downstream flow remains supersonic for weak shocks
#    - Downstream flow becomes subsonic for strong shocks
#
# 2. Weak vs Strong Solutions:
#    - Weak shock: Smaller β, usually preferred physically
#    - Strong shock: Larger β, higher entropy rise
#    - Maximum θ: Boundary between attached and detached shocks
#
# 3. Applications:
#    - Supersonic airfoils and wings
#    - Engine intakes and nozzles
#    - Compression ramps and wedges
#    - Shock-expansion theory


import numpy as np
import matplotlib.pyplot as plt
from Normal_shock import normal_shock, stagnation_temp, stagnation_pressure
from oblique_shock import oblique_shock
from Theta_beta_mach import theta_max
from Isentropic_diffuser import isentropic_diffuser

def analyze_compressible_flow(M1=3, P1=1, T1=273.15, theta1=5, theta2=4):
    """
    Comprehensive analysis of compressible flow through multiple oblique shocks,
    a final normal shock, and an isentropic diffuser.
    
    This function simulates a typical supersonic intake system where:
    1. Multiple oblique shocks compress and decelerate the flow
    2. A terminal normal shock brings flow to subsonic conditions
    3. An isentropic diffuser further decelerates the subsonic flow
    
    Parameters:
    -----------
    M1 : float
        Upstream Mach number (must be > 1 for supersonic flow)
    P1 : float
        Upstream static pressure (normalized units)
    T1 : float
        Upstream static temperature (Kelvin)
    theta1 : float
        First deflection angle in degrees
    theta2 : float
        Additional deflection angle for iterative oblique shocks
    
    Returns:
    --------
    tuple : (results, eff_initial, eff_final)
        - results: List of flow states at each process step
        - eff_initial: Initial stagnation pressure (baseline for efficiency)
        - eff_final: Final stagnation pressure after normal shock
    
    Process Flow:
    -------------
    Supersonic Inflow → Oblique Shocks → Normal Shock → Subsonic Diffuser → Outflow
    """
    # Initial stagnation conditions (reference state)
    T01 = stagnation_temp(M1, T1)
    P01 = stagnation_pressure(M1, P1)
    eff_initial = P01  # Reference stagnation pressure for efficiency calculation

    # Initialize results container
    result = []
    result.append([M1, P1, T1, T01, P01, "Initial State"])

    # First oblique shock - initial compression
    shock1 = oblique_shock(M1, P1, T1, theta1)
    if isinstance(shock1, str):
        raise RuntimeError(f"Error in first oblique shock: {shock1}")

    # Append oblique shock result (M, P, T, T0, P0, label)
    result.append([shock1[0], shock1[1], shock1[2], shock1[3], shock1[5], "Oblique Shock"])

    # Iterate through additional oblique shocks while deflection sum < theta_max
    # This simulates multiple compression ramps in a supersonic intake
    while True:
        # Check if adding theta2 exceeds maximum deflection for current Mach
        try:
            tbm = theta_max(shock1[0])
        except Exception:
            # If theta_max calculation fails, break to avoid infinite loop
            break

        if theta1 + theta2 >= tbm:
            break

        # Apply additional oblique shock
        shock1 = oblique_shock(shock1[0], shock1[1], shock1[2], theta1, theta2)
        if isinstance(shock1, str):
            print(f"Warning: oblique shock iteration stopped: {shock1}")
            break

        result.append([shock1[0], shock1[1], shock1[2], shock1[3], shock1[5], "Oblique Shock"])

    # After oblique shocks, apply terminal normal shock
    # This brings the flow to subsonic conditions for the diffuser
    shock1 = normal_shock(shock1[0], shock1[1], shock1[2])
    if isinstance(shock1, str):
        raise RuntimeError(f"Error in normal shock: {shock1}")

    a = shock1[3]   # Stagnation temperature (conserved across shocks)
    b = shock1[5]   # Stagnation pressure after normal shock
    eff_final = b   # Final stagnation pressure for efficiency calculation
    result.append([shock1[0], shock1[1], shock1[2], a, b, "Normal Shock"])

    # Apply isentropic diffuser to further decelerate subsonic flow
    # This recovers pressure and prepares flow for combustion chamber or compressor
    diffuser_out = isentropic_diffuser(shock1[0], shock1[1], shock1[2])
    if isinstance(diffuser_out, str):
        print(f"Warning: diffuser returned error: {diffuser_out}")
    else:
        # Handle different possible return formats from diffuser function
        try:
            M_d, P_d, T_d = diffuser_out[0], diffuser_out[1], diffuser_out[2]
            T0_d = diffuser_out[3] if len(diffuser_out) > 3 else a
            P0_d = diffuser_out[5] if len(diffuser_out) > 5 else b
            result.append([M_d, P_d, T_d, T0_d, P0_d, "Diffuser"])
        except Exception:
            # Fallback: append last known state but label as diffuser
            result.append([shock1[0], shock1[1], shock1[2], a, b, "Diffuser"])

    return result, eff_initial, eff_final

def print_compact_results(result, eff_initial, eff_final):
    """
    Print a compact, human-readable table showing flow properties at each process step.
    
    Parameters:
    -----------
    result : list
        List of flow states from analyze_compressible_flow
    eff_initial : float
        Initial stagnation pressure reference
    eff_final : float
        Final stagnation pressure for efficiency calculation
    """
    headers = ["Step", "Process", "M", "P (atm)", "T (K)", "T0", "P0"]
    print("\n" + "="*80)
    print("COMPRESSIBLE FLOW ANALYSIS RESULTS")
    print("="*80)
    print("{:<6s} {:<15s} {:>6s} {:>10s} {:>10s} {:>10s} {:>10s}".format(*headers))
    print("-" * 80)

    for i, row in enumerate(result):
        # Expecting row = [M, P, T, T0, P0, proc]
        try:
            M, P, T, T0, P0, proc = row
        except ValueError:
            # In case some row has different length, handle gracefully
            padded = list(row) + [""] * (6 - len(row))
            M, P, T, T0, P0, proc = padded

        print("{:<6d} {:<15s} {:6.3f} {:10.4f} {:10.2f} {:10.2f} {:10.4f}".format(
            i + 1, str(proc)[:15], float(M), float(P), float(T), float(T0), float(P0)
        ))

    # Calculate performance metrics
    efficiency = eff_final / eff_initial if eff_initial != 0 else float('nan')
    P2_P1 = result[-1][1] / result[0][1] if result[0][1] != 0 else float('nan')
    total_temp_ratio = result[-1][3] / result[0][3] if result[0][3] != 0 else float('nan')

    print("-" * 80)
    print("\nPERFORMANCE SUMMARY:")
    print(f" Overall Pressure Recovery (P0_final/P0_initial): {efficiency:.4f}")
    print(f" Total Static Pressure Ratio (P_final/P_initial):  {P2_P1:.4f}")
    print(f" Stagnation Temperature Ratio (T0_final/T0_initial): {total_temp_ratio:.4f}")
    print(f" Number of Process Steps: {len(result)}")
    print("="*80 + "\n")

def plot_results(result):
    """
    Create comprehensive plots showing flow property variations through the system.
    
    Parameters:
    -----------
    result : list
        List of flow states from analyze_compressible_flow
    """
    M = [step[0] for step in result]
    P = [step[1] for step in result]
    T = [step[2] for step in result]
    T0 = [step[3] for step in result]
    P0 = [step[4] for step in result]
    processes = [step[5] for step in result]

    x = np.arange(len(result))

    plt.figure(figsize=(14, 10))
   

    # Plot 1: Stagnation Pressure (Key performance indicator)
    plt.subplot(2, 2, 1)
    plt.plot(x, P0, marker='o', linestyle='-', color='red', linewidth=2, markersize=6)
    plt.xticks(x, processes, rotation=45, ha='right')
    plt.xlabel('Process Step')
    plt.ylabel('Stagnation Pressure (P0)')
    plt.title('Stagnation Pressure Variation\n(Measure of System Efficiency)')
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Mach Number (Shows flow deceleration)
    plt.subplot(2, 2, 2)
    plt.plot(x, M, marker='s', linestyle='-', color='blue', linewidth=2, markersize=6)
    plt.xticks(x, processes, rotation=45, ha='right')
    plt.xlabel('Process Step')
    plt.ylabel('Mach Number')
    plt.title('Mach Number Variation\n(Supersonic to Subsonic Transition)')
    plt.grid(True, alpha=0.3)
    
    # Plot 3: Static Pressure (Shows compression performance)
    plt.subplot(2, 2, 3)
    plt.plot(x, P, marker='^', linestyle='-', color='green', linewidth=2, markersize=6)
    plt.xticks(x, processes, rotation=45, ha='right')
    plt.xlabel('Process Step')
    plt.ylabel('Static Pressure (P)')
    plt.title('Static Pressure Variation\n(Compression Performance)')
    plt.grid(True, alpha=0.3)
    
    # Plot 4: Static Temperature (Shows thermal effects)
    plt.subplot(2, 2, 4)
    plt.plot(x, T, marker='d', linestyle='-', color='orange', linewidth=2, markersize=6)
    plt.xticks(x, processes, rotation=45, ha='right')
    plt.xlabel('Process Step')
    plt.ylabel('Static Temperature (K)')
    plt.title('Static Temperature Variation\n(Thermal Compression Effects)')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

def _get_input(prompt, cast_fn, default):
    """
    Utility function to get user input with default values and error handling.
    
    Parameters:
    -----------
    prompt : str
        Message displayed to user
    cast_fn : function
        Function to convert input (float, int, etc.)
    default : any
        Default value if user presses Enter
    
    Returns:
    --------
    any : Converted user input or default value
    """
    s = input(f"{prompt} [default: {default}]: ").strip()
    if s == "":
        return default
    try:
        return cast_fn(s)
    except Exception:
        print(f"Invalid input '{s}', using default {default}.")
        return default

if __name__ == "__main__":
    """
    Main executable section - provides interactive analysis of compressible flow systems.
    
    Typical Applications:
    - Supersonic aircraft intakes
    - Ramjet/Scramjet engines
    - Supersonic wind tunnel design
    - High-speed compressor analysis
    """
    print("="*60)
    print("COMPRESSIBLE FLOW SYSTEM ANALYZER")
    print("Supersonic Intake & Diffuser Performance Analysis")
    print("="*60)
    
    # Interactive inputs with sensible defaults
    print("\nEnter flow inputs (press Enter to accept default values):")
    M1 = _get_input("Free-stream Mach number M1 (must be > 1)", float, 3.0)
    if M1 <= 1.0:
        print("Warning: M1 <= 1.0 — oblique/normal shock logic expects supersonic flow. Using default M1=3.0")
        M1 = 3.0

    P1 = _get_input("Static pressure P1 (atm)", float, 1.0)
    T1 = _get_input("Static temperature T1 (K)", float, 273.15)
    theta1 = _get_input("Deflection angle theta1 (degrees)", float, 5.0)
    theta2 = _get_input("Additional deflection (degrees) to iterate (theta2)", float, 4.0)

    try:
        # Perform comprehensive flow analysis
        results, eff_init, eff_final = analyze_compressible_flow(
            M1=M1, P1=P1, T1=T1, theta1=theta1, theta2=theta2
        )
        
        # Display results
        print_compact_results(results, eff_init, eff_final)

        # Optional plotting
        do_plot = input("Show detailed plots? (y/N): ").strip().lower()
        if do_plot == 'y':
            plot_results(results)

    except Exception as e:
        print("Error during analysis:", str(e))
        print("Please check input parameters and ensure all required modules are available.")
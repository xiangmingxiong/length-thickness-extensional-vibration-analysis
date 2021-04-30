# length-thickness-extensional-vibration-analysis
Determine complex material constants of piezoelectric bars in length thickness extensional vibrations

Usage:
1. Set the density, length, width, and thickness of the piezoelectric bar in length_thickness_extensional_vibration_analysis.m.
2. Provide experimental admittance measurements in the input file admittance_measurements.txt. The data of frequency f (Hz) are in Column 1. The conductance G (S) and the susceptance B (S) are in Columns 2 and 3.
3. Run length_thickness_extensional_vibration_analysis.m in MATLAB. The optimal material constants and the minimum average relative error between the measured and best fitting admittance and impedance will be returned. The measured and best fitting admittance and impedance will be plotted in a new figure window.

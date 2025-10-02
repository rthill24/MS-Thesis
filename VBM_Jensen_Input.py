import estimate_RAO
import numpy as np

# Example boat parameters
length = 62.0  # meters
beam = 10.9     # meters
draft = 2.4    # meters
block_coefficient = 0.36  
metacentric_height = 2.47  # meters

# Example wave data
wave_data = [4.5, 6.5, 9, 180] #extreme wave height[m], wave period[s], design speed[m/s], and heading[degrees]
wave_periods = [6.5]  # seconds

boat = estimate_RAO.Boat(length, beam, draft, block_coefficient, metacentric_height)
boat_init = estimate_RAO.estimateRAO(boat, wave_data)

omega = (2*np.pi)/wave_data[1]
omegas = [(2*np.pi)/T for T in wave_periods]
VBM_WI = boat_init.estimate_VBM(omega)/1000  # Convert to kNm
#VBMs = boat_init.compute_VBM(omegas)/1000 # Convert to kNm
VBM_D = 0.5 * VBM_WI #Dynamic bending moment is 50% of the wave induced bending moment, from the ratio determined using LR calcs for CPV

print("Vertical Bending Moment (VBM):", VBM_WI)
print("Dynamic Vertical Bending Moment (VBM):", VBM_D)
#print ("Vertical Bending Moments (VBMs) for wave periods", wave_periods, ":", VBMs)
"""
Written by: Nate Clemett, 2023
This code estimates the RAOs for a given boat in a given sea state. 
It also estimates the vertical motions and vertical bending moments for a given mission profile.
"""
import numpy as np
import matplotlib.pyplot as plt

class Boat:
    """This class defines the main characteristics of the boat.
    """
    def __init__(self, length, beam, draft, block_coefficient, metacentric_height):
        """initializes the class with the given boat characteristics.

        Args:
            length (float): length between perpindiculars in meters
            beam (float): waterline beam in meters
            draft (float): draft in meters
            block_coefficient (float): block coefficient of the boat
            metacentric_height (float): transverse metacentric height in meters
        """
        self.length = length
        self.beam = beam
        self.draft = draft
        self.block_coefficient = block_coefficient
        self.metacentric_height = metacentric_height

class estimateRAO:
    """This class uses Jensen's approximation to estimate the RAOs for a given boat in a given sea state.
    """
    def __init__(self, boat, wave_n_boat, x=0):
        """initializes the class with the given boat, sea state and mission characteristics, and longitudinal position.

        Args:
            boat (class): boat class that contasins the main characteristics
            wave_n_boat (array): array that contains the significant wave height, wave period, boat velocity, and boat heading angle
            x (int, optional): longitudinal position. Defaults to 0.
        """
        # Constants
        self.g = 9.81                                                                               # Acceleration due to gravity in m/s^2
        self.rho = 1025                                                                             # Density of sea water in kg/m^3
        
        # Boat data
        self.boat = boat
        self.L = self.boat.length                                                                   # Length between perpindiculars in meters
        self.B0 = self.boat.beam                                                                    # Waterline beam in meters
        self.T = self.boat.draft                                                                    # Draft in meters
        self.Cb = self.boat.block_coefficient
        self.GM_T = self.boat.metacentric_height
        self.B = self.B0*self.Cb
        
        self.delta = self.L*self.B*self.T*self.Cb                                                   # Displacement in m^3
        self.C = 0.373 + 0.023*(self.B/self.T) - 0.043*(self.L/100)                                 # For L > 45m 
        self.T_n = 2*self.C*self.B/self.GM_T                                                        # Natural period of roll in seconds according to IMO A.685(17) resolution   
        self.x = x                                                                                  # Longitudinal position x from the center of gravity
        self.C_44 = self.g*self.delta*self.GM_T                                                     # Restoring Moment Coefficient

        # Ocean data
        self.Hs = wave_n_boat[0]                                                                     # Significant wave height in meters
        self.period = wave_n_boat[1]                                                                 # Wave period in seconds
        self.current_omega = 2*np.pi/self.period                                                    # Wave frequency in radians per second
        
        # Mission data
        self.V = wave_n_boat[2]                                                                     # Boat velocity in m/s
        self.Fn = self.V/np.sqrt(self.g*self.L)                                                     # Froude number                                                           
        #Boat heading should be relative to the wave direction
        self.Beta = wave_n_boat[3]                                                                  # Boat heading angle in degrees relative to the wave direction. Modulo wraps to 360
        self.Beta_rads = np.deg2rad(self.Beta)                                                      # Boat heading angle in radians relative to the wave direction     
        
        
    def update_parameters(self, omega):
        """Updates the parameters for the given boat in the given sea state and wave frequency.

        Args:
            omega (float): wave frequency in radians per second
        """
        self.omega = omega                                                                          # Wave frequency in radians per second
        self.k = omega**2/self.g                                                                  # Wave number in 1/m
        self.k_e = np.abs(self.k*np.cos(self.Beta_rads))                                            # Effective wave number in 1/m
        self.kappa = np.exp(-self.k_e*self.T)                                                       # Smith Correction Factor
        self.alpha = 1 - self.Fn*np.sqrt(self.k*self.L)*np.cos(self.Beta_rads)                      # Froude-Krylov force correction factor
        #self.omega_bar = self.alpha*omega                                                           # Frequency of encounter in radians per second
        self.omega_bar = omega - self.k*self.V*np.cos(self.Beta_rads)                                # Frequency of encounter in radians per second
        #self.A = 2 * np.sin(0.5*self.k*self.Beta_rads*(self.alpha**2)) * np.exp(-self.k*self.T*(self.alpha**2))
        self.A = 2*np.sin((self.omega_bar**2*self.B)/(2*self.g)) * np.exp((-self.omega_bar**2*self.T)/self.g)
        self.f = np.sqrt((1-self.k*self.T)**2 + (self.A**2/(self.k*self.B*self.alpha**3))**2)


    def estimate_vertical_motions(self, omega):
        """Uses Jensen's approximation to estimate the vertical motions for a given boat in a given sea state and wave frequency.

        Args:
            omega (float): wave frequency in radians per second

        Returns:
            three floats: vertical motions in heave, pitch, roll, for the given boat in the given sea state and wave frequency
        """
        self.update_parameters(omega)
        # Forcing functions
        self.F = self.kappa*self.f * (2/(self.k_e*self.L)) * np.sin(0.5*self.k_e*self.L)
        self.G = self.kappa*self.f * (24/((self.k_e*self.L)**2*self.L)) * ( np.sin(0.5*self.k_e*self.L) - 0.5*self.k_e*self.L*np.cos(0.5*self.k_e*self.L)) 
        # Frequency response
        self.eta = (np.sqrt((1 - 2*self.k*self.T*(self.alpha**2))**2 + ((self.A**2)/(self.k*self.B*(self.alpha**2)))**2))**-1                                                                                                    # From Book Assumption
        # Frequency response functions for heave (w) and pitch (theta)
        self.phi_heave = np.abs(self.eta*self.F)
        self.phi_pitch = np.abs(self.eta*self.G)
        self.phi_vert_motion = np.sqrt(self.phi_heave**2 + self.x**2*self.phi_pitch**2)
        self.phi_vert_acc = self.alpha**2*self.k*self.g*self.phi_vert_motion

        return self.phi_heave, self.phi_pitch, self.phi_vert_motion, self.phi_vert_acc
    
    def compute_vertical_motions(self, omegas):
        """Computes vertical motions for a range of omegas.

        Args:
            omegas (list of floats): wave frequencies in radians per second

        Returns:
            three list of floats: vertical motions in heave, pitch, roll, for the given boat in the given sea state and wave frequency
        """
        phi_heave_values = [self.estimate_vertical_motions(omega)[0] for omega in omegas]
        phi_pitch_values = [self.estimate_vertical_motions(omega)[1] for omega in omegas]
        phi_vert_motion_values = [self.estimate_vertical_motions(omega)[2] for omega in omegas]
        phi_vert_acc_values = [self.estimate_vertical_motions(omega)[3] for omega in omegas]
        return phi_heave_values, phi_pitch_values, phi_vert_motion_values, phi_vert_acc_values
    
    def estimate_roll(self, omega):
        """Uses Jensen's approximation to estimate the roll for a given boat in a given sea state and wave frequency.
        Args:
            omega (float): wave frequency in radians per second
        """
        ##### ROLL CURRENTLY NOT WORKING
        self.update_parameters(omega)
        self.phi_roll = np.abs(self.M) / np.sqrt( self.C_44*(-self.omega_bar**2*((self.T_n/(2*np.pi))**2)+1)**2 + self.omega_bar**2*self.B_44**2)
    
    def compute_all(self, omegas):
        """Computes vertical motions, vertical bending moments, and roll for a range of omegas.
        Args:
            omegas (list of floats): wave frequencies in radians per second
        Returns:
            three list of floats: vertical motions in heave, pitch, roll, for the given boat in the given sea state and wave frequency
            list of floats: vertical bending moments for the given boat in the given sea state and wave frequency
        """
        phi_heave, phi_pitch, phi_vert_motion, phi_vert_acc = self.compute_vertical_motions(omegas)
        phi_M = self.compute_VBM(omegas)
        return phi_M, phi_heave, phi_pitch, phi_vert_motion, phi_vert_acc

    
    def estimate_VBM (self,omega):
        """Uses Jensen's approximation to estimate the vertical bending moment for a given boat in a given sea state and wave frequency.
        Args:
            omega (float): wave frequency in radians per second
        Returns:
            float: vertical bending moment for the given boat in the given sea state and wave frequency
        """
        self.update_parameters(omega)
        # Correction factor for the block Coefficient
        Cb = max(0.6, self.Cb)
        gamma = 2.5*(1 - Cb)
        FcCb = (1-gamma)**2 + 0.6*self.alpha*(2-gamma)
        # Speed Correction factor
        FvFn = 1 + 3*self.Fn**2
        self.phi_M = self.kappa*(np.divide((1 - self.k*self.T),(self.k_e*self.L)**2)) * (1 - np.cos(0.5*self.k_e*self.L) - 0.25*self.k_e*self.L*np.sin(0.5*self.k_e*self.L)) * FvFn * FcCb * np.cbrt(np.abs(np.cos(self.Beta_rads))) * self.rho*self.g*self.B0*self.L**2

        return self.phi_M
    
    def compute_VBM(self, omegas):
        """Computes vertical bending moments for a range of omegas.

        Args:
            omegas (list of floats): wave frequencies in radians per second

        Returns:
            list of floats: vertical bending moments for the given boat in the given sea state and wave frequency
        """
        phi_M_values = np.array([self.estimate_VBM(omega) for omega in omegas])
        return phi_M_values
    
    def safe_divide(self, numerator, denominator, default_value=0.0):
        if denominator == 0:
            return default_value
        return numerator / denominator
    
    def VBM_coeff_sign(self, omega):
        """gets the sign of the two VBM coefficients

        Args:
            omega (float): wave frequency in radians per second
        """
        self.update_parameters(omega)
        coeff_1 = np.cos(0.5*self.k_e*self.L)
        coeff_2 = 0.25*self.k_e*self.L*np.sin(0.5*self.k_e*self.L)
        overall = 1 - coeff_1 - coeff_2
        return coeff_1, coeff_2, overall

    
    def plot_data(self, omegas, data, labels, title):
        """Plots the computed data."""
        
        if title == "Vertical Motions Estimate using Jensen's Approximation":
            fig, axes = plt.subplots(len(data), 1, figsize=(10, 15))
            
            for ax, datum, label in zip(axes, data, labels):
                ax.plot(omegas, datum, label=label)
                ax.axvline(x=self.current_omega, color='red', linestyle='--', label="Current $\omega$")
                
                # Only label x-axis for the last subplot
                if label != "Roll":
                    ax.set_xticklabels([])  # Remove x-axis tick labels for non-last subplots
                else:
                    ax.set_xlabel("$\omega$ [rad/sec]")
                
                ax.set_ylabel(f"{label} Estimate")
                ax.legend()
                ax.grid(True)
            
            plt.suptitle(title)
            plt.tight_layout()
            plt.show()
        
        else:
            plt.figure(figsize=(10, 6))
            for datum, label in zip(data, labels):
                plt.plot(omegas, datum, label=label)
            plt.axvline(x=self.current_omega, color='red', linestyle='--', label="Current $\omega$")
            plt.xlabel("$\omega$ [rad/sec]")
            plt.ylabel("Estimate")
            plt.title(title)
            plt.legend()
            plt.grid(True)
            plt.show()

"""This module plots and outputs beamforming pattern
of a circular phased array.

Author: Tom Hirshberg
""" 

import numpy as np
import matplotlib.pyplot as plt
import argparse

from beamforming_pattern_gen import *


# Get amplitude law
def get_amplitude_law(N, law = 'constant', minAmp = 1):
    """Computes an amplitude law given N (number of elements),
    law (type of law) and minAmp (minimum amplitude).
    """
    amp_law = []
    
    for n in range(N):
        if law == 'constant':
            amp_law.append(1)
        elif law == 'linear':
            beta = 0 if N%2!=0 else (1-minAmp)/(N-1)
            amp_law.append((minAmp-1-beta) * 2/(N-1) * abs(n - (N-1) / 2) + 1 + beta)
        elif law == 'log_linear':
            beta = 0 if N%2!=0 else (1-lineartodB(minAmp))/(N-1)
            amp_law.append(dBtoLinear((lineartodB(minAmp)-beta) * 2/(N-1) * abs(n-(N-1)/2) + beta))
        elif law == 'poly2':
            beta = 0 if N%2!=0 else (1-minAmp)/(N-1)
            amp_law.append((minAmp-1-beta) * (2/(N-1))**2 * (n-(N-1)/2)**2 + 1 + beta**2)
        elif law == 'poly3':
            beta = 0 if N%2!=0 else (1-minAmp)/(N-1)
            amp_law.append((minAmp-1-beta) * (2/(N-1))**3 * abs(n-(N-1)/2)**3 + 1 + beta**3)    
    
    return np.array(amp_law)


# Get phase law
def get_circular_phase_law(num_sources, alpha, R, wavelength, phi):
    phase_law = []
    for n in range(num_sources):
        phase_law.append(2 * np.pi  / wavelength * (2 * R * np.sin(n * alpha)) * np.sin(phi + n * alpha))
    return phase_law   
    

def get_circular_pattern(num_sources, R, wavelength, phi, logScale=True, ccw=True):
    alpha = 2 * np.pi / num_sources
    amp_law = np.ones(num_sources)
    phase_law = get_circular_phase_law(num_sources, alpha, R, wavelength, phi)
    
    theta = np.arange(0, 2*np.pi, np.radians(0.1))
    mag = []

    for theta_i in theta:
        im = re = 0
        # Phase shift due to off-boresight angle
        
        # Compute sum of effects of elements
        for n in range(num_sources):
            if ccw:
                factor = n * alpha
            else:
                factor = -n * alpha

            psi = np.sin(theta_i + n * alpha )
            coeff = 2 * np.pi  / wavelength * (2 * R * np.sin(n * alpha))
            # using constant amplitude law for simplicity #TODO: will be optimized in pyRoomAcoustics

            im += np.sin(coeff * (psi + np.sin(phi + factor)))
            re += np.cos(coeff * (psi + np.sin(phi + factor)))
        
        magnitude = np.sqrt(re**2 + im**2)/num_sources

        if logScale:
            magnitude = 20*np.log10(magnitude)
        mag.append(magnitude)
    
    
    return theta, mag, amp_law, phase_law    


def plot_pattern(theta, mag, amp_law, phase_law, polar=False, output_file=None):
    """Plots a magnitude pattern, amplitude law and phase law.
    Optionnally, it can export the pattern to output_file. 
    """
    # Default size & dpi
    plt.figure(figsize=(10,4),dpi=100)
    
    # Plot pattern
    if polar:
        ax = plt.subplot(131, polar=True)
        ax.plot(theta, mag)
        ax.set_theta_zero_location("N")
        ax.set_thetalim(0, 2*np.pi)
    else:
        ax = plt.subplot(131)
        ax.plot(theta,mag)
        ax.grid(True)
        ax.set_xlim([0, 2*np.pi])
        ax.set_xlabel('radinas')
    plt.ylim(top=0)
    plt.title("Antenna pattern - circular array")
    
    # Plot amplitude law
    ax = plt.subplot(132)
    ax.plot(range(len(amp_law)), amp_law, marker='o')
    ax.set_ylim([0,1])
    plt.title("Amplitude law")
    print("Amplitude law:")
    print(amp_law)
    
    # Plot phase law
    ax = plt.subplot(133)
    plt.title("Phase law")
    ax.plot(range(len(phase_law)),np.rad2deg(phase_law), marker='o')
    print("Phase law:")
    print(phase_law)
    
    # Show and save plot
    if output_file is not None:
        plt.savefig(output_file + '.png')
        np.savetxt(output_file + '.txt', np.transpose([theta,mag]),
                   delimiter='\t', header="Angle [deg]\tMagnitude [dB]")
    plt.show()


def main(args):
    # Parameters
    N = args.number_elements # Elements
    c = args.wave_celerity #m/s
    f = args.frequency #Hz
    phi = np.radians(args.steering_angle) #deg
    polar = args.polar #True=polar patter, False=cartesian pattern
    logScale = args.log_scale #True=output in dB, False=linear output
    
    output_file = 'pattern_' + str(N) if args.save_output else None
    wavelength = c/f #m

    theta, mag, amp_law, phase_law = get_circular_pattern(num_sources=N, R=0.1, wavelength=wavelength, phi=phi, logScale=logScale)

    plot_pattern(theta, mag, amp_law,phase_law,polar,output_file=output_file)


if __name__ == '__main__':
    args = get_args()
    main(args)

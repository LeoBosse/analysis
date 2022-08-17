#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
rc('axes.formatter', useoffset=False)
mpl.rcParams['font.size'] = 16
mpl.rcParams['svg.fonttype'] = 'none'

DtoR = np.pi / 180.
RtoD = 1. / DtoR
TAU = 2 * np.pi

def GetPola(V, Vcos, Vsin):
    """Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
    I0 = 2 * V
    DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V * 100
    AoLP = np.arctan2(Vsin, Vcos) / 2
    return abs(I0), abs(DoLP), AoLP

def GetV(signal):
    """return the average of a signal over 1 rotation"""
    integral = sum(signal)
    return integral / N_pts

def GetVcos(signal):
    """Return the average of V*cos(2*theta) over 1 rotation"""
    x = signal * np.cos(2 * filter_angles)
    return sum(x) / N_pts

def GetVsin(signal):
    """Return the average value of -V*sin(2*theta) over 1 rotation"""
    y = - signal * np.sin(2 * filter_angles)
    return sum(y) / N_pts


def GetStokesTime(signal):
    """return the average of a signal over N_rot rotation"""
    V = np.zeros(N_rot)
    Vcos = np.zeros(N_rot)
    Vsin = np.zeros(N_rot)
    for ir in range(N_rot):
        start_rot = ir * N_pts
        end_rot = start_rot + N_pts
        tmp_signal = signal[start_rot:end_rot]

        V[ir] = sum(tmp_signal) / N_pts
        Vcos[ir] = sum(tmp_signal * np.cos(2 * filter_angles[start_rot:end_rot])) / N_pts
        Vsin[ir] = sum(tmp_signal * np.sin(2 * filter_angles[start_rot:end_rot])) / N_pts

    return V, Vcos, Vsin

def FlatSignal(I = 1):
    """Returns a flat signal of intensity I. Numpy array of the shape of the polarising filter angle array."""
    return I * np.ones_like(filter_angles)

def CosSignal(I = 1, D = 1, freq = 2, phase = 0, filter_rot = 0):
    """Returns a cos signal of average intensity I.
    I: Base intensity (=constant). Correspond to the non-polarised light
    D: Relative amplitude of the cos wave with respect to I. Between 0 and 1.
    freq: frequency of the signal. Number of periods during one turn of the polarising filter.
    phase: in radians. Phase shift of the signal
    filter_rot: int value. Default 0. Shift the signal by a multiple of 2π. Usufull for time study of weird frequencies.
    Numpy array of the shape of the polarising filter angle array."""
    S = FlatSignal(I=I/2)
    S += I * D * np.cos(freq * (filter_angles + (filter_rot * 2 * np.pi) - phase)) / 2
    return S

def SimpleGate(width = 10*DtoR, angle = 180*DtoR):
    """Returns a signal egal to 1 only between two fixed values defined by an angle and width. Zero everywhere else.
    width: Width of the gate in radians. The gate will span width/2 on both side of the angle paramter.
    angle: Angle at the center of the gate.

    Returns a numpy array of the shape of the polarising filter angle array."""
    min = (angle - width / 2)
    max = (angle + width / 2)

    S = np.where((min < filter_angles) & (filter_angles < max), 1, 0)
    return S

def GateSignal(I = 1, D = 1, width = 45*DtoR, space = 180*DtoR, phase = 0):
    """Returns a periodic gate signal using superposition of the SimpleGate() signal at given frequency.
    I: Base intensity (=constant). The gate signal will be added to this base.
    D: Relative amplitude of the gate with respect to I. Between 0 and 1.
    width: Width of a gate in radians.
    space: Interval in radians between two gates. ~wavelength of the signal.

    Returns a numpy array of the shape of the polarising filter angle array."""
    S = FlatSignal(I=I/2)

    #In the case I is set to zero, produce a gate with zero background. D is then the absolute value of the gate, not relative to I anymore
    if I == 0:
        I = 1

    nb_gates = int(N_rot * TAU / space) + 2 #Adds 2 gates in case the first one and/or the last fall outside the intervall [0, 2π]

    a = phase - space #Angle of the first gate. (might be < 0)
    for i in range(nb_gates): # Add each gate to the signal
        S += I * D * SimpleGate(width=width, angle=a) / 2
        a += space

    return S

def NeonSignal(I = 1, D = 1, on_time = 5, off_time = 10, flash_time = 2, flash_width = 45*DtoR, flash_space = 180*DtoR, phase_time = 0):
    main_interval = (off_time + on_time + 2 * flash_time) * TAU
    S = GateSignal(I, D, width = on_time * TAU, space = main_interval, phase = phase_time * TAU)

    flash_phase = (- flash_time - on_time/2 + phase_time) * TAU
    for i in [0, 1]:
        for j in range(int(flash_time*TAU / flash_space)):
            S += GateSignal(I=0, D=I*D, width = flash_width, space = main_interval, phase = flash_phase)

            flash_phase += flash_space

        flash_phase += on_time * TAU + flash_space

    return S


def FillParameterSpace(param_name, param_space, function, function_keywords=None):
    """Returns the polarisation neasured for a given signal when varying one of the signal parameter.
    param_name: name of the paramter to vary. Must correspond to the input parameter used by the function.
    param_space: list of all values taken by the varying parameter
    function: function used to create the signal. For exemple GateSignal or CosSignal.
    function_keywords: Dictionnary of the parameters to pass to the signal function. If the varying parameter is given, it will be overwritten.

    Returns the parameter space, and the intensity, DoLP, AoLP corresponding to each value.
    """

    if not function_keywords: #Create an enpty dictionnary if the function_keywords parameter is not given
        function_keywords = {param_name:None}

    #initialize empty list for polarisation values to return
    I_list = np.empty_like(param_space)
    DoLP_list = np.empty_like(param_space)
    AoLP_list = np.empty_like(param_space)

    #Loop over all values of the varying parameter and compute the corresponding signal and polarisation
    for ip, p in enumerate(param_space):
        #Update the paramter dictionnary with the new parameter space value
        function_keywords[param_name] = p

        # Create the corresponding signal
        signal = function(**function_keywords)
        # print(signal)

        #Compute Stokes paramters equivalent for lock0in amplifier method
        V    = GetV(signal)
        Vcos = GetVcos(signal)
        Vsin = GetVsin(signal)

        # Compute polarisation values (I, DoLP, AoLP)
        I, DoLP, AoLP = GetPola(V, Vcos, Vsin)
        I_list[ip] = I
        DoLP_list[ip] = DoLP
        AoLP_list[ip] = AoLP

    return param_space, I_list, DoLP_list, AoLP_list


N_pts = 1000  #Number of points for one turn of the polarising filter
N_rot = 300
filter_angles = np.linspace(0, N_rot * 2*np.pi, N_rot * N_pts) #List of angles (rad) for the polarising filter between 0 and 2π.

x_axis = None

# plt.plot(filter_angles*RtoD, NeonSignal(I = 1, D = 1))
# plt.plot(filter_angles, NeonSignal(I = 1, D = 1))
# plt.show()

##########
# Below are different study of signals with varying paramters.
# (Un)comment the lines for quick study
##########

### NeonSignal
fig_title = "Polarisation for a perfect neon signal in time"
x_axis_label = 'Time (Rotation nb of the pola filter)'
x_axis = range(N_rot)
# signal = NeonSignal(I = 1, D = 1, on_time = 9.67, off_time = 11.33, flash_time = 4.7, flash_width = 45*DtoR, flash_space = 130*DtoR, phase_time = 15)
signal = NeonSignal(I = 1, D = 1, on_time = 10.25, off_time = 10.25, flash_time = 10, flash_width = 45*DtoR, flash_space = 180*DtoR, phase_time = 15.3)
V, Vcos, Vsin = GetStokesTime(signal)
I_list, DoLP_list, AoLP_list = GetPola(V, Vcos, Vsin)
plt.plot(filter_angles/TAU, signal)

### Cos frequency variations
# fig_title = "Polarisation for a cos signal of varying frequency with 100% DoLP"
# x_axis_label = 'Frequency (per filter rotation)'
# parameter_space, I_list, DoLP_list, AoLP_list = FillParameterSpace('freq', np.linspace(0.01, 10, 1000), CosSignal, {'I':1, 'D':1, 'phase':0})

### Cos time variations
# frequency = 3.6
# fig_title = f"Polarisation time variations for a {frequency} Hz cos signal with 100% DoLP"
# x_axis_label = 'Filter rotation number'
# # phase_shift = lambda f, p0, N: (p0 + N * 2 * np.pi * abs(1/f - 1)) % (2 * np.pi)
# parameter_space, I_list, DoLP_list, AoLP_list = FillParameterSpace('filter_rot', np.arange(0, 100), CosSignal, {'I':1, 'I':10, 'D':1, 'freq':frequency})

### Cos DoLP variations
# fig_title = "Polarisation for a 2Hz cos signal with varying amplitude"
# x_axis_label = 'Relative Amplitude'
# parameter_space, I_list, DoLP_list, AoLP_list = FillParameterSpace('D', np.linspace(0.01, 1, 1000), CosSignal)

### Gate 2Hz width variation
# fig_title = "Polarisation for a 2Hz gate signal of 100% amplitude with varying width"
# x_axis_label = 'Gate width (rad)'
# parameter_space, I_list, DoLP_list, AoLP_list = FillParameterSpace('width', np.linspace(0.01*np.pi, np.pi, 100), GateSignal, {'I': 1, 'D':1, 'space':180*DtoR, 'phase':0})

### Gate freq variation with a 45 deg width
# fig_title = "Polarisation for a 45 degree width gate signal of 100% amplitude with varying interval"
# x_axis_label = 'Interval (degree)'
# parameter_space, I_list, DoLP_list, AoLP_list = FillParameterSpace('space', np.linspace(45*DtoR, 450*DtoR, 1000), GateSignal, {'I': 1, 'D':1, 'width':45*DtoR, 'phase':0})
# # x_axis = np.log10(parameter_space)
# x_axis = parameter_space * RtoD


##########
# Plotting procedure.
# Uses the title and x-axis defined above to adapt for each paramter study
##########


n_plots = 3
f, axs = plt.subplots(n_plots, gridspec_kw = {'hspace':0}, sharex = True)
# axs.plot(filter_angles*RtoD, signal)

f.suptitle(fig_title)

if x_axis is None:
    x_axis = parameter_space

ip = 0
axs[ip].plot(x_axis, I_list)
axs[ip].set_ylabel("Intensity (AU)")
ip += 1
axs[ip].plot(x_axis, DoLP_list)
axs[ip].set_ylabel("DoLP (%)")
ip += 1
axs[ip].plot(x_axis, AoLP_list*RtoD)
axs[ip].set_ylabel("AoLP (deg)")

axs[ip].set_xlabel(x_axis_label)


plt.show()



# this is a python script for quantum photonics PIC simulation
import numpy as np
import matplotlib.pyplot as plt
import sax
from simphony.libraries import siepic, sipann
from simphony.classical import ClassicalSim
from functools import partial
import os

# first we will have a function that defines a mzi 2x2 coupler
# and then attach a phase shifter unit (ideal) to one of the arms

def pcpl(wl=1.55, phase=0) -> sax.SDict:
    """
    create a directional coupler with a phase shifter
    wl - wavelength in um
    phase - phase shift in radians (in units of 2*pi)
    e.g phase = 0.5 means pi phase shift
    """
    cir, info = sax.circuit(
            netlist={
                "instances": {
                    "stdcpl1":"standard_coupler",
                    "waveg_t": "waveguide",
                    "waveg_b": "waveguide",
                    },
                "connections": {
                    "stdcpl1,o2": "waveg_t,o0",
                    "stdcpl1,o3": "waveg_b,o0",
                    },
                "ports": {
                    "in_t":"stdcpl1,o0",
                    "in_b":"stdcpl1,o1",
                    "out_t":"waveg_t,o1",
                    "out_b":"waveg_b,o1",
                    }
                },
            models={
                "waveguide": partial(siepic.waveguide, width=500, length=10000),
                "standard_coupler":partial(sipann.standard_coupler, length=6363, horizontal=5000, vertical=2500),
                }
                )

    # need to first recover the effective index
    # Return the composite model.
    #deln = 0.000155*phase*2*np.pi*10000/wl
    deln = phase*2*np.pi
    val = cir(wl = wl)
    #print((np.cos(deln) + 1j*np.sin(deln)))
    val[('in_t', 'out_t')] = val[('in_t', 'out_t')] * (np.cos(deln) + 1j*np.sin(deln))
    val[('out_t', 'in_t')] = val[('out_t', 'in_t')] * (np.cos(deln) + 1j*np.sin(deln))
    val[('out_t', 'in_b')] = val[('out_t', 'in_b')] * (np.cos(deln) + 1j*np.sin(deln))
    val[('in_b', 'out_t')] = val[('in_b', 'out_t')] * (np.cos(deln) + 1j*np.sin(deln))
    return val



def trr(wl=1.55, theta=0, phi=0) -> sax.SDict:
    """
    create a directional coupler with a phase shifter
    wl - wavelength in um
    phase - phase shift in radians (in units of 2*pi)
    e.g phase = 0.5 means pi phase shift
    """
    su2_cir, info = sax.circuit(
            netlist = {
                "instances": {
                    "pc1": "pcpl",
                    "pc2": "pcpl",
                    },
                "connections": {
                    "pc1,out_t": "pc2,in_t",
                    "pc1,out_b": "pc2,in_b",
                    },
                "ports": {
                    "int": "pc1,in_t",
                    "inb": "pc1,in_b",
                    "outt": "pc2,out_t",
                    "outb": "pc2,out_b",
                    }},
            models={
                "pcpl": partial(pcpl, wl=wl),
                   })
    val = su2_cir(wl=wl, pc1 = {"phase":theta}, pc2 = {"phase":phi})
    return val


import numpy as np
from random import choice
from lorentz import lorentz

def create_anti(particle):
    return [1 - x for x in particle[:-1]] + [1]

def create_duplet(meson, mass_parton=0.01):
    """ Takes position, momentum, and mass from meson
    and creates parton antiparton duplet with same
    energy in LRF and same momentum of CF.
    Standard mass of parton is 0.01 GeV, colour chosen at random """

    # position for both partons is same as meson's position in cf
    pos = meson[0:3] + [1]

    # LRF Calculation of Energy
    mass_meson = meson[-1]

    r = np.sqrt((mass_meson**2 * 0.25 - mass_parton**2))
    phi = np.random.rand()*2*np.pi
    theta = np.random.rand()*2*np.pi
    

    p_parton1 = np.array([0,0,0,0], dtype=np.float32)
    p_parton1[0] = mass_meson * 0.5
    p_parton1[1] = r * np.sin(phi) * np.cos(theta)
    p_parton1[2] = r * np.sin(phi) * np.sin(theta)
    p_parton1[3] = r * np.cos(phi)

    p_parton2 = np.array([0,0,0,0], dtype=np.float32)
    p_parton2[0] = mass_meson * 0.5
    p_parton2[1] = -r * np.sin(phi) * np.cos(theta)
    p_parton2[2] = -r * np.sin(phi) * np.sin(theta)
    p_parton2[3] = -r * np.cos(phi)

    # CF Calculation of Momentum

    v_meson = np.array(meson[4:7]) / meson[3]
    p_meson = np.array(meson[3:7])


    p_parton1 = lorentz(-v_meson, p_parton1)
    p_parton2 = lorentz(-v_meson, p_parton2)
    
    # Chose color at random

    c_parton1 = choice([x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0] ] ])
    c_parton2 = create_anti(c_parton1)

    # return each parton as np.array

    parton1 = [pos, p_parton1, mass_parton, c_parton1]
    parton2 = [pos, p_parton2, mass_parton, c_parton2]

    return parton1, parton2

    

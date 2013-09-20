import numpy as np
from random import choice

def lorentz(beta, fourVector, EPS = 1e-5):
    """ Takes as input relative beta from O -> O'
    and calculates fourVector -> fourVector' """

    beta2 = np.dot(beta, beta)
    if (beta2 <= EPS):
        return fourVector
        
    gamma = 1./np.sqrt(1. - beta2)
    
    Lambda = np.array([ [gamma, -gamma*beta[0], -gamma*beta[1], -gamma*beta[2]],
                     [-gamma * beta[0], 1 + (gamma-1)*beta[0]**2/beta2, (gamma-1)*beta[0]*beta[1]/beta2, (gamma-1)*beta[0]*beta[2]/beta2],
                     [-gamma * beta[1], (gamma-1)*beta[1]*beta[0]/beta2, 1 + (gamma-1)*beta[1]**2/beta2, (gamma-1)*beta[1]*beta[2]/beta2],
                     [-gamma * beta[2], (gamma-1)*beta[2]*beta[0]/beta2, (gamma-1)*beta[2]*beta[1]/beta2, 1 + (gamma-1)*beta[2]**2/beta2]], dtype=np.float32)

    return np.dot(Lambda, fourVector)

def create_anti(particle):
    return [1 - x for x in particle[:-1]] + [1]

def create_duplet(meson, mass_parton=0.01):
    """ Takes position, momentum, and mass from meson
    and creates parton antiparton duplet with same
    energy in LRF and same momentum of CF.
    Standard mass of parton is 0.01 GeV, colour chosen at random """

    # position for both partons is same as meson's position in cf
    pos = np.array(meson[0:3] + [0], dtype=np.float32)

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

    p_parton1 = lorentz(-v_meson, p_parton1)
    p_parton2 = lorentz(-v_meson, p_parton2)
    
    p_parton1 = np.hstack((p_parton1[1:4], [mass_parton])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton2 = np.hstack((p_parton2[1:4], [mass_parton]))
    # Chose color at random

    c_parton1 = choice([x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0] ] ])
    c_parton2 = create_anti(c_parton1)
    c_parton1 = np.array(c_parton1, dtype=np.float32)
    c_parton2 = np.array(c_parton2, dtype=np.float32)

    # return each parton as np.array


    parton1 = [pos, p_parton1, c_parton1]
    parton2 = [pos, p_parton2, c_parton2]

    return parton1, parton2

def create_triplet(baryon, mass_parton=0.01):
    """ Takes position, momentum, and mass from baryon
    and creates parton triplet with same
    energy in LRF and same momentum of CF.
    Standard mass of parton is 0.01 GeV"""

    # position for partons is same as baryon's position in cf
    pos = np.array(baryon[0:3] + [0], dtype=np.float32)


    # LRF Calculation of Energy
    mass_baryon = baryon[-1]

    r = np.sqrt(((mass_baryon/3.)**2 - mass_parton**2))
    phi = np.random.rand()*2*np.pi
    theta = np.random.rand()*2*np.pi
    
    p_parton1 = np.array([0,0,0,0], dtype=np.float32)
    p_parton1[0] = mass_baryon/3.
    p_parton1[1] = r * np.sin(phi) * np.cos(theta)
    p_parton1[2] = r * np.sin(phi) * np.sin(theta)
    p_parton1[3] = r * np.cos(phi)

    r = np.sqrt(((mass_baryon/3.)**2 - mass_parton**2))
    phi = np.random.rand()*2*np.pi
    theta = np.random.rand()*2*np.pi
   
    p_parton2 = np.array([0,0,0,0], dtype=np.float32)
    p_parton2[0] = mass_baryon/3.
    p_parton2[1] = r * np.sin(phi) * np.cos(theta)
    p_parton2[2] = r * np.sin(phi) * np.sin(theta)
    p_parton2[3] = r * np.cos(phi)

    p_parton3 = -p_parton1 - p_parton2
    p_parton3[0] = mass_baryon/3. # this is wrong!

    v_baryon = np.array(baryon[4:7]) / baryon[3]

    p_parton1 = lorentz(-v_baryon, p_parton1)
    p_parton2 = lorentz(-v_baryon, p_parton2)
    p_parton3 = lorentz(-v_baryon, p_parton3)


    p_parton1 = np.hstack((p_parton1[1:4], [mass_parton])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton2 = np.hstack((p_parton2[1:4], [mass_parton])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton3 = np.hstack((p_parton3[1:4], [mass_parton])) # (E, px, py, pz) |-> (px, py, pz, m)

    # set color r,g,b (momenta at random, though no bias introduced)

    c_parton1 = np.array([0,0,1,1], dtype=np.float32)
    c_parton2 = np.array([0,1,0,1], dtype=np.float32)
    c_parton3 = np.array([1,0,0,1], dtype=np.float32)

    # return each parton as np.array

    parton1 = [pos, p_parton1, c_parton1]
    parton2 = [pos, p_parton2, c_parton2]
    parton3 = [pos, p_parton3, c_parton3]
    
    return parton1, parton2, parton3

    

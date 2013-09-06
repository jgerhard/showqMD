import numpy as np
from random import choice

def lorentz(beta, fourVector):
    """ Takes as input relative beta from O -> O'
    and calculates fourVector -> fourVector' """

    beta2 = np.dot(beta, beta)
    gamma = 1./np.sqrt(1. - beta2)
    
    Lambda = np.array([ [gamma, -gamma*beta[0], -gamma*beta[1], -gamma*beta[2]],
                     [-gamma * beta[0], 1 + (gamma-1)*beta[0]**2/beta2, (gamma-1)*beta[0]*beta[1]/beta2, (gamma-1)*beta[0]*beta[2]/beta2],
                     [-gamma * beta[1], (gamma-1)*beta[1]*beta[0]/beta2, 1 + (gamma-1)*beta[1]**2/beta2, (gamma-1)*beta[1]*beta[2]/beta2],
                     [-gamma * beta[2], (gamma-1)*beta[2]*beta[0]/beta2, (gamma-1)*beta[2]*beta[1]/beta2, 1 + (gamma-1)*beta[2]**2/beta2]], dtype=np.float32)

    return np.dot(Lambda, fourVector)

def create_meson(partonA, partonB):
    """ Takes position, velocity, and mass from two partons
    and creates meson. Potential energy between partons is converted
    to mass of meson, cumulative momentum in CF is momentum of meson  """

    # position of meson is in middle of both partons
    pos = (partonA[0:3] + partonB[0:3]) * 0.5

    # Calculation of parton 4-momenta
    velA = partonA[3:7]         # this is now (v_x, v_y, v_z, m_0)
    gammaA = 1./np.sqrt(1-(velA[0]**2 + velA[1]**2 + velA[2]**2))
    pxA = gammaA * velA[-1] * velA[0]
    pyA = gammaA * velA[-1] * velA[1]
    pzA = gammaA * velA[-1] * velA[2]
    momA = np.array(np.sqrt(vel[-1]**2 + px**2 + py**2 + pz**2), px, py, pz)

    velB = partonB[3:7]         # this is now (v_x, v_y, v_z, m_0)
    gammaB = 1./np.sqrt(1-(velB[0]**2 + velB[1]**2 + velB[2]**2))
    pxB = gammaB * velB[-1] * velB[0]
    pyB = gammaB * velB[-1] * velB[1]
    pzB = gammaB * velB[-1] * velB[2]
    momB = np.array(np.sqrt(vel[-1]**2 + px**2 + py**2 + pz**2), px, py, pz)


    

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
    pos = np.array(baryon[0:3] + [1], dtype=np.float32)


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
    p_parton3[0] = mass_baryon/3.

    # CF Calculation of Momentum

    v_baryon = np.array(baryon[4:7]) / baryon[3]
    p_baryon = np.array(baryon[3:7])

    p_parton1 = lorentz(-v_baryon, p_parton1)
    p_parton2 = lorentz(-v_baryon, p_parton2)
    p_parton3 = lorentz(-v_baryon, p_parton3)

    # set color r,g,b (momenta at random, though no bias introduced)

    c_parton1 = np.array([0,0,1,1], dtype=np.float32)
    c_parton2 = np.array([0,1,0,1], dtype=np.float32)
    c_parton3 = np.array([1,0,0,1], dtype=np.float32)

    # return each parton as np.array

    parton1 = [pos, p_parton1, c_parton1]
    parton2 = [pos, p_parton2, c_parton2]
    parton3 = [pos, p_parton3, c_parton3]

    return parton1, parton2, parton3

    

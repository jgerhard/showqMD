import numpy as np
from random import choice
from create_partons import lorentz

def create_meson(partonA, partonB, kappa=0.87):
    """ Takes position, 3momentum, and mass from two partons
    and creates meson. Potential energy between partons is converted
    to mass of meson, cumulative momentum in CF is momentum of meson  """

    # position of meson is in middle of both partons
    posA = partonA[0:3]
    posA = np.hstack(([0.0], posA)) # time component of 4-position not necessary, random number
    posB = partonB[0:3]
    posB = np.hstack(([0.0], posB))

    # Calculation of parton 4-momenta in CF
    momA = partonA[4:7]         
    massA = partonA[7] 
    EnA = np.sqrt(np.dot(momA, momA) + massA**2)
    momA = np.hstack((EnA, momA)) # valid 4-momentum
    
    momB = partonB[4:7]
    massB = partonB[7]
    EnB = np.sqrt(np.dot(momB, momB) + massB**2)
    momB = np.hstack((EnB, momB))
    
    frame_vel = ((momA + momB) / (EnA + EnB))[1:] # beta for CF -> LRF
    
    # Calculation of Hadron Energy in CMF

    lrf_momA = lorentz(frame_vel, momA)
    lrf_momB = lorentz(frame_vel, momB)

    E = (lrf_momA + lrf_momB)[0] # this is sqrt(p^2 + m^2)

    lrf_delta_pos = (lorentz(frame_vel, posA-posB))[1:] # ignore time distance (as is random)

    distance = np.sqrt(np.dot(lrf_delta_pos, lrf_delta_pos))
    E_pot = kappa * distance

    lrf_meson_mom = np.array([E + E_pot, 0, 0, 0]) # this is the momentum of the meson in LRF
    # Boosting back to CF
    meson_mom = lorentz(-frame_vel, lrf_meson_mom)
    meson_pos = 0.5 * (posA + posB)

    return np.hstack((meson_pos[1:], meson_mom))


def create_baryon(partonA, partonB, partonC, kappa=0.87):
    """ Takes position, 3momentum, and mass from three partons
    and creates baryon. Potential energy between partons is converted
    to mass of baryon, cumulative momentum in CF is momentum of baryon  """

    # position of baryon is in middle 
    posA = partonA[0:3]
    posA = np.hstack(([0.0], posA)) # time component of 4-position not necessary, random number
    posB = partonB[0:3]
    posB = np.hstack(([0.0], posB))
    posC = partonC[0:3]
    posC = np.hstack(([0.0], posC))

    # Calculation of parton 4-momenta in CF
    momA = partonA[4:7]         
    massA = partonA[7] 
    EnA = np.sqrt(np.dot(momA, momA) + massA**2)
    momA = np.hstack((EnA, momA)) # valid 4-momentum
    
    momB = partonB[4:7]
    massB = partonB[7]
    EnB = np.sqrt(np.dot(momB, momB) + massB**2)
    momB = np.hstack((EnB, momB)) # valid 4-momentum

    momC = partonC[4:7]
    massC = partonC[7]
    EnC = np.sqrt(np.dot(momC, momC) + massC**2)
    momC = np.hstack((EnC, momC)) # valid 4-momentum
    
    frame_vel = ((momA + momB + momC) / (EnA + EnB + EnC))[1:] # beta for CF -> LRF
    
    # Calculation of Hadron Energy in LRF

    lrf_momA = lorentz(frame_vel, momA)
    lrf_momB = lorentz(frame_vel, momB)
    lrf_momC = lorentz(frame_vel, momC)

    E = (lrf_momA + lrf_momB + lrf_momC)[0] # this is sqrt(p^2 + m^2)
    
    lrf_delta_pos_ab = lorentz(frame_vel, posA-posB)[1:] # ignore time distance (as is random)
    lrf_delta_pos_bc = lorentz(frame_vel, posB-posC)[1:] # ignore time distance (as is random)
    lrf_delta_pos_ca = lorentz(frame_vel, posC-posA)[1:] # ignore time distance (as is random)
    dist_ab = np.sqrt(np.dot(lrf_delta_pos_ab, lrf_delta_pos_ab))
    dist_bc = np.sqrt(np.dot(lrf_delta_pos_bc, lrf_delta_pos_bc))
    dist_ca = np.sqrt(np.dot(lrf_delta_pos_ca, lrf_delta_pos_ca))
    E_pot = 0.5 * kappa * (dist_ab + dist_bc + dist_ca) # Force between different colours is half of force between q q\bar

    lrf_baryon_mom = np.array([E + E_pot, 0, 0, 0]) # this is the momentum of the baryon in LRF

    # Boosting back to CF
    baryon_mom = lorentz(-frame_vel, lrf_baryon_mom)
    baryon_pos = 1./3. * (posA[1:] + posB[1:] + posC[1:])
    return np.hstack((baryon_pos, baryon_mom))

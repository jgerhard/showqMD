import numpy as np
from numpy import sin, cos, sqrt
from random import choice, shuffle

u, d, s = 0.0023, 0.0048, 0.095 # quark masses in GeV

MesonDict = {(-1,101):(d,u), (1,101):(u,d), (-1,104):(d,u), (1,104):(u,d), (-1,106):(s,u), (1,106):(u,s),
             (-1,-101):(u,d), (1,-101):(d,u), (-1,-104):(u,d), (1,-104):(d,u), (-1,-106):(u,s), (1,-106):(s,u)} # antimatter

DeltaDict = {-1:(d,d,d), 0:(d,d,u), 1 : (u,u,d), 2:(u,u,u)}
SigmaDict = {-1:(d,d,s), 0:(u,d,s), 1:(u,u,s)}

def Mass(fMom):
    return np.sqrt(fMom[0]**2 - fMom[1]**2 -fMom[2]**2 -fMom[3]**2)

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

def mesons_partons(meson):
    ityp, charge = meson[-3], meson[-1]
    try:
        return MesonDict[(ityp, charge)]
    except KeyError:
        if abs(ityp) >= 106:    # strangeness case
            return choice([(d,s), (s,d)])
        else:                   # no strangeness
            return choice([(u,u), (d,d)])
            
        

def baryons_partons(baryon):
    ityp, charge = baryon[-3], baryon[-1]
    
    if ityp in range(-26,27):
        (parton1, parton2, parton3) = DeltaDict[charge]
    else:
        (parton1, parton2, parton3) = SigmaDict[charge]

    return (parton1, parton2, parton3)
    # mean = 1./3. * (parton1 + parton2 + parton3)
    # return (mean, mean, mean)


def create_duplet(meson):
    """ Takes position, momentum, mass, itype, isospin, and charge from meson
    and creates parton antiparton duplet with same
    energy in LRF and same momentum of CF.
    Masses of partons are chosen according to itype and charge,
    colour chosen at random """
    
    # position for both partons is same as meson's position in cf
    pos = np.array(meson[0:3] + [0], dtype=np.float32)

    # LRF Calculation of Energy
    mass_meson = meson[-4]
    (mass_parton1, mass_parton2) = mesons_partons(meson)

    bias = (mass_meson**2 + mass_parton1**2 - mass_parton2**2)/(2*mass_meson**2) 
    # part of meson energy for different partons, such that the sum of both 3-momenta 
    # is zero and the energy is exactly the meson energy at rest
    
    phi = np.random.rand()*2*np.pi # random direction for momenta
    theta = np.random.rand()*2*np.pi
    
    r = np.sqrt(bias**2 * mass_meson**2 - mass_parton1**2)
    p_parton1 = np.array([0,0,0,0], dtype=np.float32)
    p_parton1[0] = bias * mass_meson
    p_parton1[1] = r * np.sin(phi) * np.cos(theta)
    p_parton1[2] = r * np.sin(phi) * np.sin(theta)
    p_parton1[3] = r * np.cos(phi)


    # r = np.sqrt((1-bias)**2 * mass_meson**2 - mass_parton2**2) #additional calculation not necessary, is exactly the same as above
    p_parton2 = np.array([0,0,0,0], dtype=np.float32)
    p_parton2[0] = (1-bias)*mass_meson
    p_parton2[1] = -r * np.sin(phi) * np.cos(theta)
    p_parton2[2] = -r * np.sin(phi) * np.sin(theta)
    p_parton2[3] = -r * np.cos(phi)

    # CF Calculation of Momentum

    v_meson = np.array(meson[4:7]) / meson[3]

    p_parton1 = lorentz(-v_meson, p_parton1)
    p_parton2 = lorentz(-v_meson, p_parton2)
    
    p_parton1 = np.hstack((p_parton1[1:4], [mass_parton1])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton2 = np.hstack((p_parton2[1:4], [mass_parton2]))

    # Chose color at random
    c_parton1 = choice([x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0], [1,1,0], [1,0,1], [0,1,1] ] ])
    c_parton2 = create_anti(c_parton1)
    c_parton1 = np.array(c_parton1, dtype=np.float32)
    c_parton2 = np.array(c_parton2, dtype=np.float32)

    # return each parton as np.array


    parton1 = [pos, p_parton1, c_parton1]
    parton2 = [pos, p_parton2, c_parton2]

    return parton1, parton2

def create_triplet(baryon):
    """ Takes position, momentum, mass, itype, isospin, and charge from baryon
    and creates parton triplet with same energy in LRF and same momentum of CF.
    Masses of partons are chosen according to itype and charge,
    colour chosen at random """


    # position for partons is same as baryon's position in cf
    pos = np.array(baryon[0:3] + [0], dtype=np.float32)


    # LRF calculation of 4-momentum of partons
    mass_baryon = baryon[-4]
    (mass_parton1, mass_parton2, mass_parton3) = baryons_partons(baryon)

    u = sqrt(1./9. * mass_baryon**2 - mass_parton1**2)
    d = sqrt(1./9. * mass_baryon**2 - mass_parton2**2)
    s = sqrt(1./9. * mass_baryon**2 - mass_parton3**2)
    # first calculate 2d 3-momentum vectors with sum zero and 
    # length u,d,s (kinetic energy of mass defect for each parton)
    # vectors are: pu = [ 0,  u]
    #              pd = [ x, y-u]
    #              ps = [-x,  -y]

    y = (s**2 + u**2 - d**2)/(2*u)
    x = sqrt(s**2 - y**2)

    # then do an arbitrary 3d rotation by phi, psi and theta.
    # embedd pu, pd, ps via [a,b] |-> [a,b,0]
    # this should provide isotropy!

    phi = np.random.rand()*2*np.pi 
    psi = np.random.rand()*2*np.pi 
    theta = np.random.rand()*2*np.pi

    # [0,u, 0] |-> [u * (cos(phi)sin(psi)+sin(phi)sin(theta)cos(psi)),...]
    p_parton1 = np.array([0,0,0,0], dtype=np.float32)
    p_parton1[0] = mass_baryon/3.
    p_parton1[1] = u * (cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi))
    p_parton1[2] = u * (cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi))
    p_parton1[3] = u * (-sin(phi)*cos(theta)) 


    # [x,y-u, 0] |-> ...
    p_parton2 = np.array([0,0,0,0], dtype=np.float32)
    p_parton2[0] = mass_baryon/3.
    p_parton2[1] = x * ( cos(theta)*cos(psi)) + (y-u) * (cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi))
    p_parton2[2] = x * (-cos(theta)*sin(psi)) + (y-u) * (cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi))
    p_parton2[3] = x * sin(theta) + (y-u) * (-sin(phi)*cos(theta)) 


    # [-x, -y, 0] |-> ...
    p_parton3 = np.array([0,0,0,0], dtype=np.float32) 
    p_parton3[0] = mass_baryon/3.
    p_parton3[1] = (-x) * ( cos(theta)*cos(psi)) + (-y) * (cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi))
    p_parton3[2] = (-x) * (-cos(theta)*sin(psi)) + (-y) * (cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi))
    p_parton3[3] = (-x) * sin(theta) + (-y) * (-sin(phi)*cos(theta)) 


    # boost into computational frame
    v_baryon = np.array(baryon[4:7]) / baryon[3]

    p_parton1 = lorentz(-v_baryon, p_parton1)
    p_parton2 = lorentz(-v_baryon, p_parton2)
    p_parton3 = lorentz(-v_baryon, p_parton3)

    p_parton1 = np.hstack((p_parton1[1:4], [mass_parton1])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton2 = np.hstack((p_parton2[1:4], [mass_parton2])) # (E, px, py, pz) |-> (px, py, pz, m)
    p_parton3 = np.hstack((p_parton3[1:4], [mass_parton3])) # (E, px, py, pz) |-> (px, py, pz, m)

    # set color 

    colors = [x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0] ]]
    shuffle(colors)             # isotropy in colour-space
    c_parton1 = colors[0]
    c_parton2 = colors[1]
    c_parton3 = colors[2]
    if (baryon[-3] < 0):        # anti particles
        c_parton1 = create_anti(c_parton1)
        c_parton2 = create_anti(c_parton2)
        c_parton3 = create_anti(c_parton3)

    # return each parton as np.array
    parton1 = [pos, p_parton1, c_parton1]
    parton2 = [pos, p_parton2, c_parton2]
    parton3 = [pos, p_parton3, c_parton3]
    
    return parton1, parton2, parton3

    
def old_lrf_calculation(baryon, mass_parton=0.01):
    # LRF calculation of 4-momentum of partons
    mass_baryon = baryon[-4]

    r = np.sqrt(((mass_baryon/3.)**2 - mass_parton**2))
    phi = np.random.rand()*2*np.pi # offset
    theta = np.random.rand()*2*np.pi

    alpha = phi 

    p_parton1 = np.array([0,0,0,0], dtype=np.float32)
    p_parton1[0] = mass_baryon/3.
    p_parton1[1] = r * np.sin(alpha) * np.cos(theta)
    p_parton1[2] = r * np.sin(alpha) * np.sin(theta)
    p_parton1[3] = r * np.cos(alpha)

    alpha += 2./3. *  np.pi 

    p_parton2 = np.array([0,0,0,0], dtype=np.float32)
    p_parton2[0] = mass_baryon/3.
    p_parton2[1] = r * np.sin(alpha) * np.cos(theta)
    p_parton2[2] = r * np.sin(alpha) * np.sin(theta)
    p_parton2[3] = r * np.cos(alpha)

    alpha += 2./3. * np.pi

    p_parton3 = np.array([0,0,0,0], dtype=np.float32) 
    p_parton3[0] = mass_baryon/3.
    p_parton3[1] = r * np.sin(alpha) * np.cos(theta)
    p_parton3[2] = r * np.sin(alpha) * np.sin(theta)
    p_parton3[3] = r * np.cos(alpha)


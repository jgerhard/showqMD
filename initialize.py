from OpenGL.GL import *
import random
import numpy as np
from UrQMDfilter import readfile, separate_hadrons, make_partonlist
from random import sample

def fountain_urqmd(maxnum , filename="test.f14", eventnumber=0, parton_mass = 0.01):
    """ Create partons from UrQMD inputfile "filename"
    if number exeeds maxnum only pass maxnum particles """

    events = readfile("test.f14",header=False)
    baryons, mesons = separate_hadrons(events[eventnumber]) 

    print(" ")
    print("Using UrQMD Inputfile")
    print("Inputfile has %d baryons"%len(baryons))    
    print("Inputfile has %d mesons"%len(mesons))
    print("---------------------------")

    
    potential_number = 3 * len(baryons) + 2 * len(mesons)
    if (potential_number > maxnum):
        print("Hadrons generated from UrQMD data file: %d"%potential_number)
        print("Allowed maximum of processed partons: %d"%maxnum)
        print("Reducing parton number ....")

    while (potential_number > maxnum):
        baryons = sample(baryons, 9 * len(baryons)/10)
        mesons = sample(mesons, 9*len(mesons)/10)
        potential_number = 3 * len(baryons) + 2 * len(mesons)

    
    partons = make_partonlist(baryons, mesons)

    print("Processing: %d baryons"%len(baryons))
    print("Processing: %d mesons"%len(mesons))
    print("Processing total: %d partons"%len(partons))

    print(" ")

    pos = np.ndarray((len(partons), 4), dtype=np.float32)
    mom = np.ndarray((len(partons), 4), dtype=np.float32)
    col = np.ndarray((len(partons), 4), dtype=np.float32)
    for i in range(len(pos)):
        pos[i] = partons[i][0]
        mom[i][0:3] = partons[i][1][1:4] # (E, px, py, pz) |-> (px, py, pz, 0)
        col[i] = partons[i][2]
        

    vel = np.array(map(lambda x, y: x / y, mom,  np.sqrt(parton_mass**2 + np.array(map(lambda(x): x.dot(x), mom)))), dtype=np.float32)
    vel[:,3] =  parton_mass

    return pos, col, vel

def fountain_np(num):
    """ initialize 400 baryons and 500 mesons with 300 MeV energy """

    def create_anti(particle):
        return [1 - x for x in particle[:-1]] + [1]


    # Initial position of quarks
    pos = np.ndarray((num, 4), dtype=np.float32)
    
    r = np.random.rand(num) * .5
    phi = np.random.rand(num)*2*np.pi
    theta = np.random.rand(num)*2*np.pi

    pos[:,0] = r * np.sin(phi) * np.cos(theta)
    pos[:,1] = r * np.sin(phi) * np.sin(theta)
    pos[:,2] = r * np.cos(phi)
    pos[:,3] = 1.


    # Momentum and velocity
    mass = .01                 # take 0.01 GeV/c**2 as rest mass of quark
    mom = np.ndarray((num, 4), dtype=np.float32)
    r =  0.300                 # take 300 MeV/c as momentum
    phi = np.random.rand(num) * 2 * np.pi
    theta = np.random.rand(num) * 2 * np.pi
    
    mom[:,0] = r * np.sin(phi) * np.cos(theta)
    mom[:,1] = r * np.sin(phi) * np.sin(theta)
    mom[:,2] = r * np.cos(phi)
    mom[:,3] = 0.

    vel = np.array(map(lambda x, y: x / y, mom,  np.sqrt(mass**2 + np.array(map(lambda(x): x.dot(x), mom)))), dtype=np.float32)
    vel[:,3] =  mass 

    # Color charge
    X = [x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0] ] ] # rgb
    

    baryons = np.array([random.choice(X) for i in range( (num*6) / 11)], dtype=np.float32)

    number_pions = len(vel) - len(baryons)
    pions = np.array([random.choice(X) for i in range(number_pions / 2)], dtype=np.float32)
    anti_pions = np.array([create_anti(p) for p in pions],dtype=np.float32)
    col = np.concatenate((baryons,  anti_pions, pions))
    
    return pos, col, vel
    

def fountain(num, method = fountain_urqmd):
    """Initialize position, color and velocity arrays we also make Vertex
    Buffer Objects for the position and color arrays"""

    pos, col, vel = method(num)
    
    #create the Vertex Buffer Objects
    from OpenGL.arrays import vbo 
    pos_vbo = vbo.VBO(data=pos, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)
    pos_vbo.bind()
    col_vbo = vbo.VBO(data=col, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)
    col_vbo.bind()

    return (pos_vbo, col_vbo, vel)



from OpenGL.GL import *
import random
import numpy



def create_anti(particle):
    return [1 - x for x in particle[:-1]] + [1]

def fountain_np(num):
    """ initialize 400 baryons and 500 mesons with 300 MeV energy """

    # Initial position of quarks
    pos = numpy.ndarray((num, 4), dtype=numpy.float32)
    
    r = numpy.random.rand(num) * .5
    phi = numpy.random.rand(num)*2*numpy.pi
    theta = numpy.random.rand(num)*2*numpy.pi

    pos[:,0] = r * numpy.sin(phi) * numpy.cos(theta)
    pos[:,1] = r * numpy.sin(phi) * numpy.sin(theta)
    pos[:,2] = r * numpy.cos(phi)
    pos[:,3] = 1.


    # Momentum and velocity
    mass = 10.0                 # take 10 MeV/c**2 as rest mass of quark
    mom = numpy.ndarray((num, 4), dtype=numpy.float32)
    r =  300                    # take 300 MeV/c as momentum
    phi = numpy.random.rand(num) * 2 * numpy.pi
    theta = numpy.random.rand(num) * 2 * numpy.pi
    
    mom[:,0] = r * numpy.sin(phi) * numpy.cos(theta)
    mom[:,1] = r * numpy.sin(phi) * numpy.sin(theta)
    mom[:,2] = r * numpy.cos(phi)
    mom[:,3] = 0.

    vel = numpy.array(map(lambda x, y: x / y, mom,  numpy.sqrt(mass**2 + numpy.array(map(lambda(x): x.dot(x), mom)))), dtype=numpy.float32)
    vel[:,3] =  mass 

    # Color charge
    X = [x + [1] for x in [ [0,0,1], [0,1,0], [1,0,0] ] ] # rgb
    

    baryons = numpy.array([random.choice(X) for i in range( (num*6) / 11)], dtype=numpy.float32)

    number_pions = len(vel) - len(baryons)
    pions = numpy.array([random.choice(X) for i in range(number_pions / 2)], dtype=numpy.float32)
    anti_pions = numpy.array([create_anti(p) for p in pions],dtype=numpy.float32)
    col = numpy.concatenate((baryons,  anti_pions, pions))
    
    return pos, col, vel
    

def fountain(num):
    """Initialize position, color and velocity arrays we also make Vertex
    Buffer Objects for the position and color arrays"""

    pos, col, vel = fountain_np(num)
    
    #create the Vertex Buffer Objects
    from OpenGL.arrays import vbo 
    pos_vbo = vbo.VBO(data=pos, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)
    pos_vbo.bind()
    col_vbo = vbo.VBO(data=col, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)
    col_vbo.bind()

    return (pos_vbo, col_vbo, vel)



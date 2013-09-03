import physics
import initialize
import time
from numpy import savetxt, concatenate

#max number of particles
maxnum = 6000 
#time step for integration
dt = 1e-4
#number of timesteps
timesteps = 1000

class Simulation():
    def __init__(self, maxnum=maxnum, dt=dt):
        #set up initial conditions
        (pos, col, vel) = initialize.fountain(maxnum)
        num = len(vel)
        #create our OpenCL instance
        self.cle = physics.Particles(num, dt)
        self.cle.pushData(pos, col, vel)


    def run(self, timesteps=timesteps):
        #update or particle positions by calling the OpenCL kernel
        self.cle.execute(timesteps) 
    
    def save(self, fname="output.csv"):
        (pos, col, vel) = self.cle.pullData()
        liste = []
        for i in range(len(pos)):
            current = concatenate((pos[i][0:3], vel[i], col[i][0:3]))
            liste.append(current)
        savetxt(fname, liste, delimiter=",")


if __name__ == "__main__":
    MyRun = Simulation()
    then = time.clock()
    MyRun.run()
    delta = time.clock() - then
    print("Computation took %f msec"%(delta*1e3))
    MyRun.save()



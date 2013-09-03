import physics
import initialize

#max number of particles
maxnum = 600 
#time step for integration
dt = 1e-2

class Simulation():
    def __init__(self):
        #set up initial conditions
        (pos, col, vel) = initialize.fountain(maxnum)
        num = len(vel)
        #create our OpenCL instance
        self.cle = physics.Particles(num, dt)
        self.cle.loadData(pos, col, vel)


    def run(self):
        #update or particle positions by calling the OpenCL kernel
        self.cle.execute(20) 


if __name__ == "__main__":
    MyRun = Simulation()
    MyRun.run()




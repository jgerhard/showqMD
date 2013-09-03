import physics
import initialize
import time
from numpy import savetxt, concatenate

#max number of particles
maxnum = 6000 
#time step for integration
dt = 1e-2
#number of timesteps
run_time = 5                  # run time in fm/c
save_time = 0.1               # timesteps to be saved in fm/c

class Simulation():
    def __init__(self, maxnum=maxnum, dt=dt):
        #set up initial conditions
        (pos, col, vel) = initialize.fountain(maxnum)
        num = len(vel)
        #create our OpenCL instance
        self.cle = physics.Particles(num, dt)
        self.cle.pushData(pos, col, vel)
        self.totaltime = self.cle.totaltime


    def run(self, run_time = run_time):
        #update or particle positions by calling the OpenCL kernel
        self.cle.execute(run_time) 
        self.totaltime = self.cle.totaltime
    
    def save(self, fname="output.csv"):
        (pos, col, vel) = self.cle.pullData()
        liste = []
        for i in range(len(pos)):
            current = concatenate(([self.cle.totaltime],pos[i][0:3], vel[i], col[i][0:3]))
            liste.append(current)
    
        try:
            f_handle = file(fname, 'a')
            savetxt(f_handle, liste, delimiter=",")
        except:
            savetxt(fname, liste, delimiter=",")


if __name__ == "__main__":
    MyRun = Simulation()

    while (MyRun.totaltime < run_time):
        MyRun.save()
        MyRun.run(MyRun.totaltime + save_time)
    MyRun.save()




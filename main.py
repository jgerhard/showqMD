import physics
import initialize
import time
from numpy import savetxt, concatenate, sqrt, append
from find_neighbours import create_candidates
    
#max number of particles
maxnum = 6000
#time step for integration
dt = 1e-4                      # 1e-4 leaves mass of meson at hadronization constant for all timesteps
#number of timesteps
run_time = 1                  # run time in fm/c
save_time = 0.1               # timesteps to be saved in fm/c

class Simulation():
    def __init__(self, maxnum=maxnum, dt=dt):
        #set up initial conditions
        (pos, mommass, col) = initialize.fountain(maxnum)
        num = len(pos)
        #create our OpenCL instance
        self.cle = physics.Particles(num, dt)
        self.cle.pushData(pos, mommass, col)
        self.totaltime = self.cle.totaltime


    def run(self, run_time = run_time):
        #update or particle positions by calling the OpenCL kernel
        self.cle.execute(run_time) 
        self.totaltime = self.cle.totaltime
    
    def save(self, fname="output.csv", one_file=True, step_number=0):
        """ Outputs data into fname as csv file. If not one_file
        different files with step_number in their name are created
        Data is of format (t, x,y,z, px,py,pz, m, c1,c2,c3, fx,fy,fz"""
        (pos, mom, col, force) = self.cle.pullData()
        liste = []
        for i in range(len(pos)):
            current = concatenate(([self.totaltime],pos[i][0:3], mom[i], col[i][0:3])) #, force[i][0:3]))
            liste.append(current)
        if one_file:
            f_handle = file(fname, 'a')
            savetxt(f_handle, liste, delimiter=",")
            f_handle.close()
        else:
            savetxt(fname%step_number, liste, delimiter=",")
        
    def hadronize(self):
        """ Experimental isochronal hadronization on host """
        
        (pos, mommass, col, force) = self.cle.pullData()
        partons = list(concatenate( (pos, mommass, force, col) ,1).tolist())
        baryons, mesons = create_candidates(partons)
        hadrons = baryons + mesons
        liste = []
        for hadron in hadrons:
            current = concatenate(([self.totaltime],hadron))
            liste.append(current)
        f_handle = file("hadrons.csv", 'a')
        savetxt(f_handle, liste, delimiter=",")
        f_handle.close()



if __name__ == "__main__":
    MyRun = Simulation()
    MyRun.save()
    MyRun.hadronize()
    MyRun.run(run_time)
    MyRun.save()
    MyRun.hadronize()
    
    




import physics
import initialize
import time
from numpy import savetxt, concatenate, sqrt
from create_hadrons import create_meson, create_baryon

#max number of particles
maxnum = 6000 
#time step for integration
dt = 1e-6                      # 1e-4 leaves mass of meson at hadronization constant for all timesteps
#number of timesteps
run_time = 1                  # run time in fm/c
save_time = 0.001               # timesteps to be saved in fm/c

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
            current = concatenate(([self.totaltime],pos[i][0:3], mom[i], col[i][0:3], force[i][0:3]))
            liste.append(current)
    
        if one_file:
            try:
                f_handle = file(fname, 'a')
                savetxt(f_handle, liste, delimiter=",")
            except:
                savetxt(fname, liste, delimiter=",")
        else:
            savetxt(fname%step_number, liste, delimiter=",")
    
    def hadronize(self):
        """ Experimental isochronal hadronization on host """
        
        (pos, mommass, col, force) = self.cle.pullData()
        partons = concatenate((pos, mommass, force, col),1)
        E, E_pot, hadron = create_meson(*partons)
        return E, E_pot, hadron


def run():
    print("Simulating %f fm/c"%run_time)
    MyRun = Simulation()
    Es = []
    E_pots = []
    while (MyRun.totaltime < run_time):
        MyRun.run(MyRun.totaltime + save_time)
        E, E_pot, baryon = MyRun.hadronize()
        Es.append(E)
        E_pots.append(E_pot)
        MyRun.save()
    return Es, E_pots


if __name__ == "__main__":
    Es, E_pots = run()
    import pylab as pl
    pl.plot(Es, label="Kinetic + PartonMass")
    pl.plot(E_pots, label="Potential")
    pl.plot(pl.array(Es) + pl.array(E_pots), '+')
    pl.legend()
    pl.show()
    




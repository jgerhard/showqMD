import physics
import initialize
import time
from numpy import savetxt, concatenate, sqrt
from create_hadrons import create_meson

#max number of particles
maxnum = 6000 
#time step for integration
dt = 1e-4
#number of timesteps
run_time = 1                  # run time in fm/c
save_time = 0.02               # timesteps to be saved in fm/c

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
        distance, E, E_pot, meson = create_meson(*partons)
        return distance, E, E_pot, meson


def run():
    mesons = []
    distances = []
    Es = []
    E_pots = []
    MyRun = Simulation()
    distance, E, E_pot, meson = MyRun.hadronize()
    mesons.append(meson)
    distances.append(distance)
    Es.append(E)
    E_pots.append(E_pot)
    print("Simulating %f fm/c"%run_time)
    
    while (MyRun.totaltime < run_time):
        MyRun.run(MyRun.totaltime + save_time)
        distance, E, E_pot, meson = MyRun.hadronize()
        mesons.append(meson)
        distances.append(distance)
        Es.append(E)
        E_pots.append(E_pot)
        MyRun.save()
    return distances, Es, E_pots, mesons


if __name__ == "__main__":
    dists, Es, E_pots, mesons = run()
    import pylab as pl
    pl.plot(dists, label="Distance")
    pl.plot(E_pots, label="Potential Energy")
    pl.plot(Es, label="Kinetic + Mass")
    Masses = [pl.sqrt(meson[3]**2 - meson[4]**2 - meson[5]**2 - meson[6]**2) for meson in mesons]
    pl.plot(Masses, '+')
    pl.legend(loc=0)
    pl.show()

    




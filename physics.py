import pyopencl as cl
import sys
import numpy as np


class Particles(object):
    def __init__(self, num, dt, *args, **kwargs):
        self.clinit()
        self.loadProgram("cornell.cl");
        self.totaltime = 0.0
        self.num = num
        self.num_cl = np.uint32(num)
        self.dt = np.float32(dt)
        self.force = np.ndarray((num, 4), dtype=np.float32) 


    def pushData(self, pos, col, vel):
        """ Pushes particle data from host to device """
        mf = cl.mem_flags

        self.pos_A = pos
        self.vel_A = vel
        self.col = col

        self.pos_B = pos
        self.vel_B = vel


        #pure OpenCL arrays
        self.vel_A_cl = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.vel_A)
        self.pos_A_cl = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pos_A)
        
        self.col_cl = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.col)
        
        self.vel_B_cl = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.vel_A)
        self.pos_B_cl = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pos_A)

        self.cum_force = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.force.nbytes)
        
        self.queue.finish()

    def pullData(self):
        """ Pulls back device data to host """
        cl.enqueue_copy(self.queue, self.pos_A, self.pos_A_cl)
        cl.enqueue_copy(self.queue, self.vel_A, self.vel_A_cl)
        cl.enqueue_copy(self.queue, self.col, self.col_cl)
        cl.enqueue_copy(self.queue, self.cum_force, self.force)
        return (self.pos_A, self.col, self.vel_A, self.force)
        
        

    def execute(self, run_time):

        global_size = (self.num,)
        for i in range(1,64):   # choose group size
            if (self.num % i == 0) :
                local_size_threads = i


        local_size = (local_size_threads,)

        kernelargs = (self.pos_A_cl, 
                      self.vel_A_cl,
                      self.pos_B_cl,
                      self.vel_B_cl,
                      self.col_cl,
                      self.cum_force,
                      self.dt,
                      self.num_cl)

        kernelargsT = (self.pos_B_cl, 
                       self.vel_B_cl,
                       self.pos_A_cl,
                       self.vel_A_cl,
                       self.col_cl, 
                       self.cum_force,
                       self.dt,
                       self.num_cl)
        
        while (self.totaltime < run_time):
            self.program.nbody(self.queue, global_size, local_size, *(kernelargs))
            self.program.nbody(self.queue, global_size, local_size, *(kernelargsT)) # change role of kernelargs to do double buffered calc
            self.queue.finish() # spare time by putting this after the loop -- output of simulation time depends on this
            self.totaltime += 2*self.dt
            sys.stdout.write("\r t = {0} fm/c>".format(self.totaltime))
            sys.stdout.flush() 


 
    def clinit(self):
#        plats = cl.get_platforms()
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

    def loadProgram(self, filename):
        f = open(filename, 'r')
        fstr = "".join(f.readlines())
        self.program = cl.Program(self.ctx, fstr).build()



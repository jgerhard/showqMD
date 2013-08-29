import numpy as np
from pprint import pprint

def readfile(filename, header=True):
    """ Reads all lines in UrQMD f14-file and generates
    a List of data per Event. Includes header, if header=True """
    file = open(filename, 'r')
    flist = []
    for line in file: flist.append(line)
    events=[]
    single_event=[]
    while flist:
        stri = flist.pop()
        single_event.append(stri)
        if (stri.split()[0]=="UQMD"):
            single_event.reverse()
            if header:
                events.append(single_event)
            else:
                events.append(single_event[16:])
            single_event=[]
    events.reverse()
    return events



if __name__ == "__main__":
#    hadrons = range(1,27) + [101, 104] # want N*, \Delta, \pi and \rho
    hadrons = [1,101, 104]           # consider only n, \pi^0, and \rho^0  (no charge in qMD)
    liste = readfile("test.f14",header=False)
    pprint(filter(lambda(line): int(line.split()[9])in hadrons, liste[0]))

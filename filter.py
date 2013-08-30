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



def separate_hadrons(event, meson_itypes = [101, 104], baryon_itypes = [1]):
    """ Takes one event and seperates baryons from mesons 
    output format are two lists of baryons and mesons, with each particle
    being a sub-list of form (r_x, r_y, r_z, E, p_x, p_y, p_z, m)
    """

    mesons_strings = filter(lambda(line): int(line.split()[9])in meson_itypes, event) # and int(line.split()[11])==0
    baryons_strings = filter(lambda(line): int(line.split()[9])in baryon_itypes, event) # and int(line.split()[11])==0
    
    mesons = [map(lambda(x): float(x), line.split()[1:9]) for line in mesons_strings]
    baryons = [map(lambda(x): float(x), line.split()[1:9]) for line in baryons_strings]

    return baryons, mesons
    

    
if __name__ == "__main__":
    liste = readfile("test.f14",header=False)
    baryons, mesons = separate_hadrons(liste[0]) # take only first event
    pprint(len(baryons))
    pprint(len(mesons))
    pprint(len(liste[0]))
        

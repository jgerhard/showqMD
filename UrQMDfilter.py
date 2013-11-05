from create_partons import create_triplet, create_duplet
import numpy as np

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



def separate_hadrons(event, meson_itypes =range(-140,-100)+range(100,141), baryon_itypes = range(1,56)):
    """ Takes one event and seperates baryons from mesons 
    output format are two lists of baryons and mesons, with each particle
    being a sub-list of form (r_x, r_y, r_z, E, p_x, p_y, p_z, m, ityp, 2*I_z, charge)
    """

    mesons_strings = filter(lambda(line): int(line.split()[9])in meson_itypes, event) # and int(line.split()[11])==0
    baryons_strings = filter(lambda(line): int(line.split()[9])in baryon_itypes, event) # and int(line.split()[11])==0
    
    mesons = [map(lambda(x): float(x), line.split()[1:12]) for line in mesons_strings]
    baryons = [map(lambda(x): float(x), line.split()[1:12]) for line in baryons_strings]
    mesons = filter(lambda(x): x[-2] != 0, mesons) # UrQMD denotes photons as mass=0 mesons
    return baryons[:2], mesons[:1]
    
def make_partonlist(baryons, mesons):
    """ Takes list of baryons and mesons and creates list of partons """
    partons = []
    for baryon in baryons:
        partons.extend(create_triplet(baryon))
    for meson in mesons:
        partons.extend(create_duplet(meson))
    return partons


    
if __name__ == "__main__":
    np.random.seed(42)
    events = readfile("test.f14",header=False)
    baryons, mesons = separate_hadrons(events[0]) # take only first event
    liste = make_partonlist(baryons, mesons)
    np.savetxt("baryons.csv", baryons, delimiter=",")
    np.savetxt("mesons.csv", mesons, delimiter=",")

import numpy as np
import pprint

def readfile(filename, header=True):
    """ Reads all lines in UrQMD f14-file and generates
    a List of data per Event. Includes header, if header=True """
    file = open(filename, 'r')
    flist = []
    for line in file: flist.append(line)
    events=[]
    single_event=[]
    numofevents = 0
    while flist:
        stri = flist.pop()
        single_event.append(stri)
        if (stri.split()[0]=="UQMD"):
            single_event.reverse()
            events.append(single_event)
            single_event=[]
            numofevents += 1
    events.reverse()
    return numofevents, events



if __name__ = "__main__":
    num, liste = readfile("test.f14")
    print("Number of events %d", num)
    pprint(liste)

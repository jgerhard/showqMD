from numpy import sqrt
from create_hadrons import *
from itertools import combinations
def create_measure(parton):
    """ creates a distance to given parton function """
    def delta(other):
        return sqrt( (parton[0]-other[0])**2 + 
                     (parton[1]-other[1])**2 +
                     (parton[2]-other[2])**2
                 )
    
    def cmp(X, Y):
        return(1 if (delta(X) < delta(Y)) else -1)
    return delta, cmp



def is_white(*partons):

    total_colour = [0,0,0]
    for parton in partons:
        total_colour[0] += parton[12]
        total_colour[1] += parton[13]
        total_colour[2] += parton[14]
    return (total_colour == [1,1,1]) or (total_colour == [2,2,2]) # 222 is the anti_baryon case
    

def create_candidates(all_partons, max_dist = 1.0, MAXITER = None):
    """ Takes list of all partons and 
    subdivides to meson and baryon candidates.
    Firstly excludes all partons with more 
    than 1.0 fm/c distance"""
    
    print
    print("--------------------")
    print("Hadronization...")
    print("Method: nearest neighbour")
    print("--------------------")
    print("    1 fm/c distance max")
    mesons = []
    baryons = []
    i = 0
    if (not MAXITER):
        MAXITER = 10 * len(all_partons)

    while ( (all_partons) and (i < MAXITER) ):
        i += 1
        X = all_partons.pop(0)
        dist, compare = create_measure(X)
        candidates = filter( lambda z: (dist(z) <= max_dist), all_partons)
        candidates.sort(cmp=compare)

        if candidates and (is_white(X, candidates[0])): # meson case
            mesons.append(create_meson(X, candidates[0]))
            all_partons.remove(candidates[0])
        elif len(candidates)>1 and (is_white(X, candidates[0], candidates[1])): # baryon case
            baryons.append(create_baryon(X, candidates[0], candidates[1]))
            all_partons.remove(candidates[0])
            all_partons.remove(candidates[1])
        else:
            all_partons.append(X)   # did not find correct neighbours yet
   
    if all_partons:
        print("    Partons left over: %d" %len(all_partons))
        print("------------------------")

    i = 0
    max_dist = 5.0
    MAXITER = 10 * len(all_partons)
    if all_partons:
        print("     5 fm/c distance max")
    while ( (all_partons) and (i < MAXITER) ):
        i += 1
        X = all_partons.pop(0)
        dist, compare = create_measure(X)
        candidates = filter( lambda z: (dist(z) <= max_dist), all_partons)
        candidates.sort(cmp=compare)

        if candidates and (is_white(X, candidates[0])): # meson case
            mesons.append(create_meson(X, candidates[0]))
            all_partons.remove(candidates[0])
        elif len(candidates)>1 and (is_white(X, candidates[0], candidates[1])): # baryon case
            baryons.append(create_baryon(X, candidates[0], candidates[1]))
            all_partons.remove(candidates[0])
            all_partons.remove(candidates[1])
        else:
            all_partons.append(X)   # did not find correct neighbours yet
    
    if all_partons:
        print("    Partons left over: %d" %len(all_partons))
        print("------------------------")

    if all_partons:
        print("    ignore max distance")

    i = 0
    MAXITER = 10 * len(all_partons)

    while ( (all_partons) and (i < MAXITER) ):
        i += 1
        X = all_partons.pop(0)
        dist, compare = create_measure(X)
        candidates = filter(lambda(x): True, all_partons)
        candidates.sort(cmp=compare)
        if candidates and (is_white(X, candidates[0])): # meson case
            mesons.append(create_meson(X, candidates[0]))
            all_partons.remove(candidates[0])
        elif len(candidates)>1 and (is_white(X, candidates[0], candidates[1])): # baryon case
            baryons.append(create_baryon(X, candidates[0], candidates[1]))
            all_partons.remove(candidates[0])
            all_partons.remove(candidates[1])
        else:
            all_partons.append(X)   # did not find correct neighbours yet

    if all_partons:
        print("    Partons left over: %d" %len(all_partons))
        print("------------------------")
    
    if all_partons:
        print("Method: Combinatorical search")

    while(all_partons):
        meson_candidates = combinations(all_partons, 2)
        for a, b in meson_candidates:
            if is_white(a,b):
                mesons.append(create_meson(a,b))
                all_partons.remove(a)
                all_partons.remove(b)
                break
        baryon_candidates = combinations(all_partons, 3)
        for a,b,c in baryon_candidates:
            if is_white(a,b,c):
                baryons.append(create_baryon(a,b,c))
                all_partons.remove(a)
                all_partons.remove(b)
                all_partons.remove(c)
                break
    if all_partons:
        print("    Partons left over: %d" %len(all_partons))
        print("------------------------")
    else:
        print("------------------------")
        print("|  ALL PARTONS GROUPED |")
        print("------------------------")
    return baryons, mesons
    



        
    
    

from numpy import sqrt

def create_measure(parton):
    """ creates a distance to given parton function """
    def delta(other):
        return sqrt( (parton[0]-other[0])**2 + 
                     (parton[1]-other[1])**2 +
                     (parton[2]-other[2])**2
                 )



def create_candidates(all_partons, max_dist = 1.0):
    """ Takes list of all partons and 
    subdivides to meson and baryon candidates.
    Firstly excludes all partons with more 
    than 1.0 fm/c distance"""

    X = all_partons.pop()
    dist = create_measure(X)
    candidates = filter(dist(z) <= max_dist, all_partons)
    
    print candidates
    
    

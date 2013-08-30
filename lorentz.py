from numpy import array, dot, sqrt

def lorentz(beta, fourVector):
    """ Takes as input relative beta from O -> O'
    and calculates fourVector -> fourVector' """

    beta2 = dot(beta, beta)
    gamma = 1./sqrt(1. - beta2)
    
    Lambda = array([ [gamma, -gamma*beta[0], -gamma*beta[1], -gamma*beta[2]],
                     [-gamma * beta[0], 1 + (gamma-1)*beta[0]**2/beta2, (gamma-1)*beta[0]*beta[1]/beta2, (gamma-1)*beta[0]*beta[2]/beta2],
                     [-gamma * beta[1], (gamma-1)*beta[1]*beta[0]/beta2, 1 + (gamma-1)*beta[1]**2/beta2, (gamma-1)*beta[1]*beta[2]/beta2],
                     [-gamma * beta[2], (gamma-1)*beta[2]*beta[0]/beta2, (gamma-1)*beta[2]*beta[1]/beta2, 1 + (gamma-1)*beta[2]**2/beta2]])

    return dot(Lambda, fourVector)

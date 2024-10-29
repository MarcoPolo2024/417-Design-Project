import numpy as np

def HouseholdP2(x_intial:float,scheme_function:function,scheme_prime:function,scheme_double_prime:function)->float:
  max_iterations = 1000
  while abs(scheme_function(x_intial)) > 1e-8:
    x_intial = x_intial - ((2*scheme_function(x_intial))/(scheme_prime(x_intial) - (scheme_prime(x_intial)**2-scheme_function(x_intial)*scheme_double_prime(x_intial))**0.5))
    max_iterations -=1
    if max_iterations ==0:
      print('The scheme didn\'t converge')
      break
  return x_intial

def Householder(x_position:float,section_supersonic:bool,area_function:function)->float:
  gamma = 1.148
  P = 2/(gamma+1)
  Q = 1-P
  if section_supersonic==False:
    R = (area_function(x_position))**2
    a = P**(1/Q)
    r = (R-1)/(2*a)
    x_intial = 1 / ((1+r)+np.sqrt(r*(r+2)))
    f = lambda X : (P+Q*X)**(1/Q) - R*X
    f_prime = lambda X: (P+Q*X)**((1/Q)-1) - R
    f_double_prime = lambda X: P*(P+Q*X)**((1/Q)-2)
    x_final = HouseholdP2(x_intial,f,f_prime,f_double_prime)
    return (x_final)**0.5
  if section_supersonic == True:
    R = (area_function(x_position))**(2*Q/P)
    a = Q**(1/P)
    r = (R-1)/(2*a)
    x_intial = 1 / ((1+r)+np.sqrt(r*(r+2)))
    f = lambda X : (P*X+Q)**(1/P) - R*X
    f_prime = lambda X: (P*X+Q)**((1/P)-1) - R
    f_double_prime = lambda X: Q*(P*X+Q)**((1/P)-2)
    x_final = abs(HouseholdP2(x_intial,f,f_prime,f_double_prime))
    return 1/(x_final)**0.5

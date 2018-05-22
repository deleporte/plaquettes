import random
import scipy.linalg
import numpy as np

def binary(n,d):
    if n>= 2**d:
        print "overflow while converting integer"
        digits=[0]*d
    else:
        digits=[]
        for k in range(d):
            digits.append(n%2)
            n=n/2
        return digits

def decimal(digits):
    cdigits=digits[:]
    n=0
    while len(cdigits)>0:
       n*=2
       n+= cdigits.pop()
    return n

def generate_C(d,m=1,M=2):
    """
    Generates a bunch of randomly chosen coefficients, far from zero
    """
    C=[]
    for k in range(2**d):
        C.append(random.random()*(M-m)+m)
    return C

def generate_biaised_C(d,m=1,M=2):
    """
    Generates a bunch of randomly chosen coefficients, far from zero
    """
    C=[]
    for k in range(2**d):
        C.append((random.random()*(M-m)+m)*(0.1+k%2))
    return C

def generate_deterministic_C(d,m=1,M=2):
    C=[]
    for k in range(2**d):
        C.append(m+(k%2)*(M-m))
    return C
 
def chain(C=generate_C(4)):
    #creates the transition matrix of a Markov chain with random plaquette weights, then prints its spectral gap
    #the transition matrix
    M=np.zeros((2**(d-1),2**(d-1)))
    for i in range(2**(d-1)):
        M[i,2*(i%(2**(d-2)))]=C[2*i]/(C[2*i]+C[2*i+1])
        M[i,2*(i%(2**(d-2)))+1]=C[2*i+1]/(C[2*i]+C[2*i+1])
    Sp,vects=scipy.linalg.eig(M,left=True,right=False)
    mu=vects.T[0]
    print Sp
    maxval=0.
    for val in Sp:
        if abs(val)>maxval and abs(1-val)>1.0e-14:
            maxval=abs(val)
    print np.log(maxval)
    return Sp,mu
    

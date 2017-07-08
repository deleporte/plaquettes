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
 
def chain(d=4):
    #creates the transition matrix of a Markov chain with random plaquette weights, then prints its spectral gap
    #random weights away from zero
    C=[]
    for k in range(2**d-1):
        C.append(random.random()+0.1)
    C.append(0.)
    #the transition matrix
    M=np.zeros((2**(d-1),2**(d-1)))
    for i in range(2**(d-1)):
        M[i,2*(i%(2**(d-2)))]=C[2*i]/(C[2*i]+C[2*i+1])
        M[i,2*(i%(2**(d-2)))+1]=C[2*i+1]/(C[2*i]+C[2*i+1])
    Sp=scipy.linalg.eigvals(M)
    print Sp
    maxval=0.
    for val in Sp:
        if abs(val)>maxval and abs(1-val)>1.0e-14:
            maxval=abs(val)
    print np.log(maxval)
    return Sp
    

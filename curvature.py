import numpy as np

def generate_C(d,m=1,M=2):
    C=[]
    for k in range(2**d):
        C.append(random.random()*(M-m)+m)
    return C

def sparse_mult(mat,C):
    """returns an efficient multiplication of mat by the transfer matrix.
    Runs in quadratic time.
    """
    prod=np.zeros((len(C)/2,len(C)/2),dtype=float)
    for i in range(len(C)/2):
        for j in range(len(C)/2):
            prod[i,j]=mat[i,(2*j)%(len(C)/2)]*C[2*j]
            prod[i,j]+=mat[i,(2*j+1)%(len(C)/2)]*C[2*j+1]
            prod[i,j]/=C[2*j]+C[2*j+1]
    return prod
    

def curvature(C,ERROR=1.0e-14):
    """Computes, via a Markov chain, the curvature tensor
    for constant plaquettes up to some prescribed error."""
    #compute the sequence (M**i) until convergence
    powerseq=[np.eye(len(C)/2)]
    converged=False
    while not converged:
        current=powerseq[-1]
        new=sparse_mult(current,C)
        powerseq.append(new)
        dist=0.
        for k in range(len(C)/2):
            dist+=abs(new[k,0]-new[k,1])
        if dist<ERROR:
            converged=True
        print dist
    print "Markov chain convergence after",len(powerseq),"steps"
    #the lines of the last matrix are the equilibium measure
    eq_m=new.T[0]
    print "Equilibrium measure: ", eq_m
    G=np.zeros((len(C)/2,len(C)/2),dtype=float)
    for i in range(len(C)/2):
        for j in range(len(C)/2):
            for k in range(len(powerseq)):
                current=powerseq[k]
                G[i,j]+=eq_m[j]*(current[i,j]-eq_m[i])
            for k in range(len(powerseq)-1):
                current=powerseq[k+1]
                G[i,j]+=eq_m[i]*(current[j,i]-eq_m[j])
    return G

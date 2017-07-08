def double_transfer_powers(C,kmax):
    """It should return a matrix of size len(C)*(len(C)**2).
    dT_k[a,b] is the probability, if we are in a at 0, to be in b at k-d."""
    Tfor=np.zeros((len(C),len(C)),dtype=float)
    Tbac=np.zeros((len(C),len(C)),dtype=float)
    for i in range(len(C)):
        j=((2*i)%len(C))
        Tfor[i,j]=C[j]/(C[j]+C[j+1])
        Tfor[i,j+1]=C[j+1]/(C[j]+C[j+1])
        j=i/2
        Tbac[i,j]=C[j]/(C[j]+C[j+len(C)/2])
        Tbac[i,j+len(C)/2]=C[j+len(C)/2]/(C[j]+C[j+len(C)/2])
    dT=[]
    dT0=np.zeros((len(C),len(C)**2),dtype=float)
    for i in range(len(C)**2):
        dT0[i%len(C),i]=1.
        for k in range(len(C)-1):
            dT0[i%len(C),i]*=Tbac[(i/2**k)%len(C),(i/2**(k+1))%len(C)]
    dT.append(dT0)
    for k in range(kmax):
        dT_next=np.zeros((len(C),len(C)**2),dtype=float)
        #We implement a fast algorithm for dT_next=np.dot(dT[k],Tfor)
        for i in range(len(C)):
            for j in range(len(C)**2):
                dT_next[i,j]=dT[k][i,j/2]+dT[k][i,j/2+len(C)**2/2]
                dT_next[i,j]*=Tfor[(j/2)%len(C),j%len(C)]
        dT.append(dT_next)
    return dT
        

def Ising_Transverse(J,C):
    """Returns the tangent vector closest to [H(J)|\phi(C)>], and the angle."""
    

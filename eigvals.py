import numpy as np
import scipy.linalg

def eigvals(A,wr,wi,vecl,vecr):
    dim=len(wr)
    pA=np.array(A)
    pA=np.reshape(pA,(dim,dim))
    #diag
    pw,pvecr=scipy.linalg.eig(pA,overwrite_a=True)
    pvecl=scipy.linalg.inv(pvecr)
    wr=pw.real.tolist()[:]
    wi=pw.imag.tolist()[:]
    #under default python behaviour eigenvalues appear in pairs
    for i in range(dim):
        if(abs(wi[i])<0.000001):
            wi[i]=0;
            vecl[i*dim:(i+1)*dim]=pvecl.T[i].real.tolist()[:]
            vecr[i*dim:(i+1)*dim]=pvecr.T[i].real.tolist()[:]
           # for j in range(dim): #a more pythonic way ?
                #vecl[i*dim+j]=real(pvecl[j,i])
                #vecr[i*dim+j]=real(pvecr[j,i])
        else:
            vecl[i*dim:(i+1)*dim]=pvecl.T[i].real.tolist()[:]
            vecr[i*dim:(i+1)*dim]=pvecr.T[i].real.tolist()[:]
            vecl[(i+1)*dim:(i+2)*dim]=pvecl.T[i].imag.tolist()[:]
            vecr[(i+1)*dim:(i+2)*dim]=pvecr.T[i].imag.tolist()[:]
            i+=1
    return (A,wr,wi,vecl,vecr)
    

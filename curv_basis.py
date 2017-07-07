# a tentative basis

def testbasis(d=2):
    if d<=2:
        basis=np.array([[1,1],[1,-1]])
    else:
        basis=np.zeros((2**(d-1),2**(d-1)),dtype=float)
    for k in range(2**(d-2)-1):
        basis[k,k]+=1.
        basis[k,2**(d-2)+k]+=1.
        basis[k,2*k]+=-1.
        basis[k,2*k+1]+=-1.
        basis[2**(d-2)+k,k]+=1.
    basis[2**(d-2)-1]=1.
    basis[-1,-1]=1.
    return basis

def blindsearch(G,d):
    zeros=[]
    for k in range(3**(2**(d-1))):
        v=[]
        kcopy=k
        for i in range(2**(d-1)):
            if kcopy%3==2:
                v.append(-1.)
            else:
                v.append(float(kcopy%3))
            kcopy=kcopy/3
        if sum(abs(np.dot(v,G)))<1.0e-10:
            print v
    return 0

def lshift(c,d):
    return (2*c)%(2**d)+c/(2**(d-1))

def rtldigits(c,d):
    """Returns the reversed digits of c in base 2, of length d"""
    ccopy=c
    dig=[]
    for j in range(d):
        dig.append(ccopy%2)
        ccopy=ccopy/2
    return np.array(dig)

def rotation(d=3):
    """Rotates from the base along z, to the base along x."""
    U=np.zeros((2**d,2**d))
    for c in range(2**d):
        for k in range(2**d):
            U[c,k]=(-1)**(sum((1-rtldigits(c,d))*(1-rtldigits(k,d))))/2**(0.5*d)
    return U

def fullbasis(d=3):
    U=np.zeros((2**d,2**d))
    line=0
    #First we look at the zero vectors
    #they look like 1c-c1 for some c
    #but we have to give an orthogonal basis !
    for j in range(d-2):
        for c in range(2**j):
            #we consider the vectors from 1...10c0 to 0c01.....1 of size d-1
            #here they form a Toeplitz matrix of rank d-2-j,
            #which can easily be diagonalized
            for k in range(d-2-j):
                for l in range(d-2-j):
                    vect=2**(d-1)-2**(l+j+2)+c*2**(l+1)+2**l-1
                    val=np.sin(np.pi*(k+1)*(l+1)/(d-1-j))
                    U[line+k,2**(d-1)+vect]+=val
                    U[line+k,2*vect+1]-=val
                U[line+k]/=np.sqrt(np.dot(U[line+k],U[line+k]))
            line+=d-2-j
    #then we only forgot the vectors 1.....10 to 01......1
    for k in range(d-1):
        for l in range(d-1):
            vect=2**(d-1)-1-2**l
            val=np.sin(np.pi*(k+1)*(l+1)/d)
            U[line+k,2**(d-1)+vect]+=val
            U[line+k,2*vect+1]-=val
        U[line+k]/=np.sqrt(np.dot(U[line+k],U[line+k]))
    #The vector 1...1 is itself orthogonal
    line+=d-1
    U[line,2**d-1]=1.
    line+=1
    #There are several types of orthogonal (nonzero) vectors:
    #All 0c0 vectors
    for c in range(2**(d-2)):
        U[line+c,2*c]=1
    line+=2**(d-2)
    #All 10c0+0c01 vectors
    for c in range(2**(d-3)):
        U[line+c,2*c+2**(d-1)]=1./np.sqrt(2)
        U[line+c,4*c+1]=1./np.sqrt(2)
    line += 2**(d-3)
    #Some weird vectors
    for c in range(2**(d-3)):
        current=2*c+2**(d-1)+2**(d-2)
        shiftlist=[current]
        while current>=2**(d-1):
            current=lshift(current,d)
            shiftlist.append(current)
        for shifted in shiftlist:
            U[line+c,shifted]=1./np.sqrt(len(shiftlist))
    return U

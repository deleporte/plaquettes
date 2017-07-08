def inverse_transition_matrix(C,ERROR=1.0e-14):
    """Computes in small time the pseudo inverse of the I-transition matrix.
    This matrix should be applied only to vectors with zero sum.

    This returns S such that S(I-T)=(I-T)S is the projector on the set
    of vectors of zero sum, along the equilibrium matrix.
    """    
    powerseq=Markov_powers(C,ERROR)
    inverse=np.zeros((len(C),len(C)),dtype=float)
    for power in powerseq:
        inverse+=power
    return inverse

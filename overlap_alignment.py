import affine

def overlap_align(x, y, submatrix, g):
    """Computes an optimal overlap pairwise alignment from sequence x to sequence y
    with an affine gap scoring function.
        
    Args:
        x: a string representing the first sequence
        y: a string representing the second sequence
        submatrix: a substitution matrix that also contains character-specific space scores
        g: the gap existence score for gaps
    Returns:
        A tuple, (score, alignment), where score is a numeric value giving the score of the
        alignment and alignment is a list of two strings
    """
    
    # First, initialize matrices
    M, Ix, Iy = affine.initialize_matrices(x, y, submatrix, g)
    
    # Next, fill out the rest of the matrices, including tracebacks
    M, Ix, Iy, M_T, Ix_T, Iy_T = affine.fill_matrices(M, Ix, Iy, x, y, submatrix, g)
    
    
    # BEGIN TRACEBACK
    
    # get starting cell of traceback
    start = affine.get_start(x, y, M, Ix, Iy)
    
    # Get the score
    score = start[2]
    
    # get the alignment using the traceback algorithm
    alignment = affine.traceback(x, y, start, M_T, Ix_T, Iy_T)
    
    
    return (score, alignment)
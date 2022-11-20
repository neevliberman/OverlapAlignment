NEGATIVE_INFINITY = float("-inf")

def traceback(x, y, start, M_T, Ix_T, Iy_T):
    """Use the traceback algorithm for affine gap to find the final alignment 
    using the traceback matrices initialized in initialize_tracebacks.
    
    Args:
        x: a string representing the first sequence
        y: a string representing the second sequence
        start: a tuple representing (matrix_name, index_of_max, max_value). i.e., ('M', 2, 20), for which to start the traceback from
        M_T: the traceback matrix for the M matrix
        Ix_T: the traceback matrix for the Ix matrix
        Iy_T: the traceback matrix for the Iy matrix
    Returns: a list of two sequences aligned using the traceback algorithm; e.g., ['TAT ', ' ATA']
    """
    
    
    # initalize the final alignment which will be returned
    alignment = ['', '']
    
    
    # keep track of current j to break out of loop when the first column is reached
    curr_j = start[1]
    
    # keep track of current i to access x's characters
    curr_i = len(x)
    
    curr_matrix = start[0]
    
    # initialize second alignment, based on y, with any characters that won't be aligned (i.e., past the current j)
    # then, add spaces to first alignment, based on x, according to the number of characters just added to the second alignment.
    try:
        alignment[1] = y[curr_j:]
        alignment[0] = ' '*len(alignment[1])
    except:
        pass
    
    while curr_j > 0:
        
        if curr_matrix == "M":
            
            # get next current matrix based on last iterations i and j values, to examine 
            # whatever matrix the curr cell in T matrix points to.
            curr_matrix = M_T[curr_i][curr_j]
            
            # decrement i and j, since curr matrix was M
            curr_j -= 1
            curr_i -= 1
            
            # add next character from x to first alignment sequence.
            alignment[0] = x[curr_i] + alignment[0]
            
            # add next character from y to second alignment sequence
            alignment[1] = y[curr_j] + alignment[1]


            
        elif curr_matrix == "Ix":
            
            # get next current matrix based on last iterations i and j values, to examine 
            # whatever matrix the curr cell in T matrix points to.
            curr_matrix = Ix_T[curr_i][curr_j]
            
            # decrement only i, since curr matrix was Ix
            curr_i -= 1
            
            # add next character from x to first alignment sequence.
            alignment[0] = x[curr_i] + alignment[0]
            
            # add a dash character to match first alignment sequence's character
            alignment[1] = '-' + alignment[1]
            
            
            

            
        elif curr_matrix == "Iy":
            
            # get next current matrix based on last iterations i and j values, to examine
            # whatever matrix the curr cell in T matrix points to.
            curr_matrix = Iy_T[curr_i][curr_j]
            
            # decrement only j, since curr matrix was Iy
            curr_j -= 1
            
            # add a dash character to match second alignment sequence's character
            alignment[0] = '-' + alignment[0]
            
            # add next character from y to second alignment sequence.
            alignment[1] = y[curr_j] + alignment[1]
            
            
            
               
    # if loop condition is unsatisfied (i.e. current column is the first column), a ppend rest of x to alignment[0].
    # Also, add spaces to y to match the # of characters in x.
    alignment[0] = x[:curr_i] + alignment[0]
    alignment[1] = ' '*len(x[:curr_i]) + alignment[1]
    
    return alignment
            

def get_start(x, y, M, Ix, Iy):
    """Gets the starting cell of the traceback.
    
    Args:
        x: a string representing the first sequence
        y: a string representing the second sequence
        M: the M matrix, where values are the best scores given that x[i] is aligned to y[j]
        Ix: the Ix matrix, where values are the best scores given that x[i] is aligned to a gap
        Iy: the Iy matrix, where values are the best scores given that y[j] is aligned to a gap
    Returns: a tuple representing (matrix_name, index_of_max, max_value). i.e., ('M', 2, 20)
    """
    
    # Get maximum of the last row of each matrix, prioritizing larger j indices. 
    # Each max is structured as a tuple (matrix_name, index_of_max, max_value)
    Ix_max = ("Ix", len(y) - Ix[-1][::-1].index(max(Ix[-1])), max(Ix[-1]))
    M_max = ("M", len(y) - M[-1][::-1].index(max(M[-1])), max(M[-1]))
    Iy_max = ("Iy", len(y) - Iy[-1][::-1].index(max(Iy[-1])), max(Iy[-1]))
    
    # Get maximum of these maximums, prioritizing Ix <- M <- Iy
    start = max([Ix_max, M_max, Iy_max], key=lambda x:x[2])
    
    return start

def initialize_tracebacks(x, y):
    """Initializes all three traceback matrices' yet-to-be calculated cells to None, except for Iy_T, 
    where each successive Iy value in the first row of the matrix refers to the previous Iy value.
    
    Args:
        x: a string representing the first sequence
        y: a string representing the second sequence
    Returns: a tuple of three initalized 2-dimensional traceback matrices (M_T, Ix_T, Iy_T).
    """
    
    M_T = [[None]*(len(y)+1) for z in range(len(x)+1)]
    
    Ix_T = [[None]*(len(y)+1) for z in range(len(x)+1)]
     
    # Iy's in first row just refer to Iy's to the left
    Iy_T = [[None]*(len(y)+1) for z in range(len(x)+1)]
    for j in range(1, len(Iy_T[0])):
        Iy_T[0][j] = "Iy"
        
    return M_T, Ix_T, Iy_T



def initialize_matrices(x, y, submatrix, g):
    """Initializes the first column and row of all three matrices using ruleset 
    as explained in the notebook. Sets all other yet-to-be calculated values to None.
    
    Args:
        x: a string representing the first sequence
        y: a string representing the second sequence
        submatrix: a substitution matrix that also contains character-specific space scores
        g: the gap existence score for gaps
    Returns: a tuple of three initalized 2-dimensional matrices (M, Ix, Iy).
    """
    
    M = [[0]+[NEGATIVE_INFINITY]*len(y), 
         *[[NEGATIVE_INFINITY] + [None]*len(y) for i in range(len(x))]]
    
    Ix = [[0]+[NEGATIVE_INFINITY]*len(y),
          *[[0]+[None]*len(y) for i in range(len(x))]]
    
    
    Iy = [[g]+[None]*len(y),
          *[[NEGATIVE_INFINITY] + [None]*len(y) for i in range(len(x))]]
    
    for j in range(1, len(y)+1):
        Iy[0][j] = Iy[0][j-1] + submatrix[(y[j-1], '-')]
    
    return M, Ix, Iy



def fill_matrices(M, Ix, Iy, x, y, submatrix, g):
    """Fills the rest of the rows and columns, after the the first row and column were filled.
    Fills according to ruleset described in the notebook, typically filling based off of 
    the previous cell in one of 3 directions.
    
    Args:
        M: the M matrix, where values are the best scores given that x[i] is aligned to y[j]
        Ix: the Ix matrix, where values are the best scores given that x[i] is aligned to a gap
        Iy: the Iy matrix, where values are the best scores given that y[j] is aligned to a gap
        x: a string representing the first sequence
        y: a string representing the second sequence
        submatrix: a substitution matrix that also contains character-specific space scores
        g: the gap existence score for gaps
    Returns: a tuple of six 2-dimensional matrices: (M, Ix, Iy), the filled matrices, 
    and (M_T, Ix_T, Iy_T), the traceback matrix for each respective matrix.
    """
    

    # Here, we initialize the traceback matrices we will use to keep track of which index to jump back to, 
    # using strings like "M", "Ix", and "Iy".
    # We use one traceback matrix for each matrix M, Ix, and Iy.
    
    M_T, Ix_T, Iy_T = initialize_tracebacks(x, y)
    
    # Note: when accessing submatrix, we must use previous indices instead of current 
    # (i.e. y[j-1] instead of j)
    # This is because while the matrix starts at index zero, the sequences are lined up along the gridlines
    # meaning that their index is alwyas 1 less than the current matrix cell.
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            
            ###### BEGIN M-STEP
            # list to store each potential max value in the [i-1][j-1] position for curr M cell
            M_vars = [Ix[i-1][j-1] + submatrix[x[i-1], y[j-1]],
                M[i-1][j-1] + submatrix[x[i-1], y[j-1]],
                Iy[i-1][j-1] + submatrix[x[i-1], y[j-1]]]
            
            # get max of M_vars
            M[i][j] = max(M_vars)
            
            # based on what this max is, record in traceback. 
            # Even though its impossible for the pos to be -inf, consistency is good.
            if M[i][j] != NEGATIVE_INFINITY:
                M_T[i][j] = ["Ix", "M", "Iy"][M_vars.index(M[i][j])]
            ###### END M-STEP
        
            
            ###### BEGIN IX-STEP
            # list to store each potential max value in the [i-1][j] position for curr Ix cell
            Ix_vars = [Ix[i-1][j] + submatrix[x[i-1], '-'],
                       M[i-1][j] + g + submatrix[x[i-1], '-']]
            
            # get max of Ix_vars
            Ix[i][j] = max(Ix_vars)
            
            # based on what this max is, record in traceback. Do not point back to negative infinitys!!!!
            if Ix[i][j] != NEGATIVE_INFINITY:
                Ix_T[i][j] = ["Ix", "M"][Ix_vars.index(Ix[i][j])]
            ###### END IX-STEP
            
            
            ###### BEGIN IY-STEP
            # list to store each potential max value in the [i][j-1] position for curr Iy cell
            Iy_vars = [M[i][j-1] + g + submatrix[y[j-1], '-'],
                       Iy[i][j-1] + submatrix[y[j-1], '-']]
            
            # get max of IY_vars
            Iy[i][j] = max(Iy_vars)
            
            # based on what this max is, record in traceback. Do not point back to negative infinitys!!!!
            if Iy[i][j] != NEGATIVE_INFINITY:
                Iy_T[i][j] = ["M", "Iy"][Iy_vars.index(Iy[i][j])]
            ###### END IX-STEP
            
            
    return M, Ix, Iy, M_T, Ix_T, Iy_T
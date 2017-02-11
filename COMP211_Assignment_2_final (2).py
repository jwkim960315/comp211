# COMP211_Assignment_2
# Sean Kim, Brian Oh, Kim Pham
# Due Nov 8 2016


#1
def num_alignments(n,m):
    """
    Input: lengths of two strings
    Output: total number of possible alignments
    """
    memo = {}
    for i in range(0,m+1):
        memo[(i,0)] = 1
    for j in range(0,n+1):
        memo[(0,j)] = 1
    for i in range(1,m+1):
        for j in range(1,n+1):
            memo[(i,j)] = memo[(i-1,j-1)] + memo[(i,j-1)] + memo[(i-1,j)]
    return memo[(i,j)]

#2

MATCH = 5
MISMATCH = -2
GAP = -6

def seq_align_score_matrix(str1,str2):
    """
    Input: two DNA sequences in strings
    Output: a dictionary (representing a score matrix)
    """
    memo = {}
    memo[(0,0)] = 0
    for i in range(1,len(str1)+1):  
        memo[(i,0)] = i*GAP             # creating the first column
    for j in range(1,len(str2)+1):
        memo[(0,j)] = j*GAP             # creating the first row
    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if str1[i-1] == str2[j-1]:  
                memo[(i,j)] = memo[(i-1,j-1)] + MATCH       # adding MATCH from left diagonal score
            else:
                memo[(i,j)]=max((memo[(i-1,j-1)] + MISMATCH,        
                                   memo[(i,j-1)] + GAP,
                                   memo[(i-1,j)] + GAP))
    return memo

def seq_align_score(str1,str2):
    """
    Input: two DNA sequences in strings
    Output: the score of an optimal alignment
    """
    return seq_align_score_matrix(str1,str2)[(len(str1),len(str2))]



#3 
def possible_alignment(str1,str2):
    """
    Input: two DNA sequences in strings
    Output: an output of possible_alignment_hlpr,
            which is an alignmentof the best-score
            -sequences
    """
    memo = seq_align_score_matrix(str1,str2)
    str1 = str1[::-1]       # reverse string 1
    str2 = str2[::-1]       # reverse string 2
    return possible_alignment_hlpr(str1,str2,memo)

def possible_alignment_hlpr(str1,str2,memo):
    """
    Input: - two DNA sequences in strings
           - memo: a dictionary (representing score matrix)
    Output: a list containing sequence alignments
    """
    if str1 == '' and str2 == '':
        return []
    elif str1 == '':
        return possible_alignment_hlpr(str1,str2[1:],memo) + [('_',str2[0])]       # moving along columns(leftward) of the first row
    elif str2 == '':
        return possible_alignment_hlpr(str1[1:],str2,memo) + [(str1[0],'_')]       # moving along rows(upward) of the first column
    elif str1[0] == str2[0]:
        return possible_alignment_hlpr(str1[1:],str2[1:],memo) + [(str1[0],str2[0])]       # reduce matrix(1row,1column) & adding current string pairs
    elif memo[(len(str1),len(str2))] == memo[(len(str1)-1,len(str2)-1)] + MISMATCH:
        return possible_alignment_hlpr(str1[1:],str2[1:],memo) + [(str1[0],str2[0])]       # reduce matrix(1row,1column) & adding current string pairs
    elif memo[(len(str1),len(str2))] == memo[(len(str1),len(str2)-1)] + GAP:
        return possible_alignment_hlpr(str1,str2[1:],memo) + [('_',str2[0])]       # reduce matrix(1column) & inserting gap for str1
    elif memo[(len(str1),len(str2))] == memo[(len(str1)-1,len(str2))] + GAP:
        return possible_alignment_hlpr(str1[1:],str2,memo) + [(str1[0],'_')]       # reduce matrix(1row) & inserting gap for str2

def print_possible_alignment(alignment):
    """
    Input: an alignment of two DNA sequences
    Output: prints two separate sequences
    """
    seq1 = ''
    seq2 = ''
    for tup in alignment:
        seq1 += tup[0]      # adds first element from each base pairs
        seq2 += tup[1]      # adds second element from each base pairs
    print seq1
    print seq2
        


#4

def alignment_score(alignment):
    """
    Input: a possible alignment
    Output: a score for possible alignment
    """
    score = 0
    for tup in alignment:
        if tup[0] != tup[1]:
            if tup[0] == '_' or tup[1] == '_':
                score += GAP
            else:
                score += MISMATCH
        if tup[0] == tup[1]:
            score += MATCH
    return score
        

# Edit Distance for checking #2

def DPED2(str1,str2):
    '''Returns edit distance between str1 and str2 - linear space
       dynamic programming version'''

    memo = []
    for i in range(len(str1)+1):
        memo = memo + [i]
    memo2 = [1]
    for i in range(1,len(str1)+1):
        memo2 = memo2 + [0]

    for i in range(1,len(str2)+1):
        memo2[0] = i
        for j in range(1,len(str1)+1):
            if str1[j-1] == str2[i-1]:
                memo2[j] = memo[j-1]
            else:
                memo2[j] = 1 + min((memo[j-1],memo[j],memo2[j-1]))
        for k in range(len(str1)+1):
            memo[k] = memo2[k]

    return memo[len(str1)]


def test_num_alignments():
    """
    Input: None
    Output: prints a statement that "num_alignments" function passed all tests
    """
    assert num_alignments(0,0) == 1
    assert num_alignments(0,1) == 1
    assert num_alignments(1,0) == 1
    assert num_alignments(1,2) == 5
    assert num_alignments(2,1) == 5
    assert num_alignments(20,20) == 260543813797441
    print "Passed all tests"
    

def test_seq_align_score():
    """
    Input: None
    Output: prints a statement that "seq_align_score" function passed all tests
    """
    assert seq_align_score('ACGT','CGTA') == 3
    print "Passed all tests"

def test_possible_alignment():
    """
    Input: None
    Output: prints a statement that "possible_alignment" function passed all tests
    """
    assert possible_alignment('ACGT','CGTA') == [('A', '_'), ('C', 'C'), ('G', 'G'), ('T', 'T'), ('_', 'A')]
    assert possible_alignment('','') == []
    assert possible_alignment('CGTA','') == [('C', '_'), ('G', '_'), ('T', '_'), ('A', '_')]
    assert possible_alignment('','GCAT') == [('_', 'G'), ('_', 'C'), ('_', 'A'), ('_', 'T')]
    assert possible_alignment('AAAC','AGC') == [('A', '_'), ('A', 'A'), ('A', 'G'), ('C', 'C')]
    assert possible_alignment('TTCATA','TGCTCGTA') == [('T', 'T'), ('_', 'G'), ('_', 'C'), ('T', 'T'), ('C', 'C'), ('A', 'G'), ('T', 'T'), ('A', 'A')]
    print "Passed all tests"
    
def test_seq_align_score():
    """
    Input: None
    Output: prints a statement that "seq_align_score" function passed all tests
    """
    global MATCH,MISMATCH,GAP
    MATCH = 1
    MISMATCH = -1
    GAP = -2
    assert seq_align_score('ACGT','CGTA') == -1
    MATCH = 5
    MISMATCH = -2
    GAP = -6
    assert seq_align_score('ACGT','CGTA') == 3
    MATCH = 0
    MISMATCH = -1
    GAP = -1
    assert seq_align_score('ACGT','CGTA') == -DPED2('ACGT','CGTA')
    print "Passed all tests"


def test_alignment_score():
    """
    Input: None
    Output: prints a statement that "alignment_score" function passed all tests
    -Assuming MATCH,MISMATCH,and GAP are defined by solution #2
    """
    global MATCH,MISMATCH,GAP
    MATCH = 5
    MISMATCH = -2
    GAP = -6
    assert alignment_score(possible_alignment('a','_')) == -6
    assert alignment_score(possible_alignment('a','b')) == -2
    assert alignment_score(possible_alignment('a','a')) == 5
    assert alignment_score(possible_alignment('_','b')) == -6
    assert alignment_score(possible_alignment('abc','cba')) == 1
    assert alignment_score(possible_alignment('ab','bc')) == -4
    print "Passed all tests"


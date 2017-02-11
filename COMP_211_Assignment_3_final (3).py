#Sean Kim, Kim Pham, Brian Oh
#COMP 211 Assignment 3
#Nov 17 2016

# Global Variables
MATCH = 5
MISMATCH = -2
GAP = -6


#1
def openFASTA():
    """
    Input: None
    Output: returns a list of tuples, each containing a name and a sequence
    """
    tup_lst = []
    elem_1 = ''
    elem_2 = ''
    f = open('mitochondrial_dna.fasta','r')
    for line in f:
        if '>' in line:
            tup_lst.append((elem_1,elem_2))
            elem_1,elem_2 = '',''
            elem_1 += line[1:].strip('\n')
        else:
            elem_2 += line.strip('\n')
    tup_lst.append((elem_1,elem_2))
    return tup_lst[1:]

#2
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


def dist(seq1,seq2):
    """
    Input: two DNA sequences in strings
    Output: the score of two DNA sequences
    """
    return (len(seq1)*MATCH + len(seq2)*MATCH)/2. - seq_align_score(seq1,seq2)


def matrix_of_distance(tup_lst):
    '''
    Input: list of tuples of the form (name,sequence)
    Output: matrix of distances between the species
    '''
    lst = []
    for row_tup in tup_lst:
        sub_lst = []
        for colmn_tup in tup_lst:
            sub_lst += [dist(row_tup[1],colmn_tup[1])]
        lst += [sub_lst]
    return lst


def test_matrix_of_distance():
    """
    Input: None
    Output: checks all the cases below for matrix_of_distance and prints
            "passed all tests"
    """
    assert matrix_of_distance([('a',''),('b','')]) == [[0.0, 0.0], [0.0, 0.0]]
    assert matrix_of_distance([('a','a'),('b','')]) == [[0.0, 8.5], [8.5, 0.0]]
    assert matrix_of_distance([('a',''),('b','a')]) == [[0.0, 8.5], [8.5, 0.0]]
    assert matrix_of_distance([('a','A'),('b','A')])==[[0.0, 0.0], [0.0, 0.0]]
    assert matrix_of_distance([('a','AAAA'),('b','AAAA'),('c','AAAA')])==[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    assert matrix_of_distance([('a','A'),('b','G')])==[[0.0, 7.0], [7.0, 0.0]]
    assert matrix_of_distance([('a','G'),('b','A')])==[[0.0, 7.0], [7.0, 0.0]]
    assert matrix_of_distance([('a','A'),('b','G'),('c','C')])==[[0.0, 7.0, 7.0], [7.0, 0.0, 7.0], [7.0, 7.0, 0.0]]
    assert matrix_of_distance([('a','C'),('b','G'),('c','A')])==[[0.0, 7.0, 7.0], [7.0, 0.0, 7.0], [7.0, 7.0, 0.0]]
    assert matrix_of_distance([('a','AG'),('b','A')])==[[0.0, 8.5], [8.5, 0.0]]
    assert matrix_of_distance([('a','A'),('b','AG')])==[[0.0, 8.5], [8.5, 0.0]]
    assert matrix_of_distance([('a','AGT'),('b','TGA')])==[[0.0, 14.0], [14.0, 0.0]]
    assert matrix_of_distance([('a','AGTC'),('b','CTGA')])==[[0.0, 28.0], [28.0, 0.0]]
    assert matrix_of_distance([('a','AAAA'),('b','GGGG'),('c','CCCC'),('d','AGCT')])==[[0.0, 28.0, 28.0, 21.0], [28.0, 0.0, 28.0, 21.0], [28.0, 28.0, 0.0, 21.0], [21.0, 21.0, 21.0, 0.0]]
    print "passed all tests"


#3

def average_distance(tup_lst):
    '''
    Input: list of tuples of the form (name,sequence)
    Output: the average of all entries to right (or left) of the 0 diagonal
    '''
    lst = matrix_of_distance(tup_lst)
    sumlst = 0
    for sublst in lst:
        sumlst += sum(sublst)
    sumlst /= ((len(lst)**2.)-len(tup_lst))
    return sumlst

def test_average_distance():
    """
    Input: None
    Output: checks all the cases below for average_distance and prints
            "passed all tests"
    """
    assert average_distance([('a',''),('b','')])==0.0
    assert average_distance([('a','a'),('b','')])==8.5
    assert average_distance([('a',''),('b','a')])==8.5
    assert average_distance([('a','A'),('b','A')])==0.0
    assert average_distance([('a','AAAA'),('b','AAAA'),('c','AAAA')])==0.0
    assert average_distance([('a','A'),('b','G')])==7.0
    assert average_distance([('a','G'),('b','A')])==7.0
    assert average_distance([('a','A'),('b','G'),('c','C')])==7.0
    assert average_distance([('a','C'),('b','G'),('c','A')])==7.0
    assert average_distance([('a','AG'),('b','A')])==8.5
    assert average_distance([('a','A'),('b','AG')])==8.5
    assert average_distance([('a','AGT'),('b','TGA')])==14.0
    assert average_distance([('a','AGTC'),('b','CTGA')])==28.0
    assert average_distance([('a','AAAA'),('b','GGGG'),('c','CCCC'),('d','AGCT')])==24.5
    print "passed all tests"
    
              
#4
import random

def random_seq_lst(tup_lst):
    """
    Input: a list of tuples, each containing a name and a sequence
    Output: a list of tuples with randomized sequences
    """
    seq = ''
    for i in range(len(tup_lst)):
        while len(seq) != len(tup_lst[i][1]):
            seq += random.choice('ACGT')
        tup_lst[i] = (tup_lst[i][0],seq)
        seq = ''
    return tup_lst

def average_random_seq_lst():
    '''
    Input: None
    Output: the average of all entries to right (or left) of the 0 diagonal
    '''
    tup_lst = random_seq_lst(openFASTA())
    return average_distance(tup_lst)



#5
def random_seq_lst_prob(tup_lst):
    """
    Input: a list of tuples, each containing a name and a sequence
    Output: a list of tuples with randomized, but identical distributions
            as the real sequences, sequences
    """
    lst = []
    for i in range(len(tup_lst)):
        tup_into_lst = list(tup_lst[i][1])
        random.shuffle(tup_into_lst)
        lst += [(tup_lst[i][0],''.join(tup_into_lst))]
    return lst

#6
def average_distance_random_prob():
    """
    Input: None
    Output: the average of all entries to right (or left) of the 0 diagonal
    """
    tup_lst = random_seq_lst_prob(openFASTA())
    return average_distance(tup_lst)

#7 (Comparison)
def main():
    """
    Input: None
    Output: prints three distances (real, unequal distribution, equal
            distribution)
    """
    print 'avg. distance of real sequences: ' + str(average_distance(openFASTA()))
    print 'avg. distance of randomly distributed sequences: ' + str(average_random_seq_lst())
    print 'avg. distance of randomly mixed but equal distribution sequences: ' + str(average_distance_random_prob())

#7
"""
The average distance of real sequences is lower than those of equally
and unequally distributed sequences. From this result, we could
conclude that on average, the species from the real set are
further apart in time than those from other sequences. We think that the real
sequences are related in some way because our average score of 1342 from
the real set could be a smaller number, which would imply that the sequences with
an average score of 1342 are closer than sequences with an average score of 50.
"""

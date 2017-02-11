#Sean Kim, Kim Pham, Brian Oh
#COMP_211_Assignment_4
#Dec 1 2016

#Global Variables
MATCH = 5
MISMATCH = -2
GAP = -6

#import from assignment#3

def seq_align_scoring_matrix(str1,str2):
    
    '''Returns the dynamic programming scoring matrix for the global
       alignment of str1 and str2 assuming the constant costs associated
       with MATCH, MISMATCH, and GAP are defined as global variables.'''

    memo = {}
    memo[(0,0)] = 0
    for i in range(1,len(str1)+1):
        memo[(i,0)] = i*GAP
    for j in range(1,len(str2)+1):
        memo[(0,j)] = j*GAP
    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if str1[i-1] == str2[j-1]:
                match_mismatch = MATCH
            else:
                match_mismatch = MISMATCH
            memo[(i,j)]=max((memo[(i-1,j-1)]+match_mismatch,
                             memo[(i,j-1)]+GAP,
                             memo[(i-1,j)]+GAP))
    return memo


def seq_align_score(str1,str2):
    
    '''Returns the maximum alignment score for str1 and str2 assuming
       the constant costs associated with MATCH, MISMATCH, and GAP are
       defined as global variables.'''
    
    return seq_align_scoring_matrix(str1,str2)[(len(str1),len(str2))]

##def dist(seq1,seq2):
##    """
##    Input: two DNA sequences in strings
##    Output: the score of two DNA sequences
##    """
##    return (len(seq1)*MATCH + len(seq2)*MATCH)/2. - seq_align_score(seq1,seq2)

def compute_distance_matrix(seq_list):
    
    '''Takes a list of tuples (name,sequence) and returns
       a symmetric matrix of distances (derived from scores) for
       each pair of sequences. Note: diagonal is all zero.'''

    n = len(seq_list)
    dist_matrix = [] #Initialize as nxn zero matrix
    for i in range(n):
        dist_matrix.append(n*[0.0])
    for i in range(n-1):
        for j in range(i+1,n):
            dist_matrix[i][j] = ((MATCH/2.)*(len(seq_list[i][1]+seq_list[j][1])) -
                                 seq_align_score(seq_list[i][1],seq_list[j][1]))
            dist_matrix[j][i] = dist_matrix[i][j]
    return dist_matrix


def get_distances(matrix):
    
    '''Returns a list of the distances found in a symmetric distance matrix.'''

    n = len(matrix)
    scores = []
    for i in range(n-1):
        for j in range(i+1,n):
            scores.append(matrix[i][j])
    return scores

def dist_avg(distances):
    
    '''Returns the average of a list of distances'''
    
    return sum(distances)/(1.0*len(distances))


#1
class FullBinaryTree(object):
    '''Implements a full binary tree; each node has exactly two children,
       left and right. For internal nodes left and right are other internal
       nodes. For leaves, they are both None. Tree must contain at least
       one node.'''

    def __init__(self,left=None,right=None,parent=None):

        '''Constructor: left,right,parent are trees;
           default creates tree with a single node'''

        self.left = left
        self.right = right
        self.parent = parent
    
    def is_leaf(self):

        '''Returns True if Tree is a leaf'''

        return not self.left and not self.right

    def size(self):

        '''Returns the size (number of nodes) of tree'''

        if self.is_leaf():
            return 1
        else:
            return 1 + self.left.size() + self.right.size()

    def height(self):

        '''Returns the height (longest root to leaf path) of tree'''

        if self.is_leaf():
            return 0
        else:
            return 1 + max((self.left.height(),self.right.height()))

    def get_parent(self):
        
        '''
        Returns a parent tree object of the tree
        '''
        
        return self.parent

    def set_parent(self,new_parent):
        
        '''
        Given a tree object, the function sets that tree object as a parent tree
        '''
        self.parent = new_parent

    def is_root(self):
        '''
        Returns True if the tree is a root;
        Otherwise, returns False
        '''
        if not self.parent:
            return True
        return False


#2

class PhyloTree(FullBinaryTree):

    def __init__(self,left=None,right=None,parent=None,name='random_species',time=0.0):
        '''Constructor: left,right,parent are trees;
                        time is the passed time from the point where it splited
                        name is the species name
                        
           
           default creates tree with a single node'''
        FullBinaryTree.__init__(self,left,right,parent)
        self.name = name
        self.time = time
        if not self.is_leaf():
            self.left.set_parent(self)
            self.right.set_parent(self)
    #a)
    def __str__(self):
        
        '''
        Returns a Newick Representation of a tree
        '''
        
        if self.is_leaf():
            return str(self.name)
        else:
            return '('+str(self.left)+','+str(self.right)+')'
    #b)
    def parent_list(self):
        
        '''
        Returns a list of all the parents of the tree
        '''
        
        if self.parent == None:
            return [self]
        else:
            return [self] + self.parent.parent_list()
    
    
    def lca(self,other):

        '''
        Returns a least common ancestor of two trees
        '''
        
        self_lst = self.parent_list()
        other_lst = other.parent_list()
        for parent in self_lst:
            if parent in other_lst:
                return parent
        return None

    #c)
    def list_of_leaves(self):

        '''
        Returns a list of leaves of the tree
        '''
        
        if self.is_leaf():
            return [self]
        else:
            return self.left.list_of_leaves()+self.right.list_of_leaves()
    #d)
    def get_time(self):

        '''
        Returns the passed time from the point it splited
        '''
        
        return self.time
        
    #e)
    def get_species(self,name):
        '''
        Given a name of the tree, search for that tree;
        Return True if that tree exists in the current tree;
        Return None otherwise
        '''
        if self.name == name:
            return self
        elif self.left == None and self.right == None:
            return None
        else:
            return (self.left.get_species(name) or
                   self.right.get_species(name))
#3
def UPGMA(lst,matrix):
    '''
    Input: lst = list of names
           matrix = list of lists containing distances
    Output: a root of all the PhyloTrees
    '''
    
    lst_trees = []
    for species in lst:
        lst_trees.append(PhyloTree(name=species))
    return reduce_matrix([],[],matrix,lst_trees)


def reduce_matrix(matrix_hlpr,new_matrix,matrix,lst_trees):
    '''
    Input: matrix_hlpr == list used to append to new_matrix (initially empty)
           new_matrix == list that eventually become a matrix (initially empty)
           matrix == an original distance matrix
           lst_trees == a list of tree objects
    Output: returns a root of all the PhyloTrees
    '''
    if len(matrix) == 2: ## [[0,#],[#,0]]
        print matrix
        return PhyloTree(lst_trees[0],lst_trees[1],name='root',time=matrix[0][1]/2)
    else:
        lst_trees,row_num,column_num,min_distance = update_tree(matrix,lst_trees) #Call update_tree(matrix,lst_trees)
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if (i != row_num and i != column_num) and (j != row_num and j != column_num): #entries not in two rows and two columns
                    matrix_hlpr.append(matrix[i][j])
                    if len(matrix_hlpr) == len(matrix)-2: #new row has been filled
                        new_matrix.append(matrix_hlpr)
                        matrix_hlpr = []
        new_matrix = make_new_matrix(matrix,new_matrix,lst_trees[row_num],lst_trees[column_num],row_num,column_num)
        lst_trees.remove(lst_trees[row_num])
        lst_trees.remove(lst_trees[column_num-1])
        matrix = new_matrix[:]
        return reduce_matrix([],[],matrix,lst_trees)


def update_tree(matrix,lst_trees):
    '''
    Input: distance matrix
           list of tree objects
    Ouptut: an updated list of trees
            row number of a minimum value
            column number of a minimum value
            the distance between two tree objects
    '''
    min_value,row_num,column_num = min_matrix(matrix) #call min_matrix(matrix)
    min_distance = min_value/2.
    lst_trees.append(make_parent(lst_trees,row_num,column_num,min_distance))
    return lst_trees,row_num,column_num,min_distance

def min_matrix(matrix): 
    '''
    Input: a distance matrix
    Output: a minimum entry value of the given matrix
    '''
    min_matrix = matrix[0][1] #make the minimum just an arbitrary value for now
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] < min_matrix and matrix[i][j] != 0:
                min_matrix = matrix[i][j]
                row_num = i #store the row index for minimum entry
                column_num = j #store the column index for minimum entry
    return min_matrix,row_num,column_num

def make_parent(lst_trees,row_num,column_num,min_distance):
    '''
    Input: list of tree objects
           a row number of the mininum value
           a column number of the minimum value
           a distance between two tree objects
    Output: a parent tree object of the two tree objects
    '''
    return PhyloTree(lst_trees[row_num],lst_trees[column_num],
                     name='parent of '+lst_trees[row_num].name+' and '+lst_trees[column_num].name,
                     time=min_distance)



def make_new_matrix(matrix,new_matrix,tree1,tree2,row_num,column_num):
    '''
    Input: matrix == a previous matrix
           new_matrix == a new matrix that lacks a row and a column of the minimum value
           tree1,tree2 == two tree objects that are going to be joined
           row_num == a row number for the minimum value from a previous matrix
           column_num == a column number for the minimum value from a previous matrix
    Output: a computed-reduced matrix with joined row and column
    '''
    matrix_hlpr = []
    matrix_dist = make_dist_entry(row_num,column_num,matrix)
    for i in range(len(new_matrix)):
        for j in range(len(new_matrix)):
            if j == len(matrix)-3:
                new_entry = new_entry_comp(tree1,tree2,matrix_dist[0][i],matrix_dist[1][i])
                new_matrix[i].append(new_entry)
                matrix_hlpr.append(new_entry)
    new_matrix.append(matrix_hlpr+[0])
    return new_matrix
                

def make_dist_entry(row_num,column_num,matrix):
    
    '''
    Input: row_num == a row number for the minimum value from a previous matrix
           column_num == a column number for the minimum value from a previous matrix
           matrix == a previous matrix
    Output: return a matrix of distances between each of the two closest pair and the rest
            of the species
    '''
    
    matrix_dist = []
    matrix_hlpr = []
    index = row_num
    counter = 0
    while index == column_num or index == row_num:
        for j in range(len(matrix)):
            if j != row_num and j != column_num:
                matrix_hlpr.append(matrix[index][j])
        matrix_dist.append(matrix_hlpr)
        matrix_hlpr = []
        index = column_num
        counter += 1
        if counter == 2:
            index = 'stop'
    return matrix_dist


    
def new_entry_comp(tree1,tree2,dist1,dist2):
    '''
    Input: two tree objects
           two distances between two closest pair and a selected species
    Output: returns the average distance between the union of two closest pair and a selected species
    '''
    return (len(tree1.list_of_leaves())*float(dist1)+len(tree2.list_of_leaves())*dist2)/(len(tree1.list_of_leaves())+len(tree2.list_of_leaves()))


#Test for UPGMA (#4 and #5)

#a) fasta reading function (returns a list of names)
def read_fasta_names(fasta_file):
    """
    Input: fasta_file
    Output: returns a list of names
    """
    lst = []
    f = open(fasta_file,'r')
    for line in f:
        if '>' in line:
            lst.append(line[1:].strip('\n'))
    f.close()
    return lst

#b) fasta reading function (returns a list of tuples)
def read_fasta(fasta_file):
    
    '''Reads in a file in fasta format and returns list of tuples
       (name,sequence)'''

    sequence_list = []
    file_in = open(fasta_file, 'r')
    line = file_in.readline().strip()
    while line != '':
        if line[0] == '>':
            name = line[1:]
            seq = ''
            line = file_in.readline().strip()
            while line != '' and line[0] != '>':
                seq = seq + line
                line = file_in.readline().strip()
            sequence_list.append((name,seq))
        else:
            print "Error occurred while reading file."
    file_in.close()
    return sequence_list

lst = read_fasta_names('mitochondrial_dna.fasta')
tup_lst = read_fasta('mitochondrial_dna.fasta')
matrix = compute_distance_matrix(tup_lst)
root_tree = UPGMA(lst,matrix)

#4

def main():
    '''
    Input: None
    Output: Prints the answers to #4
    '''
    #a)
    print root_tree     #((Frog,(Carp,Loach)),(Chicken,(Human,((Mouse,Rat),(Seal,(Cow,Whale))))))

    #b)
    lca_human_cow = root_tree.get_species('Human').lca(root_tree.get_species('Cow'))
    print lca_human_cow         #(Human,((Mouse,Rat),(Seal,(Cow,Whale))))

    #c)
    print lca_human_cow.get_time()          #685.4

    #d)
    '''
    Overall the resulting tree represents a reasonable reconstruction of the phylogeny of
    the given species because one side includes mammals, while the other side includes
    aquatic animals (carp,loach,frog). It makes sense that carps and loaches are closer
    because they are freshwater fish, while frogs are branched off because they're
    amphibians. Additionally, it makes sense that mice and rats are closer to each other
    and to humans because they're land animals. What doesn't make sense is that the cows
    are closer to whales than whales are to seals. We expected the cows to have the same
    LCA as humans, rats, and mice.
    '''




#5
def ultrametric_matrix(root_tree):
    '''
    Input: root of all the PhyloTrees
    Output: a matrix of distances of an ultrametric tree
    '''
    UPGMA_matrix = []
    matrix_hlpr = []
    tree_lst = root_tree.list_of_leaves()
    for tree1 in tree_lst:
        for tree2 in tree_lst:
            if tree1 == tree2:
                matrix_hlpr.append(0)
            else:
                matrix_hlpr.append((tree1.lca(tree2).get_time())*2)
        UPGMA_matrix.append(matrix_hlpr)
        matrix_hlpr = []
    return UPGMA_matrix

ult_avg_dist = dist_avg(get_distances(ultrametric_matrix(root_tree)))
ori_avg_dist = dist_avg(get_distances(matrix))

print ult_avg_dist,ori_avg_dist         #1342.08888889 , 1342.08888889

#Comment
    '''
    We computed the sames averages for the orginal matrix and the UPGMA
    matrix. This is not surprising since the sum of the entries in the
    ultrametric matrix is the same as the original matrix. In UPGMA, averages
    are weighted by the number of sizes (ie.size(A)+size(B)), so each distance
    contributes equally to the result.

    '''





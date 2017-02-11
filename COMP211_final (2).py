#COMP211 Final
#Sean Kim
#Dec 17 2016


#import modules
import heapq, random, math, zlib, os, json, filecmp


#import FullBinaryTree class

class FullBinaryTree(object):

    '''Implements a full binary tree; each node should have exactly two children,
       left and right, and one parent. For interal nodes left and right are
       are other internal nodes. For leaves, the are both None. All nodes
       have a parent that is an internal node except the root whose parent
       is None. Tree must contain at least one node.'''

    def __init__(self,left=None,right=None,parent=None):

        '''Constructor creates a single node tree as default. Sets
           parent relation if left and right are given.'''

        self.left = left
        self.right = right
        self.parent = parent
        if self.left:
            self.left.set_parent(self)
        if self.right:
            self.right.set_parent(self)

    def set_parent(self,tree):

        '''Sets a given tree as a parent node'''

        self.parent = tree

    def get_parent(self):

        '''Returns a parent node'''

        return self.parent

    def is_leaf(self):

        '''Returns true iff node is a leaf'''

        return not self.left and not self.right

    def is_root(self):

        '''Returns true iff node is the root'''

        return not self.parent

    def size(self):

        '''Returns the size of the tree'''

        if self.is_leaf():
            return 1
        else:
            return 1 + self.left.size() + self.right.size()

    def height(self):

        '''Returns the height of the tree'''

        if self.is_leaf():
            return 0
        else:
            return 1 + max((self.left.height(),self.right.height()))

    def lca(self,tree):

        '''Returns the least common answer of self and tree'''

        my_anc = self.list_of_ancestors()
        tree_anc = tree.list_of_ancestors()
        i=0
        while  i<len(my_anc) and i<len(tree_anc) and my_anc[i] == tree_anc[i]:
            i = i+1
        if my_anc[i-1] == tree_anc[i-1]:
            return my_anc[i-1]
        else:
            return None


    def contains(self,tree):

        '''Returns true iff self contains tree as a subtree'''

        if self == tree:
            return True
        elif self.is_leaf():
            return False
        else:
            return self.left.contains(tree) or self.right.contains(tree)

    def list_of_ancestors(self):
        '''Returns list of ancestors including self'''

        if self.is_root():
            return [self]
        else:
            return self.parent.list_of_ancestors() + [self]

    def list_of_leaves(self):

        '''Returns a list of all of the leaves of tree'''

        if self.is_leaf():
            return [self]
        else:
            return self.left.list_of_leaves()+self.right.list_of_leaves()


#1

class HuffmanTree(FullBinaryTree):
    
    def __init__(self,left=None,right=None,parent=None,symbol=None,prob=None,code=''):
        '''Constructor for a single node HuffmanTree as a default.
           Sets a parent relation if left and right are given.'''
        FullBinaryTree.__init__(self,left,right,parent)
        self.symbol = symbol
        self.prob = prob
        self.code = code

    #(a)
    def __cmp__(self,other):
        
        '''Compares two trees according to their probabilities'''
        
        return cmp(self.prob,other.prob)

    #(b)
    def get_codeword(self):
        
        '''Returns the binary string (made up of 0s and 1s)
           created by concatenating the code values on the path from
           the root to self (not including the root code value).'''
        
        if self.is_root():
            return ''
        else:
            return self.parent.get_codeword() + self.code

    #(c)
    def get_symbol(self,symbol):
        
        '''Returns the leaf node in the tree containing 
           the given symbol if such a leaf exists.'''
        
        for tree in self.list_of_leaves():
            if symbol == tree.symbol:
                return tree
        return None
    
    #Extra Methods
    def set_code(self,cw):
        '''Sets code of a tree with cw (codeword)'''
        self.code = cw
    
    def get_prob(self):
        '''Returns probability of a tree'''
        return self.prob
        
#2

#Compression
def encode(infile,outfile):
    
    '''Creates an outfile, which is a compressed version of infile'''
    
    ori_f = open(infile,'r') #original file
    comp_file = open(outfile,'w')
    content = ori_f.read()
    symbol_freq_d = make_symbol_freq_dic(content)
    symbol_codeword_d = huffman_code_dic(huffman_code_tree(symbol_freq_d),symbol_freq_d)
    json.dump(symbol_codeword_d,comp_file)
    comp_file.write("\x00")     #dictionary,content separator
    comp_file.write(binary2char(make_binary_str(content,symbol_codeword_d)))
    ori_f.close()
    comp_file.close()

def make_symbol_freq_dic(content):
    
    '''Input: a string that has been read from encode function
       Output: a dictionary in the form of {symbol:frequency}'''
       
    symbol_freq_d = {}
    if len(content) == 0:
        return {None:1}     # infile is empty
    else:
        for symbol in content:
            if not symbol in symbol_freq_d:
                symbol_freq_d[symbol] = 1
            else:
                symbol_freq_d[symbol] += 1
    return symbol_freq_d

#(a) pt.1
def huffman_code_tree(d):
    
    '''Input: a dictionary in the form of {symbol:frequency}
       Output: a list of a HuffmanTree root'''
    
    pq = []
    for symbol,prob in d.items():
        pq.append(HuffmanTree(symbol=symbol,
                                    prob=prob))
    while len(pq) != 0:
        if len(pq) == 1:
            return pq[0]
        else:
            tree1 = heapq.heappop(pq)
            tree2 = heapq.heappop(pq)
            tree1.set_code('0')
            tree2.set_code('1')
            parent = HuffmanTree(left=tree1,right=tree2,prob=tree1.get_prob()+tree2.get_prob())
            heapq.heappush(pq,parent)
#(a) pt.2
def huffman_code_dic(huffmanTree,symbol_freq_d):
    
    '''Input: Huffman Tree object, a dictionary in the form of {symbol:frequency}
       Output: dictionary in a form of {symbol:codeword}'''
        
    symbol_codeword_d = {}
    if huffmanTree.size() == 1:     #file is empty or one symbol
        return {symbol_freq_d.keys()[0]: '0'}    #Arbitrarily sets '0' for codeword
    else:
        for symbol in symbol_freq_d:
            symbol_codeword_d[symbol] = huffmanTree.get_symbol(symbol).get_codeword()
        return symbol_codeword_d

#(b)        
def make_binary_str(string,dic):
    
    '''Returns a string of prefix codes (binary numbers) according to the infile characters'''
    
    binary_str = ''
    for symbol in string:
        binary_str += dic[symbol]
    return binary_str


def binary2char(string):

    '''Returns character encoded version of a binary string.
       Note: padded to be divisible by 8 with pad length as first char.'''

    pad = 8 - len(string)%8
    string = string+pad*'0'
    out = str(pad)+''.join([chr(int(string[i:i+8],2))
                            for i in range(0,len(string),8)])
    return out   

#Decompression
def decode(infile,outfile):
    
    '''Given a compressed file, creates a decompressed file'''
    
    comp_file = open(infile,'r')
    decomp_file = open(outfile,'w')
    jsonobj = ''
    content = comp_file.read()
    while jsonobj[:-2:-1] != "\x00":     #json object + \x00
        jsonobj += content[0] 
        content = content[1:]
    symbol_codeword_dic = json.loads(jsonobj[:len(jsonobj)-1]) #returns only the dictionary
    decomp_file.write(binary_to_symbol(char2binary(content),symbol_codeword_dic))
    comp_file.close()
    decomp_file.close()


def char2binary(string):

    '''Returns binary string represented by a character string.
       Assumes first char represents number of pad bits.'''

    pad = int(string[0])
    out = ''.join([(10-len(bin(ord(char))))*'0' + bin(ord(char))[2:] for
                    char in string[1:]])
    return out[:-1*pad]


#(c)
def binary_to_symbol(binary_string,dic):
    
    '''Returns the original string of text (consisting
       of the symbols given by the codewords)'''
       
    codeword = ''
    original_text = ''
    dic_reverse = {value: key for key, value in dic.iteritems()}
    for binary in binary_string:
        codeword += binary
        if codeword in dic_reverse:
            original_text += dic_reverse[codeword]
            codeword = ''
    return original_text


#3
def create_f(n,alphabet,outfile):
    
    '''Creates a file consisting of n symbols chosen 
       uniformly at random from a given alphabet'''
       
    f = open(outfile,'w')
    string = ''
    while len(string) != n:
        string += random.choice(alphabet)
    f.write(string)
    f.close()


#4
def cmp_files(file1,file2):
    
    '''Compares two files and returns true if they have identical
       content. Otherwise, it returns False'''
    
    return filecmp.cmp(file1,file2)


#5
def calc_entropy(content):
    
    '''A function that computes the entropy of a string from the relative
       frequencies of the symbols it contains'''
       
    symbol_freq_dic = make_symbol_freq_dic(content)
    entropy = 0.0
    total_freq = sum(symbol_freq_dic.values())
    for symbol,prob in symbol_freq_dic.items():
        prob /= float(total_freq)
        entropy += prob*math.log(1/(prob),2)
    return entropy



#Creaing outfile by using zlib
def zlib_f_creator(infile,outfile):
    
    '''Creates a file, which is a compressed version (by using zlib) of infile.
       Returns a content of an infile'''
    
    f = open(infile,'r')
    zlib_comp = open(outfile,'w')
    content = f.read()
    zlib_comp.write(zlib.compress(content))
    zlib_comp.close()
    f.close()
    return content

# Test Cases

def main():
    
    '''Prints out the result similar to the table in final.pdf'''
    
    #sample.txt
    encode('sample.txt','sample_comp.txt')
    decode('sample_comp.txt','sample_decomp.txt')
    content = zlib_f_creator('sample.txt','zlib_comp_sample.txt')
    print 'File: sample.txt'
    print 'Original Size: ' + str(os.path.getsize('sample.txt')) #46975
    print 'Compressed Size: ' + str(os.path.getsize('sample_comp.txt')) #30708
    print 'Fully recovered: ' + str(cmp_files('sample.txt','sample_decomp.txt')) #True
    print 'Zip Size: ' + str(os.path.getsize('zlib_comp_sample.txt')) #18968
    print 'Entropy: ' + str(calc_entropy(content)) #4.70968590567
    print '----------------------'

    
    #fasta.txt
    encode('fasta.txt','fasta_comp.txt')
    decode('fasta_comp.txt','fasta_decomp.txt')
    content = zlib_f_creator('fasta.txt','zlib_comp_fasta.txt')
    print 'File: fasta.txt'
    print 'Original Size: ' + str(os.path.getsize('fasta.txt')) #7060
    print 'Compressed Size: ' + str(os.path.getsize('fasta_comp.txt')) #2848
    print 'Fully recovered: ' + str(cmp_files('fasta.txt','fasta_decomp.txt')) #True
    print 'Zip Size: ' + str(os.path.getsize('zlib_comp_fasta.txt')) #2073
    print 'Entropy: ' + str(calc_entropy(content)) #2.11422675223
    print '----------------------'
    
    #random.txt
    create_f(20000,'abcdefghijklmnopqrstuvwxyz ','random.txt')
    encode('random.txt','random_comp.txt')
    decode('random_comp.txt','random_decomp.txt')
    content = zlib_f_creator('random.txt','zlib_comp_random.txt')
    print 'File: random.txt' 
    print 'Original Size: ' + str(os.path.getsize('random.txt')) #20000
    print 'Compressed Size: ' + str(os.path.getsize('random_comp.txt')) #12389 (varies a little everytime)
    print 'Fully recovered: ' + str(cmp_files('random.txt','random_decomp.txt')) #True
    print 'Zip Size: ' + str(os.path.getsize('zlib_comp_random.txt')) #12623 (varies a little everytime)
    print 'Entropy: ' + str(calc_entropy(content)) #4.75369482888
    print '----------------------'
    
    
    #allAs.txt
    create_f(20000,'A','allAs.txt')
    encode('allAs.txt','allAs_comp.txt')
    decode('allAs_comp.txt','allAs_decomp.txt')
    content = zlib_f_creator('allAs.txt','zlib_comp_allAs.txt')
    print 'File: allAs.txt'
    print 'Original Size: ' + str(os.path.getsize('allAs.txt')) #20000
    print 'Compressed Size: ' + str(os.path.getsize('allAs_comp.txt')) #2513
    print 'Fully recovered: ' + str(cmp_files('allAs.txt','allAs_decomp.txt')) #True
    print 'Zip Size: ' + str(os.path.getsize('zlib_comp_allAs.txt')) #43
    print 'Entropy: ' + str(calc_entropy(content)) #0.0
    print '----------------------'
    
    #empty.txt
    create_f(0,'','empty.txt')
    encode('empty.txt','empty_comp.txt')
    decode('empty_comp.txt','empty_decomp.txt')
    content = zlib_f_creator('empty.txt','zlib_comp_empty.txt')
    print 'File: empty.txt'
    print 'Original Size: ' + str(os.path.getsize('empty.txt')) #0
    print 'Compressed Size: ' + str(os.path.getsize('empty_comp.txt')) #16
    print 'Fully recovered: ' + str(cmp_files('empty.txt','empty_decomp.txt')) #True
    print 'Zip Size: ' + str(os.path.getsize('zlib_comp_empty.txt')) #8
    print 'Entropy: ' + str(calc_entropy(content)) #0.0
    print '----------------------'
    
#6

'''
According to the result, it looks like lower entropy means better compression 
(except empty.txt). I was surprised with the result at the end because I never 
realized that the compressed size could be greater than the original size. The 
result of my algorithm is worse than that of professor's algorithm except allAs.txt 
compression (barely better). I think the reason is that I dumped dictionary 
itself, but the professor might have a different way to store it, which takes 
up less byte.
'''



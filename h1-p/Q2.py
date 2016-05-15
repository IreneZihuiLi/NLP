
# coding: utf-8

# In[84]:

#! /usr/bin/python

__author__="Daniel Bauer <bauer@cs.columbia.edu>"
__date__ ="$Sep 12, 2011"

import sys
from collections import defaultdict
from operator import itemgetter
import math

"""
Count n-gram frequencies in a data file and write counts to
stdout. 
"""

def simple_conll_corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuples. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            fields = line.split(" ")
            ne_tag = fields[-1] # tag
            word = " ".join(fields[:-1])
            yield word, ne_tag
        else: # Empty line
            yield (None, None)                        
        l = corpus_file.readline()

def sentence_iterator(corpus_iterator):
    """
    Return an iterator object that yields one sentence at a time.
    Sentences are represented as lists of (word, ne_tag) tuples.
    """
    current_sentence = [] #Buffer for the current sentence
    for l in corpus_iterator:    
#             print ("***NOW",l) # Every line (every word) l -->  ('nucleotidase', 'I-GENE')
            if l==(None, None):
                if current_sentence:  #Reached the end of a sentence
                    yield current_sentence
                    current_sentence = [] #Reset buffer
                else: # Got empty input stream
                    sys.stderr.write("WARNING: Got empty input file/stream.\n")
                    raise StopIteration
            else:
                current_sentence.append(l) #Add token to the buffer

    if current_sentence: # If the last line was blank, we're done
#         print ("Current ",current_sentence)
        yield current_sentence  #Otherwise when there is no more token
                                # in the stream return the last sentence.
        #Current_sentence  [('When', 'O'), ('CSF', 'O'), ('joj', 'O'),...]  (read in as a key-value pairs collection)

def get_ngrams(sent_iterator, n):
    """
    Get a generator that returns n-grams over the entire corpus,
    respecting sentence boundaries and inserting boundary tokens.
    Sent_iterator is a generator object whose elements are lists
    of tokens. 每个句子前追加2个*，后追加一个STOP标志(w_boundary)
    """
    for sent in sent_iterator:
         # sent (list) --> [('omparison', 'O'), ('with', 'O'),  ...] (Tupels in every sentence)
         #Add boundary symbols to the sentence
         w_boundary = (n-1) * [(None, "*")]
         w_boundary.extend(sent)
         w_boundary.append((None, "STOP"))
         #Then extract n-grams
        
         ngrams = (tuple(w_boundary[i:i+n]) for i in range(len(w_boundary)-n+1))
         for n_gram in ngrams: #Return one n-gram at a time
#             print (n_gram)
            yield n_gram      

#         ngrams
#         (1, 2, 3)
#         (2, 3, 4)
#         (3, 4, 5)
#         (4, 5, 6)
#         (5, 6, 7)
#         >>> 类似这种
#  返回的是这样：
# ((None, '*'), (None, '*'), ('omparison', 'O'))
# ((None, '*'), ('omparison', 'O'), ('with', 'O'))
# (('omparison', 'O'), ('with', 'O'), ('alkaline', 'I-GENE'))


class Hmm(object):
    """
    Stores counts for n-grams and emissions. 
    """

    def __init__(self, n=3):
        assert n>=2, "Expecting n>=2."
        self.n = n
        self.emission_counts = defaultdict(int)
        self.emision_parameters = defaultdict(float)
        self.rare = list()

        self.ngram_counts = [defaultdict(int) for i in range(self.n)]
        self.all_states = set()
#       stores tragram q(C|A,B) = count (A,B,C)/count(B,C)   --> tuple (A,B,C):float (q)
        self.trigram_param = defaultdict(float)
        self.bigram = defaultdict(int)
        self.trigram = defaultdict(int)
#       rare word list
#         self.rare = list()
#       argmax word list :    word([tag,em],[tag,em],...)
#         self.argmax_list = defaultdict(list)
#         self.argmax_value = defaultdict(str)
        self.possible_tag = defaultdict(str)

   
    def write_counts(self, output, printngrams=[1,2,3]):
        """
        Writes counts to the output file object.
        Format:

        """
        # First write counts for emissions
        for word, ne_tag in self.emission_counts:            
            output.write("%i WORDTAG %s %s\n" % (self.emission_counts[(word, ne_tag)], ne_tag, word))


        # Then write counts for all ngrams
        for n in printngrams:            
            for ngram in self.ngram_counts[n-1]:
                ngramstr = " ".join(ngram)
                output.write("%i %i-GRAM %s\n" %(self.ngram_counts[n-1][ngram], n, ngramstr))

    def read_counts(self, corpusfile):

        self.n = 3
        self.emission_counts = defaultdict(int)
        self.ngram_counts = [defaultdict(int) for i in range(self.n)]
        self.all_states = set()
        self.bigram = defaultdict(int)
        self.trigram = defaultdict(int)

        for line in corpusfile:
            parts = line.strip().split(" ")
            #if this line is about '3-GRAM'
            if '3-GRAM' in parts :
                key_tuple = tuple(parts[2:])
                self.trigram [key_tuple] = int (parts[0])
            elif parts[1] == "WORDTAG":
                ne_tag = parts[2]
                word = parts[3]
                self.emission_counts[(word.lower(), ne_tag)] = int (parts[0])
                self.all_states.add(ne_tag)
            elif '2-GRAM' in parts :
                key_tuple = tuple(parts[2:])
                self.bigram [key_tuple] = int (parts[0])
            elif parts[1].endswith("GRAM"):
                n = int(parts[1].replace("-GRAM",""))
                ngram = tuple(parts[2:])
                self.ngram_counts[n-1][ngram] = int (parts[0])
                
        #print(self.trigram,'\n')
        #print(self.bigram,'\n')
        print (self.all_states)
    
    def get_emission_param(self):
        #('resuscitation', 'O'): 2.0
        self.emision_parameters = defaultdict(float)
        for wordNtag,counts in self.emission_counts.items():
            tag = wordNtag[1] #get tag of the tuple
            tagCount = self.ngram_counts[0][(tag,)]
            emission = counts / tagCount
            #change word to lower case
            new_word_n_tag = tuple((wordNtag[0].lower(),wordNtag[1]))
            self.emision_parameters[new_word_n_tag]=emission
            if 'BTK' in new_word_n_tag:
                print ('BTK',new_word_n_tag, emission)
            
        
        #print ('emi\n',self.emision_parameters,'\n')
    
    def get_trigram_param(self):
    
        self.trigram_param = defaultdict(float)
        for tag_tuple,counts in self.trigram.items():
            subtag = tuple(tag_tuple[1:]) #get tag of the tuple
            tagCount = self.bigram[subtag]
            if tagCount:
                q = counts / tagCount
            else: 
                # if 分母为0
                q = 0.
            self.trigram_param[tag_tuple]=q
        
        #print (self.trigram_param,'\n')
        
        

    def replaceRare(self,train,trainingChange_file):
        for wordNtag,counts in self.emission_counts.items():
            if counts < 5:
                self.rare.append(wordNtag[0])
                

        i = 0
        for line in train:
            i=i+1
            if i%10 ==0:
                print ('Now',i)
            l = line.strip()
            if l:
                fields = l.split(" ")
                ne_tag = fields[-1] # tag
                word = fields[0]
                if word in self.rare:
                    #newline = "%s %s\n" % ('_RARE_',ne_tag)
                    newline = "%s\n" % ('_RARE_')

                    trainingChange_file.write(newline)
                else: 
#                     print (l)
                    trainingChange_file.write(line)
            else:
                trainingChange_file.write('\n')
                
        
        trainingChange_file.close()


    def get_argmax(self):
#         buid argmax list
        for wordNtag,values in self.emision_parameters.items():
#             if wordNtag[0]=='_RARE_':
                
            self.argmax_list[wordNtag[0]].append((wordNtag[1],values))
#         build argmax real tag
        for word,values in self.argmax_list.items():
            tag = max(values,key=itemgetter(1))[0]
            self.argmax_value[word]=tag
            
        print(self.argmax_value)
        
        
        
    def write_prediction(self,testing,output):
        # read in word first
        for line in testing:
            word = line.strip()
            if word:
                tag=self.argmax_value[word]
                if tag:
                    output.write("%s %s\n" % (word, tag))
                else:
                    tag=self.argmax_value['_RARE_']
#                     print ('New word',tag)
                    output.write("%s %s\n" % (word, tag))
            else:
                output.write("\n")
        
                
#         # write (word tag) 
#         for word, tag in self.argmax_value.items():            
#             output.write("%s %s\n" % (word, tag))
        print ("Finish!")
        output.close()
        
        
        
        
        
    #Viterbi Algo         
    def viterbi(self,sentence):
        
        def findSet(index):
            if index in range(1,len(sentence)+1):
                return self.all_states
            elif index == 0 or index == -1:
                return {'*'}
                
        
        
        #S = self.all_states

        #stores (word:tag) in this whole sentence
        sentence_with_tag = defaultdict(str)
        
            
        #inner function to commpute pi values--start
        def pi_viterbi(k,u,v,sentence):
            prob = defaultdict(float)
            #print ('...getting...',k,u,v)
            #initial
            if k==0 and u == '*' and v == '*':
            #if k==0:
                return (1.,'*')
            else:
                for w in findSet(k-2):
                    prev = pi_viterbi(k-1,w,u,sentence)[0]
                    #tuple((w,u,v))
                    q = self.trigram_param[tuple((w,u,v))]
                    e = self.emision_parameters[tuple((sentence[k-1].lower(),v))]
                    probability = prev * q * e
                    prob [tuple((w,u))] = probability
                    #if probability==0:
                        #print ('Got-> prev ',k-1,w,u,':',prev,' q of ',v,w,u,':',q,' e of ',sentence[k-1],v,':',e,' so:',probability)
                max_tuple = max(prob.items(), key=lambda x: x[1])
                
                #print (max_tuple[1],max_tuple[0][0])
                return max_tuple[1],max_tuple[0][0]
            
        #inner function to commpute pi values--end
        
        
        sentence_with_tag= list()
        backpointer=defaultdict(str)
        tags = defaultdict(str)
        k = len(sentence)
        u_glob = ''
        v_glob = ''
        glob=0.
        #k=4
        for i in range(1,k+1):
            prob = defaultdict(float)
            for u in findSet(i-1):
                for v in findSet(i):
                    value,w=pi_viterbi(i,u,v,sentence)
                    prob [tuple((i,u,v))] = value
                    backpointer[tuple((i,u,v))]=w
            max_tuple = max(prob.items(), key=lambda x: x[1])
            #print ('Found',sentence[i-1],max_tuple)
            #if i == 1:
            
            backpointer[tuple((i,max_tuple[0][1],max_tuple[0][-1]))] = max_tuple[0][1] # bp (k,u,v)= tag w
            
            #sentence_with_tag.append(max_tuple[0][-1])
            u_glob = max_tuple[0][-2]
            v_glob = max_tuple[0][-1]
            glob = max_tuple[1]
            print ('Max',max_tuple)
        tags[k-1]=u_glob
        tags[k]=v_glob
        
        for i in range((k-2),0,-1):
            tag = backpointer[tuple(((i+2),tags[i+1],tags[i+2]))]
            tags[i]=tag
        
        tag_list=list()
        for i in range(1,len(tags)+1):
            tag_list.append(tags[i])
            
        #tag list as results
        return tag_list
        
        
    
            
            
    
    
#     for word,values in self.argmax_list.items():
#             tag = max(values,key=itemgetter(1))[0]
#             self.argmax_value[word]=tag
            
    
        
    def parse_sentence(self,testing,output):
        #read in sentences
        sentence = list()
        for line in testing:
            line = line.strip()
            if line:
                #get a sentence
                word = line.split(" ")
                if word:
                    sentence = sentence + word
            else:
                print ('Found a sentence!')
                
                
                #for each sentence
                #start Viterbi here
                
                tags=self.viterbi(sentence)

                
                #stop Viterbi here
                
                # write results to new file
                for index,value in enumerate(tags):
                    print ('wrote: ',(sentence[index], value))
                    output.write("%s %s\n" % (sentence[index], value))
                sentence.clear()
                
        output.close()
                
        #testing code 
        
        #sentence = ['BACKGROUND', ':', 'Ischemic', '_RARE_', '_RARE_', 'is', 'the', '.']
        #sentence = ['STAT5A', 'mutations', 'in', 'the', 'Src', 'homology', '2', '(', '_RARE_', ')', 'and', 'SH3', 'domains', 'did', 'the', '_RARE_', '-', '_RARE_', 'tyrosine', '_RARE_', '.']
        #self.viterbi(sentence)
    
  
            
            
    
        


def usage():
    print ("""
    python count_freqs.py [input_file] > [output_file]
        Read in a gene tagged training input file and produce counts.
    """)

if __name__ == "__main__":

    try:
        # Read in a count file
        input = open("gene.counts.replaced",'r')
        testing = open("mygenebig.dev.replace",'r')
        #newfile = open("mygenebig.dev.replace",'w')
        result = open("prediction.replace",'w')
        
    except IOError:
        sys.stderr.write("ERROR: Cannot read inputfile %s.\n" % arg)
        sys.exit(1)
    
    # Initialize a trigram counter
    counter = Hmm(3)
    
    #Q2: Read the counts
    counter.read_counts(input)
    counter.get_trigram_param()
    counter.get_emission_param()
    
    counter.parse_sentence(testing,result)
#     counter.replaceRare(testing,newfile)
    
    


# In[ ]:




# In[ ]:


# 47
# down vote
# accepted
# Use max():

 
# Using itemgetter():

# In [53]: lis=[(101, 153), (255, 827), (361, 961)]

# In [81]: from operator import itemgetter

# In [82]: max(lis,key=itemgetter(1))[0]    #faster solution
# Out[82]: 361


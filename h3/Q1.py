
# coding: utf-8

# In[16]:

import sys
from collections import defaultdict
import math,random,re

"""
IBM Model 1: estimate t parameters
EM Algorithm: finding the alignments
"""



class IBM_Model_1(object):
    """
    Stores counts for n-grams and emissions. 
    """

    def __init__(self, n=3):
        #sentences <--bad solution
        self.source = list()
        self.target = list()
        
        self.source_corpos = defaultdict(int)
        self.target_corpos = defaultdict(int)
        #counts
        self.source_target_occurance = defaultdict(int) # c(e,f)
        self.source_occurance = defaultdict(int) # c (e)
        self.conditioned_occurance = defaultdict(int) # c(j|i,l,m)
        self.joint_occurance = defaultdict(int) # c (i,l,m)
        #term-frenquency?
        self.t = defaultdict(float) # t (f|e)
        self.q = defaultdict(float)
        #delta
        self.delta = defaultdict(float)

    def counter(self,source_file,target_file):
        for line in source_file:
            line = line.strip().split(' ')
            
            for word in line:
                self.source_corpos[word] +=1
            # save only words and digis
            #line = [x for x in line if re.match("^[A-Za-z0-9_-]*$", x)]
            line.insert(0,'NULL')
            self.source.append(line)
        for line in target_file:
            line = line.strip().split(' ')
            #line = [x for x in line if re.match("^[A-Za-z0-9_-]*$", x)]
            self.target.append(line)
            for word in line:
                self.target_corpos[word] +=1
        print ('Load done!')
        print ('Found %i words in source'%len(self.source_corpos))
        print ('Found %i words in target'%len(self.target_corpos))
        print (len(self.source))
        print (len(self.target))

    #def sentences(self,source_file,target_file):
        #load in a pair of sentences at one time
    
    def em(self):
        
        # all counts are zeros (down by defaultdict), t (f|e ) are random values range in [.0,1.]
        for k in range(1, len(self.source)):
            for i in range (1, len(self.target[k])):
                for j in range (1,len(self.source[k])):
                    self.t[tuple((self.target[k][i],self.source[k][j]))] = 1/self.source_corpos[self.source[k][j]]

    
        # normally should be about 10-20
        for s in range(1,11):
            # iterate every sentence pairs
            for k in range(0, len(self.source)):
                m = len(self.target[k])
                l = len(self.source[k])
                for i in range (0, m):
                    down = 0.0
                    for j in range (0,l):
                        # 分母
                        down += self.t[tuple((self.target[k][i],self.source[k][j]))] 
                    for j in range (0,len(self.source[k])):
                        if down == 0.0:
                            self.delta [tuple((k,i,j))] = 0.0
                        else:
                            self.delta [tuple((k,i,j))] = self.t[tuple((self.target[k][i],self.source[k][j]))] /down
                        self.source_target_occurance [tuple((self.source[k][j],self.target[k][i]))] += self.delta [tuple((k,i,j))]
                        self.source_occurance [self.source[k][j]] += self.delta [tuple((k,i,j))]
                        #self.conditioned_occurance [tuple((j,i,l,m))] += self.delta [tuple((k,i,j))]
                        #self.joint_occurance [tuple((i,l,m))] += self.delta [tuple((k,i,j))]
            
            # Update t
            sums = 0
            for k,v in self.t.items():
                if self.source_target_occurance [k] * self.source_occurance[k[1]] >0 : 
                    result = self.source_target_occurance [k]/self.source_occurance[k[1]]
                    sums +=1
                else: 
                    result = 0.0
                self.t[k] = result
#             sumss = 0
#             for k,v in self.q.items():
#                 if self.conditioned_occurance [tuple((j,i,l,m))] * self.joint_occurance [tuple((i,l,m))] >0 : 
#                     result = self.conditioned_occurance [tuple((j,i,l,m))] / self.joint_occurance [tuple((i,l,m))]
#                     sumss +=1
#                 else: 
#                     result = 0.0
#                 self.q[tuple((j,i,l,m))] = result
            
            
            print ('Now..',s)
       
                    
    def write_t (self, output):
         
        for k,v in self.t.items():
            if v >0:
                output.write('%s %s %f\n'%(k[0],k[1],v))
    
        output.close()
        print ('Done!')
        
        
    def get_align(self, source, target, output):
        #sentences <--bad solution
        self.source = list()
        self.target = list()
        
        for line in source:
            line = line.strip().split(' ')
            
            #line.insert(0,'NULL')
            self.source.append(line)
        for line in target:
            line = line.strip().split(' ')
            self.target.append(line)
           
        
        
        for k in range(1, len(self.source)):
            m = len(self.target[k])
            l = len(self.source[k])
            # for each sentence
            for i in range (1,len(self.source)):
                source = self.source[i]
                target = self.target[i]
                
                # for each word in target sentence
                for (index,f) in enumerate(target):
                    align = list ()
                    # indicate which sentence
                    align.append(i)
                    align.append(index+1)
                    alignments = defaultdict(float)
                    for (indexe, e) in enumerate(source):
                        if self.t[tuple((f,e))] > 0:
                            alignments[indexe] = self.t[tuple((f,e))]
                    if len(alignments) == 0:
                        a = 0
                    else :
                        # find the max
                        a = max(alignments.items(), key=lambda a: a[1])[0]
                        align.append(a+1)
                        output.write('%i %i %i\n' %(align[0],align[1],align[2]))

                    
                    align.clear()
        output.close()
        print ('Finish Writing!')

def usage():
    print ("""
    
        Read in both English and Spanish files
    """)

if __name__ == "__main__":

    
    source_file = open (r'corpus.en','r')
    target_file = open (r'corpus.es','r')
    output_file = open ('small_output','w')
    output_key = open ('small_keys','w')
    source = open (r'dev.en','r')
    target = open (r'dev.es','r')
    # Initialize a trigram counter
    my_model = IBM_Model_1()
    my_model.counter(source_file,target_file)
    my_model.em()
    my_model.write_t(output_file)
    my_model.get_align(source,target,output_key)
    


# In[ ]:




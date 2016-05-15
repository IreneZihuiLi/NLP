
# coding: utf-8

# In[ ]:

#! /usr/bin/python

__author__="Alexander Rush <srush@csail.mit.edu>"
__date__ ="$Sep 12, 2012"

import sys, json
from collections import defaultdict

"""
Count rule frequencies in a binarized CFG.
Get emissions
Parse the sentence by using a CKY algo:
 
input:  What was the monetary value of the Nobel Peace Prize in 1989 ?
output: ["SBARQ", ["WHNP+PRON", "What"], ["SBARQ", ...["NP+NUM", "1989"]]]]]], [".", "?"]]]

By using pretty_print_tree.py [parsed file name] to show the result
"""
# 1. write as json file
# 2. algorithm


class Counts:
  def __init__(self):
    self.unary = {}
    self.binary = {}
    self.nonterm = {}
    # rare words counts and rare words set
    self.count_word = defaultdict(int)
    self.rare=set()
    self.binary_emission = defaultdict(float) # q(X->Y1,Y2)
    self.unary_emission = defaultdict(float)  # q(X->w)
    
    
    
    
  def show(self):
#     for symbol, count in self.nonterm.items():
#       print (count, "NONTERMINAL", symbol)  # Count(X)

#     for (sym, word), count in self.unary.items():
#       print (count, "UNARYRULE", sym, word.lower())  # Count(X->w)

    for (sym, y1, y2), count in self.binary.items():
      print (count, "BINARYRULE", sym, y1, y2)  # Count(X->Y1,Y2)
    

  def count(self, tree):
    """
    Count the frequencies of non-terminals and rules in the tree.
    """
    if isinstance(tree, str): return

    # Count the non-terminal symbol. 
    symbol = tree[0]
    self.nonterm.setdefault(symbol, 0)
    self.nonterm[symbol] += 1
    
    if len(tree) == 3:
      # It is a binary rule.
      y1, y2 = (tree[1][0], tree[2][0])
      key = (symbol, y1, y2)
      self.binary.setdefault(key, 0)
      self.binary[(symbol, y1, y2)] += 1
      
      # Recursively count the children.
      self.count(tree[1])
      self.count(tree[2])
    elif len(tree) == 2:
      # It is a unary rule.
      y1 = tree[1]
      key = (symbol, y1.lower())
      self.unary.setdefault(key, 0)
      self.unary[key] += 1

  def check_rare(self):
    for (sym, word), count in self.unary.items():
        self.count_word[word.lower()]+=count

    for (word ,time) in self.count_word.items():
        if time < 5:
            self.rare.add(word.lower())
    
    
    
  def replace(self,infile,outfile):
    output = open (outfile,'w')
    for l in open(infile):
      t = json.loads(l)
      adict = self.replace_word(t)
      #print (type(adict))
      json.dump(adict,output,separators=(',',','))
      output.write('\n')
    output.close()
    print ('Done')
     
    
  

  def replace_word(self, tree):
    """
    Count the frequencies of non-terminals and rules in the tree.
    """
    if isinstance(tree, str): return 

    # Count the non-terminal symbol. 
    symbol = tree[0]
    
    if len(tree) == 3:
      # It is a binary rule.
      y1, y2 = (tree[1][0], tree[2][0])
      # Recursively count the children.
      self.replace_word(tree[1])
      self.replace_word(tree[2])
    elif len(tree) == 2:
      # It is a unary rule.
      y1 = tree[1]
      for word in self.rare:
        if y1 == word:
          tree [1] = '_RARE_'
    return tree

  def replaceRare(self, oldfile,testing):
        print (len(self.rare))
        testing_file = open(testing,'w')
        old = open(oldfile,'r')
        i = 0
        for line in old:
            i=i+1
            if i%40 ==0:
                print ('Now',i)
            l = line.strip()
            if l:
                fields = l.split(" ")
                newline = list()
                for item in fields:
                    if item.lower() in self.rare:
                        newline.append('_RARE_')
                        print ('Found rare')
                    else: 
                        newline.append(item)
#                 for word in newline:
#                     testing_file.write(word)
#                     testing_file.write()
                testing_file.write(str(' '.join(newline)))
                testing_file.write('\n')

        testing_file.close()
    
  def get_emission(self):
    for (sym, word), count in self.unary.items():
      # Count(X->w)
      self.unary_emission[tuple((sym, word.lower()))] = count / self.nonterm[sym]
#       if word.lower() =='what':
#           print (sym,self.unary_emission[tuple((sym, word.lower()))])
      
    for (sym, y1, y2), count in self.binary.items():
      # Count(X->Y1,Y2)
      self.binary_emission[tuple((sym, y1, y2))] = count / self.nonterm[sym]
    
    #print(self.unary_emission)
    #print (self.binary_emission)
        
  
      
  def cky(self):
    
    
                         
    line = '_RARE_ world _RARE_'
    #line = 'How much did Manchester'
    sentence = line.split(' ')
    n=len(sentence)
    
    def pi(i,j,sym):
        if i==j:
            result= self.unary_emission[tuple((sym, sentence[i-1].lower()))]
            if result:
                #print (i,j,sym,sentence[i-1],result)
                return result
            else:
                return 0.0
        else:
            pi_values = list ()
            pi_values.append(0.0)
            for l in range(1,n):
                for i in range(1,n):
                    j=i+l
                    rules = list ()
                    # find all rules start with sym
                    for (rule,value) in self.binary_emission.items():
                        if rule[0] == sym:
                            rules.append(rule)
                    #print (len(rules))
                    for rule in rules:
                        for s in range(i,j):
                            product = self.binary_emission[rule] * pi(i,s,rule[1]) * pi(s+1,j,rule[2])
                            pi_values.append(product)
            #print (len(pi_values),max(pi_values))
            
            return max(pi_values)
                    
        
            
   
    print ('res:',pi(1,n,'NP'))


def main(parse_file,oldfile,newfile):
  counter = Counts() 
  for l in open(parse_file):
    t = json.loads(l)
    counter.count(t)
    
  #counter.show()
  #print ('wegwef')
  counter.check_rare()
  counter.get_emission()
  counter.cky()

  #counter.replaceRare(oldfile,newfile)

def usage():
    sys.stderr.write("""
    Usage: python count_cfg_freq.py [tree_file]
        Print the counts of a corpus of trees.\n""")

    
if __name__ == "__main__":
    sys.setrecursionlimit(1000) 

    parse_file = 'parse_train.dat'
    replace_counts ='cfg.replace.counts'
    replace_file ='replace.dat'
    old_file='parse_dev.dat'
    newfile='newfile.test'

    main (replace_file,old_file,newfile)
    # will generate cfg.replace.counts as new count
  


# In[ ]:




# In[ ]:




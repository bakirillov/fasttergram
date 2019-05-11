import numpy as np
from itertools import groupby


class KMerList():
    
    def __init__(self, fasta, ks):
        self.fasta = groupby(fasta.seq, key=lambda x: int(x != "N"))
        self.ks = ks
        
    def __iter__(self):
        return(self)
    
    def __next__(self):
        is_not_n, non_overlapping = self.fasta.__next__()
        if is_not_n:
            overlapping = KMerList.overlap(
                "".join(non_overlapping), self.ks
            )
        else:
            overlapping = self.__next__()
        return(overlapping)
    
    @staticmethod
    def overlap(s, ks):
        kmers = []
        for a in range(len(s)):
            kmers.append(s[a:a+np.random.choice(ks)])
        return(np.array(kmers[0:len(kmers)-min(ks)+1]))
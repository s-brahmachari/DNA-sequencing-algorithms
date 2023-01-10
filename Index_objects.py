import bisect
from typing import List, Tuple

class Index:
    def __init__(self, t: str, k: int):
        """ Create index from all substrings of size 'length' """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p: str) -> List[int]:
        """ Return index hits for first k-mer of P """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(pattern: str, text: str, index: Index) -> List[int]:
    k = index.k
    offsets = []
    for i in index.query(pattern):
        if pattern[k:] == text[i + k:i + len(pattern)]:
            offsets.append(i)
    return offsets

class SubseqIndex:
    """ Holds a subsequence index for a text T """

    def __init__(self, t: str, k: int, ival: int):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i + self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p: str) -> List[int]:
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def querySubseqIndex(pattern: str, text: str, index: SubseqIndex) -> List[int]:
    k = index.k
    offsets = []
    for i in index.query(pattern):
        if pattern == text[i:i + len(pattern)]:
            offsets.append(i)
    return offsets

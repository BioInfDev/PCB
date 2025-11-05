from collections.abc import Iterator
from itertools import groupby
from typing import List, Tuple

from scipy.stats import entropy

from . import seq2smth

class ShannonDetector:
    '''
    Identifies regions of low complexity in a protein sequence. 
    Algorithm:
        0. Slicing protein sequence to subsequences using the sliding window method
        1. Calculates the Shannon entropy of subsequences. 
        2. Finds subsequences with entropy lower than random — low complexity points (LCPs).
        3. Low-complexity points are converted into a sequence of amino acids (low complexity blobs, LCBs) and trimmed.
        4. Low-complexity regions (LCRs) are formed by merging intersecting LCBs.
    '''
    def __init__(self, uid: str, sequence: str) -> None:
        '''
        Initializes the Shannon entropy-based low complexity region detector.
        
        Args:
            uid: unique identifier for the protein sequence
            sequence: amino acid sequence as a string
            
        Attributes:
            window: sliding window size for entropy calculation (default: 3)
            rand_entropy: reference entropy value for random sequence of window size
            lcr: detected low complexity regions (initialized as None)
        '''
        self.uid = uid
        self.sequence = sequence
        self.window = 3
        self.rand_entropy = entropy([1/self.window]*self.window)
        self.lcr = None
    
    def determine_lcp(self, entropies: Iterator[float]) -> Iterator[int]:
        '''
        Identifies Low Complexity Points (LCPs) based on Shannon entropy threshold.
        
        A position is considered a low complexity point if its local entropy
        (calculated over the sliding window) is lower than the expected entropy
        of a random amino acid distribution.
        
        Args:
            entropies: iterator of Shannon entropy values for each window position
            
        Returns:
            Iterator of indices where entropy falls below the random sequence threshold
        '''
        return (i for i, aa_entropy in enumerate(entropies) if aa_entropy < self.rand_entropy)
                
    def determine_lcb(self, lcp: Iterator[int]) -> Iterator[Tuple[str, int, int]]:
        '''
        Converts contiguous Low Complexity Points into Low Complexity Blobs (LCBs).
        
        Groups consecutive LCPs and extends them by window size to form continuous
        amino acid sequences. Each LCB represents a candidate low-complexity region.
        
        Args:
            lcp: iterator of low complexity point indices
            
        Returns:
            Iterator of low complexity blobs as (sequence, start, end) tuples
        '''
        for _, group in groupby(enumerate(lcp), lambda x: x[0] - x[1]):
            lcp_group = [pos for _, pos in group]
            start, end = lcp_group[0], lcp_group[-1] + self.window - 1
            lcb = self.sequence[start: end + 1]
            yield self.trimming_lcb(lcb, start , end)
    
    def trimming_lcb(self, lcb: str, start: int, end: int) -> Tuple[str, int, int]:
        '''
        Trims Low Complexity Blobs by removing flanking unique residues.
        
        Args:
            lcb: low complexity blob sequence
            start: starting position in original sequence
            end: ending position in original sequence
            
        Returns:
            Tuple of (trimmed_sequence, new_start, new_end)
        '''
        if lcb.count(lcb[0]) == 1:
            lcb = lcb[1:]
            start += 1
        if lcb.count(lcb[-1]) == 1:
            lcb = lcb[:-1]
            end -= 1
        return lcb, start, end
    
    def determine_lcr(self, lcb: Iterator[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
        '''
        Merges overlapping low complexity blobs into final low complexity regions.
        
        Combines LCBs that are directly adjacent (end of one blob equals start of next)
        to form larger, continuous low-complexity regions, eliminating artificial
        fragmentation.
        
        Args:
            lcb: iterator of low complexity blobs as (sequence, start, end) tuples
            
        Returns:
            List of merged low complexity regions as (sequence, start, end) tuples
        '''
        lcr = list(lcb)
        i = 0
        while i < len(lcr) - 1:
            if lcr[i][2] - lcr[i + 1][1] == 0:
                lcr[i] = self.sequence[lcr[i][1]: lcr[i + 1][2] + 1], lcr[i][1], lcr[i + 1][2]
                del lcr[i + 1]
            else:
                i += 1
        return lcr

    def run(self) -> None:
        '''
        Executes the complete low complexity region detection pipeline:
            1. Calculate Shannon entropy for all sliding window positions
            2. Identify low complexity points based on entropy threshold
            3. Convert points to low complexity blobs with trimming
            4. Merge adjacent blobs into final low complexity regions
        
        Returns: 
            None

        Side effects:
            Updates self.lcr with detected low complexity regions
        """
        '''
        entropies = seq2smth.seq2entropy(self.sequence, self.window)
        lcp = self.determine_lcp(entropies)
        lcb = self.determine_lcb(lcp)
        self.lcr = self.determine_lcr(lcb)
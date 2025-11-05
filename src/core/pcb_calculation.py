import statistics
from collections.abc import Iterator
from itertools import groupby
from typing import List, Tuple

import numpy as np
from sklearn.metrics import root_mean_squared_error

from . import seq2smth

class BlockinessCalculator:
    '''
    Calculates charge blockiness patterns in protein sequences.

    Charge blockiness identifies contiguous regions with homogeneous charge distribution,
    which correlates with protein phase separation propensity.

    Algorithm:
    1. Cutoff frequency determination via binary search until RMSE between
    LPF-filtered full sequence and low-complexity regions converges below threshold
    2. Apply ideal low-pass filter to charge sequence using determined cutoff
    3. NCPR-based thresholding to classify residues as positive/negative/neutral
    4. Detect contiguous blocks and assign average block charge to each residue

    The approach smooths local charge fluctuations while preserving biologically
    relevant charge patterns.
    '''
    def __init__(self, uid: str, sequence: str, lcr: Iterator[Tuple[str, int, int]], rmse: float = 0.1) -> None:
        '''
        Initialize the BlockinessCalculator.

        Args:
            uid: unique identifier for the protein sequence
            sequence: amino acid sequence as a string
            lcr: iterator of low-complexity regions as (lcr_sequence, start_index, end_index) tuples
            rmse: root mean square error
        '''
        self.uid = uid
        self.sequence = sequence
        self.seq_charge = tuple(seq2smth.seq2charge(self.sequence))
        self.threshold_seq = None

        self.lcr = lcr
        self.lcr_charge = tuple(seq2smth.lcr2charge(self.sequence, self.lcr))
        self.rmse = rmse
        self.fft_size = len(sequence)
        self.cutoff_freq = None
        self.lpf_seq, self.lpf_lcr = None, None

        self.ncpr = None
        self.blockiness = None

    def ideal_lpf(self, seq_charge: Iterator[int], cutoff_freq: int) -> np.ndarray:
        '''
        Applies an ideal low-pass filter to the charge sequence.

        Args:
            seq_charge: input charge sequence as numeric values
            cutoff_freq: cutoff frequency for the filter
            
        Returns:
            Filtered signal in time domain as numpy array
        '''
        freqs = np.fft.fftfreq(self.fft_size, 1/self.fft_size)
        mask = np.abs(freqs) <= cutoff_freq
        fft = np.fft.fft(seq_charge)
        LPF_spectrum = fft * mask
        LPF_signal = np.real(np.fft.ifft(LPF_spectrum))
        return LPF_signal

    def determine_cutoff_freq(self) -> None:
        '''
        Determines optimal cutoff frequency using binary search based on RMSE convergence.
        Uses iterative frequency adjustment until RMSE between filtered sequence
        and low-complexity regions meets the threshold.
        
        Side effects:
            Updates cutoff_freq, lpf_seq, and lpf_lcr attributes
        '''
        cutoff_freq = self.fft_size
        delta = cutoff_freq // 2
        
        while True:
            lpf_seq = self.ideal_lpf(self.seq_charge, cutoff_freq)
            lpf_lcr = self.ideal_lpf(self.lcr_charge, cutoff_freq)

            lpf_seq_regs = (lpf_seq[start: end + 1] for _, start, end in self.lcr)
            lpf_lcr_regs = (lpf_lcr[start: end + 1] for _, start, end in self.lcr)

            rmse = (root_mean_squared_error(seq_reg, lcr_reg) for seq_reg, lcr_reg in zip(lpf_seq_regs, lpf_lcr_regs))
            avg_rmse = statistics.mean(rmse)

            if avg_rmse > self.rmse:
                cutoff_freq -= delta
            else:
                cutoff_freq += delta - delta//2
                delta //=2
            
            if delta == 0:
                self.cutoff_freq = cutoff_freq
                self.lpf_seq, self.lpf_lcr = lpf_seq, lpf_lcr
                break

    def ncpr_threshold_filter(self) -> np.ndarray:
        '''
        Applies NCPR-based thresholding to filtered sequence.
        Decisive rule:
            - 1 if charge > NCPR threshold
            - -1 if charge < -NCPR threshold  
            - 0 otherwise
        
        Returns:
            Integer array of classified charges
        '''
        self.ncpr = np.abs(sum(self.seq_charge)/len(self.seq_charge))
        return np.where(self.lpf_seq > self.ncpr, 1, np.where(self.lpf_seq < -self.ncpr, -1, 0))

    def determine_blockiness(self) -> List[float]:
        '''
        Groups consecutive residues with same threshold classification and
        assigns the average charge of each block to all residues in that block.

        Returns:
            List of charge values where each residue gets its block's average charge
        
        Special case: single-residue blocks are assigned the global NCPR value
        rather than their individual charge to avoid over-representing isolated residues.    
        '''
        blockiness = []
        for delta, group in groupby(enumerate(self.threshold_seq), lambda x: x[1]):
            block = [i for i, charge in group]
            start, end = block[0], block[-1] + 1

            avg_block_charge = self.ncpr if len(block) == 1 else sum(self.seq_charge[start: end])/len(block)
            block_charge = [avg_block_charge]*len(block)
            blockiness.extend(block_charge)
        return blockiness

    def run(self) -> None:
        '''
        Executes the complete blockiness calculation pipeline:
            1. Determine optimal cutoff frequency
            2. Apply NCPR threshold filtering  
            3. Calculate block-wise charge densities
        '''
        self.determine_cutoff_freq()
        self.threshold_seq = self.ncpr_threshold_filter()
        self.blockiness = self.determine_blockiness()
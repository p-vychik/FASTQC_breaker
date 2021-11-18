
'''
Fastq file example:

@EAS54_6_R1_2_1_413_324
CCCTTCTTGTCTTCAGCGTTTCTCC
+
;;3;;;;;;;;;;;;7;;;;;;;88
@EAS54_6_R1_2_1_540_792
TTGGCAGGCCAAGGCCGATGGATCA
+
;;;;;;;;;;;7;;;;;-;;;3;83
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+EAS54_6_R1_2_1_443_348
;;;;;;;;;;;9;7;;.7;393333
'''

import numpy as np
import sys


# -----------------------------------------------------------------------------
class Sequences:
    '''
        Representation of set of multiple reads as numpy 2D array:
        It allows fast processing big data.
        Each matrix line contents array of sequence qualities.
    '''

    def __init__(self, seq_len, read_num):
        # print(f' read_num={read_num}, seq_len={seq_len}')
        self.count = 0
        self.max_q = 0
        self.max_meanq = 0
        self.min_meanq = 100
        self.max_len = seq_len
        self.seq_mat = np.empty((read_num, seq_len), dtype=np.int8)
        self.qual_mat = np.empty((read_num, seq_len), dtype=np.int8)
        # self.tile_arr = np.empty((read_num), dtype=np.int32)
        self.meanq_arr = np.zeros((94), dtype=np.uint32)

    def print(self):
        # print(self.tile_arr, self.tile_arr.shape)
        print(self.seq_mat, self.seq_mat.shape)
        print(self.qual_mat, self.qual_mat.shape)
        print(self.meanq_arr, self.meanq_arr.shape)
        print(f'sz={sys.getsizeof(self.seq_mat)} bytes in memory')

    def add(self, descr_line, seq_line, qual_line):
        lst = [ord(ch) for ch in seq_line]
        idx = self.count
        sl = len(seq_line)
        # make array if tiles
        if 0:
            tile = int(descr_line.split(sep=':')[4])
            self.tile_arr[idx] = tile
        # make matrix of ascii
        if 0:  # I don't need seq yet
            lst = [ord(ch) for ch in seq_line]
            for i in range(sl, self.max_len):
                lst.append(-1)
            self.seq_mat[idx] = lst
        # make matrix of qual
        if 1:
            lst = [ord(ch)-ord('!') for ch in qual_line]
            mean_q = int(sum(lst) / len(lst))
            for i in range(sl, self.max_len):
                lst.append(-1)
            self.qual_mat[idx] = lst
        # increment value in mean-quility array
        self.meanq_arr[mean_q] += 1
        self.max_meanq = max(mean_q, self.max_meanq)
        self.min_meanq = min(mean_q, self.min_meanq)
        # save max quality
        mx = max(lst)
        self.max_q = max(self.max_q, mx)
        self.count += 1

    def get_sequences(self):
        return self.sequences

    def max_length(self):
        return self.max_len

    def max_qual(self):
        return self.max_q

    def get_qual_column(self, n):
        return self.qual_mat[:, n]

    def get_transpose_qual_matrix(self):
        return self.qual_mat.transpose()

    def get_meanq_arr(self):
        return self.meanq_arr, self.min_meanq, self.max_meanq

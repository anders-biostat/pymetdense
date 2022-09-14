import numpy as np
import os
import bisect
import warnings
import pandas as pd
from typing import Sequence, Union, Optional
    

class MetDenseFile:
    """ Class to load data from a .metdense file.
        The upper and lower bounds of position ranges are always treated as inclusive.
        
        filename: path to the file
        verbose: verbosity of the initialization, in [0,3]
        
        Loading data:
            read_data(chromosome: Union[int, str], lower: int=0, upper: int=1e30)
            read_positions(chromosome: Union[int, str], lower: int=0, upper: int=1e30)
            read_data_positions(chromosome: Union[int, str], ranges: Sequence[Sequence[int]], mode='block')
        
        Relevant helpers:
            get_chromosome_length(chromosome: Union[int, str])
            set_cell_mask(mask: Sequence[Union[int, bool]])
            get_cell_names()
        
        
        TODO:
            - Add split loading and masking if requested row number is very high!
        
    """
    def __init__(self, filename: str, verbose: int=0):
        self.filename = filename # Filename
        self.size = os.path.getsize(filename) # File size in bytes
        self.size_in_blocks = self.size//4 # in 4-blocks
        self.file = open(self.filename, "rb")
        self.verbose = verbose
        if self.verbose>0: print("Load file",self.filename)
        self.endianess = "little" # Should also work with "big", but untested!
        self.nptype = "<u4" if self.endianess == "little" else ">u4"
        self.cells_masked = False
        self.raise_position_errors = False # Whether errors should be thrown when a requested region doesn't contain any reads
        
        self.__read_intro()
        self.__read_chromosomes_map()
    
    ###### Basic Helper Functions
    def __len__(self):
        """ Length, thinks in 4-blocks.
        """
        return self.size_in_blocks
    
    def __getitem__(self, i: int):
        """ To run binary search directly on this object, thinks in 4-blocks.
            
            i: address of the value
        """
        if not i < self.size_in_blocks: raise IndexError("Index is beyond end of file.")
        return self.__get_int(4*i+self.nonconform_offset)
    
    def __get_int(self, position: Optional[int]=None, size: Optional[int]=4):
        """ Read a single number from file.
            
            position: defaults to current, otherwise byte position in file
            size: number of bytes to read into int, by default 4 (uint32)
        """ # 
        if position==None:
            return int.from_bytes(self.file.read(size),self.endianess)
        else:
            self.file.seek(position)
            return int.from_bytes(self.file.read(size),self.endianess)
    
    def __get_lines(self, n: int):
        """ Read n lines from current file position.
        """
        lines = list()
        for i in range(n):
            lines.append(self.file.readline()[:-1].decode())
        return lines
    
    ###### Read initial stuff
    def __read_intro(self):
        """ Read the intro part of the metdense file.
        """
        self.file.seek(8)
        self.vMajor, self.vMinor = [self.__get_int() for i in range(2)] # Version number
        self.positions_bytesize = 4 if (self.vMajor==0 and self.vMinor==0) else 8 # uint64 for positions from 0.1 onward
        self.data_start, self.chromosome_start = [self.__get_int(size=self.positions_bytesize) for i in range(2)]
        # Check if the metdense file conforms to the specification, not starting the data block at a multiple of four is easy to screw up. Should still work.
        if self.data_start % 4 != 0: warnings.warn("The file does not conform to the specification! Data block starts at "+str(self.data_start)+". Parser should still work.")
        self.nonconform_offset = self.data_start % 4
        self.n_cells = self.__get_int()
        self.cell_names = self.__get_lines(self.n_cells)
        self.data_row_size = int(np.ceil(self.n_cells/16)) # Row size in uint32
        self.data_row_size_bytes = 4*self.data_row_size # Row size in bytes
        if self.verbose>0: print("The file follows version {}.{} and contains data for {} cells.".format(self.vMajor,self.vMinor,self.n_cells))
     
    def __read_chromosomes_map(self):
        """ Read the chromosome map.
        """
        self.file.seek(self.chromosome_start)
        self.n_chromosomes = self.__get_int()
        self.positions_start = [self.__get_int(size=self.positions_bytesize) for i in range(self.n_chromosomes)]
        self.positions_end = self.positions_start[1:]+[self.chromosome_start]
        self.data_block_size = self.positions_start[0]-self.data_start
        self.chromosome_names = self.__get_lines(self.n_chromosomes)
        self.n_positions = (self.positions_end[-1]-self.positions_start[0])//4
        self.n_positions_per_chromosome = [(self.positions_end[i] - self.positions_start[i])//4 for i in range(self.n_chromosomes)]
        if self.data_block_size/self.data_row_size_bytes != self.n_positions: warnings.warn("The file does not conform to the specification! Data block size, number of cells and number of CpGs does not match.")
        if self.verbose>1: print("Contains Chromosomes:",self.chromosome_names)
        if self.verbose>2: print("Chromosomes have lengths:",self.n_positions_per_chromosome,"with total length",self.n_positions)
    
    ###### Data access
    #### Helpers
    def get_chromosome_length(self, chromosome: Union[int, str]):
        """ Returns the number of available positions on a chromosome.
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        pos = self.chromosome_names.index(str(chromosome))
        return self.n_positions_per_chromosome[pos]
    
    def set_cell_mask(self, mask: Sequence[Union[int, bool]]):
        """ Set cell mask to be applied to all returned data.
            Can either be a boolean array with same dimension as self.n_cells,
            or a list of indices in arbitrary order.
        """
        self.cells_masked = True
        self.cell_mask = mask
    
    def get_cell_names(self):
        """ Get cell names.
        """
        return self.cell_names
    
    def __get_position_range(self, chromosome: Union[int, str], lower: Optional[int]=0, upper: Optional[int]=1e30):
        """ Returns the position range in the positions block for a chromosome CpG range.
            
            chromosome: name of the chromosome
            lower: lower range on the chromosome (inclusive)
            upper: upper range on the chromosome (inclusive)
            
            lower_index: index of the first position in the range
            upper_index: index of the first position above the range (or on the next chromosome)
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        if self.raise_position_errors and not lower <= upper: raise ValueError("Upper range is not larger than lower range.")
        pos = self.chromosome_names.index(str(chromosome))
        lower_bound, upper_bound = self.positions_start[pos]//4, self.positions_end[pos]//4
        lower_bound_pos, upper_bound_pos = self[lower_bound], self[upper_bound-1]
        if self.raise_position_errors and not lower_bound_pos <= upper: raise ValueError("Upper bound does not include anything. Chromosome "+str(chromosome)+", between "+str(lower)+" and "+str(upper)+".")
        if self.raise_position_errors and not upper_bound_pos >= lower: raise ValueError("Lower bound does not include anything. Chromosome "+str(chromosome)+", between "+str(lower)+" and "+str(upper)+".")
        lower_index = 4*bisect.bisect_left(self, lower, lo=lower_bound, hi=upper_bound) if lower > lower_bound_pos else 4*lower_bound
        if lower == upper:
            if not self.__get_int(lower_index) == lower: raise ValueError("A single position was requested, but does not exist.")
            return lower_index, lower_index+4
        upper_index = 4*bisect.bisect_right(self, upper, lo=lower_bound, hi=upper_bound) if upper < upper_bound_pos else 4*upper_bound
        return lower_index, upper_index # upper includes the index of the first CpG that is no longer included!
    
    def __position_to_row(self, position: int):
        """ Turns a position in the position block into the row number in the data block.
        """
        row = (position - self.positions_start[0])//4
        return row
    
    def __position_to_row_on_chromosome(self, chromosome: Union[int, str], position: int):
        """ Turns a position in the position block into the row number in the data block of the chromosome.
        """
        if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
        pos = self.chromosome_names.index(str(chromosome))
        row = (position - self.positions_start[pos])//4
        return row
    
    def __position_to_row_start(self, position: int):
        """ Turns a position in the position block into the position in the file where the corresponding row starts.
        """
        row = self.__position_to_row(position)
        row_position = self.data_start + row*self.data_row_size_bytes
        return row_position
    
    def __row_to_row_start(self, row: int, chromosome: Union[int, str]):
        """ Turns a row in the data block into the position in the file where the row starts.
        """
        start = self.data_start + row*self.data_row_size_bytes
        if chromosome:
            if not str(chromosome) in self.chromosome_names: raise ValueError("Requested chromosome "+str(chromosome)+" does not exist.")
            pos = self.chromosome_names.index(str(chromosome))
            start += (self.positions_start[pos] - self.positions_start[0])//4*self.data_row_size_bytes
        return start
    
    def __transform_position_ranges(self, chromosome: Union[int, str], ranges: Sequence[Sequence[int]]):
        """ Transform a list of ranges on the chromosome, into a list of the corresponding row ranges in the data block of that chromosome.
        """
        CpG_ranges = np.asarray([self.__get_position_range(chromosome, lower=min(ran[0],ran[1]), upper=max(ran[0],ran[1])) for ran in ranges], dtype=object)
        return self.__position_to_row_on_chromosome(chromosome, CpG_ranges), CpG_ranges.min(), CpG_ranges.max()
        
    def __transform_position_ranges_to_indices(self, chromosome: Union[int, str], ranges: Sequence[Sequence[int]]):
        """ Transform a list of ranges on the chromosome, into ...
        """
        row_ranges, pos_min, pos_max = self.__transform_position_ranges(chromosome, ranges)
        rmin, rmax = row_ranges.min(), row_ranges.max()
        return np.asarray([np.linspace(l, u-1, u-l, dtype=int) for l,u in row_ranges],dtype=object), rmin, rmax, pos_min, pos_max
    
    def __build_load_mask(self, ranges: Sequence[Sequence[int]], length: int):
        """ Build mask from ranges.
        """
        mask = np.full((length),False)
        for ran in ranges:
            mask[ran] = True
        return mask
    
    def __mask_ranges_transform(self, mask: Sequence[bool], ranges: Sequence[Sequence[int]]):
        """ Transform indices in ranges, taking into account that some parts in between will get deleted.
        """
        key = (np.linspace(0,mask.shape[0]-1,mask.shape[0])-np.invert(mask).cumsum()).astype(int)
        ranges_tf = list()
        for ran in ranges:
            ranges_tf.append(key[ran])
        return ranges_tf
    
    #### Read Data
    def __read_data_rows(self, lower: int, upper: int, mask = None):
        """ Read a range of data rows from the data block.
            Loading is almost instantaneous (for small ranges), most of the time is spent on the shift and the concatenation in __process_data_rows.
            (equal parts, no idea why concatenation is so slow; allocating the memory beforehand and setting parts of the new array is even slower)
            
            lower: lower row number (inclusive)
            upper: upper row number (exclusive)
        """
        self.file.seek(0)
        data = np.reshape(np.fromfile(self.file,self.nptype,offset=lower+self.nonconform_offset ,count=(upper-lower)//4),(-1,self.data_row_size))
        if not mask is None:
            data = data[mask].copy()
        return self.__process_data_rows(data)
    
    def __process_data_rows(self, data):
        """ Process the output of __read_data_rows.
            Takes the binary data, turns it into a numpy array of shape (positions, cells).
            If a cell_mask is set, this is also applied directly to the cell axis here.
        """
        result = list()
        for i in range(16):
            result.append(np.bitwise_and(data,0x03).astype(np.uint8))
            data = np.right_shift(data,2)
        del data # No reason to keep it in the memory
        result = [np.array_split(res,res.shape[1],axis=1) for res in result] # Reuse the variable, or just delete the old one afterwards
        result = list(map(list, zip(*result)))
        result = np.concatenate([el for ls in result for el in ls],axis=1)
        result = result[:,:self.n_cells]
        if self.cells_masked:
            result = result[:,self.cell_mask]
        return result
    
    def read_data(self, chromosome: Union[int, str], lower: int=0, upper: int=1e30, mask: Sequence[Union[int, bool]]=None):
        """ Read all data within a single position range on a chromosome.
            Can apply a mask along the position axis to the data before processing;
            this is faster than applying it afterwards.
            
            chromosome: which chromosome
            lower: lower position on the chromosome (inclusive)
            upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.__get_position_range(chromosome, lower, upper) # Position to index
        row_lower, row_upper = self.__position_to_row_start(lower_index), self.__position_to_row_start(upper_index) # Index to row in data block
        return self.__read_data_rows(row_lower, row_upper, mask) # Load, process rows
    
    #### Read Positions
    def read_positions(self, chromosome: Union[int, str], lower: int=0, upper: int=1e30):
        """ Read all available positions within a range on a chromosome.
            
            - chromosome: which chromosome
            - lower: lower position on the chromosome (inclusive)
            - upper: upper position on the chromosome (inclusive)
        """
        lower_index, upper_index = self.__get_position_range(chromosome, lower, upper)
        self.file.seek(0)
        positions = np.fromfile(self.file,self.nptype,offset=lower_index+self.nonconform_offset ,count=(upper_index-lower_index)//4)
        return positions
    
    #### Read Data & Positions
    def __read_data_positions_from_ranges_asblock(self, chromosome: Union[int, str], ranges: Sequence[Sequence[int]]):
        """ Read all data for multiple position ranges on a chromosome.
            Loads all data between the lowest and highest position in ranges,
            and selects the relevant parts before processing further.
            Fastest option if the ranges are close or overlapping, but
            potentially very large RAM requirements.
            
            chromosome: which chromosome
            ranges
        """
        range_indices, rmin, rmax, pos_min, pos_max = self.__transform_position_ranges_to_indices(chromosome, ranges)
        range_indices -= rmin
        mask = self.__build_load_mask(range_indices, rmax-rmin)
        
        data = self.__read_data_rows(self.__row_to_row_start(rmin, chromosome), self.__row_to_row_start(rmax, chromosome), mask)
        range_indices = self.__mask_ranges_transform(mask, range_indices)
        final = [data[ran] for ran in range_indices]
        
        self.file.seek(0)
        positions = np.fromfile(self.file,self.nptype,offset=pos_min+self.nonconform_offset ,count=(pos_max-pos_min)//4)[mask]
        pos = [positions[ran] for ran in range_indices]
        
        return final, pos
    
    def __read_data_positions_from_ranges_separately(self, chromosome: Union[int, str], ranges: Sequence[Sequence[int]]):
        """ Read all data for multiple position ranges on a chromosome.
            Load data separately for every range, fastest option if the
            ranges are very far apart from each other.
            
            chromosome: which chromosome
            ranges
        """
        data = []
        pos = []
        for ran in ranges:
            pos.append(self.read_positions(chromosome, ran[0], ran[1]))
            data.append(self.read_data(chromosome, ran[0], ran[1]))
        
        return data, pos
    
    def read_data_positions(self, chromosome: Union[int, str], ranges: Sequence[Sequence[int]], mode='block'):
        """ Read all data for multiple position ranges on a chromosome.
            
            chromosome: which chromosome
            ranges
            mode: 'block'       Loads all data between the lowest and highest position in ranges,
                                and selects the relevant parts before processing further.
                                Fastest option if the ranges are close or overlapping, but
                                potentially very large RAM requirements.
            mode: 'separate'    Load data separately for every range, fastest option if the
                                ranges are very far apart from each other.
        """
        if not mode in ['block', 'separate']: raise ValueError("Requested data load mode "+mode+" does not exist.")
        if mode == 'block':
            return self.__read_data_positions_from_ranges_asblock(chromosome, ranges)
        else:
            return self.__read_data_positions_from_ranges_separately(chromosome, ranges)
    






import pickle
import numpy as np
import pyBigWig as pybw
from utilities import chrom_sizes

"""
Constructor: creates an instance of a genome database with modular tracks
at a user-specified resolution

Inputs:
resolution - the resolution of the genome
"""

class Track:

    def __init__(self, name, resolution: int, sizes_file):
        self.chrom_data = {}
        self.name = name
        self.resolution = int(resolution)
        self.sizes = chrom_sizes(sizes_file)

        for chrom in self.sizes:
            self.chrom_data[chrom] = np.zeros((self.sizes[chrom] // self.resolution + 1,)).astype(object)

    def get(self, chrom, start, end=None):
        start = int(start // self.resolution)
    
        if end is None or end == start:
            end = start // self.resolution + 1
        else:
            end = int(end // self.resolution + 1)

        return self.chrom_data[chrom][start:end]

    def from_bed(self,bed_file):
        f = open(bed_file,'r')

        for line in f:
            data = line.split()
            chrom = data[0]
            start = int(data[1]) // self.resolution
            end = int(data[2]) // self.resolution
            try:
                value = float(data[3])
            except:
                value = data[3]
            

            self.chrom_data[chrom][start:end] = value

        return self

    def from_bigWig(self,bigWig_file):
        bw = pybw.open(bigWig_file)
        for chrom in self.chrom_data:
            for i in range(len(self.chrom_data[chrom])):
                start = i * self.resolution
                end = np.min([(i + 1) * self.resolution,self.sizes[chrom]])
                
                try:
                    value = bw.stats(chrom,start,end)[0]
                except RuntimeError:
                    value = 0

                if value is None:
                    value = 0

                self.chrom_data[chrom][i] = value
        
        return self

    """
    Input: gtf file in Ensembl V95 format or hg19 knownGene format on UCSC
    genome browser
    """
    def from_gtf(self,gtf_file):
        pass

class Sequential:
    def __init__(self,name,tracks):
        self.name = name
        self.tracks = {}

        for track in tracks:
            self.tracks[track.name] = track

    def get(self, chrom, start, end=None):
        data = {}

        for track_name in self.tracks:
            data[track_name] = self.tracks[track_name].get(chrom,start,end)

        return data

class Genome:
    def __init__(self, resolution: int):
        self.resolution = int(resolution)
        self.tracks = {}

    def add_track(self, track):
        track_name = track.name

        self.tracks[track_name] = track

    def remove_track(self, track):
        if type(track) == Track:
            track_name = track.name
        elif type(track) == str:
            track_name = track
        else:
            return TypeError("Invalid track type - must be either Track or string")

        if track_name in self.tracks:
            del self.tracks[track_name]
        else:
            return NameError("Track name %s not found in dataset" % track.name)


    """
    Queries data in all tracks in chromosome <chrom> within interval [<start>,<end>]

    <chrom> must be formatted as 'chrN' where N is the chromosome number/symbol
    <start> will be rounded down to the nearest multiple of the genome resolution
    <end> will be rounded up
    if <start> and <end> are equal or end is not specified, <end> is set to
    <start> + self.resolution

    Outputs a dict where each key is the track name and values are the data
    returned by each track.
    """
    def get(self, chrom, start, end=None):
        
        data = {}

        for track_name in self.tracks:
            data[track_name] = self.tracks[track_name].get(chrom,start,end)

        return data

    """
    Queries data in one specific track in chromosome <chrom> within interval
    [<start>,<end>]

    Outputs data returned by one track.
    """
    def get_track(self, track_name, chrom, start, end):

        return self.tracks[track_name].get(chrom,start,end)

    """Saves genome into a pickle-able file format"""
    def save(self, fname):
        pickle.dump(self, open(fname,'wb'))
import pickle
import h5py
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

    def __init__(self, name, resolution: int, sizes):
        self.chrom_data = {}
        self.name = name
        self.resolution = int(resolution)
        self.sizes = sizes

        for chrom in self.sizes:
            self.chrom_data[chrom] = np.zeros((self.sizes[chrom] // self.resolution + 1,))

    def get(self, chrom, start, end=None):
        start = int(start // self.resolution)
    
        if end is None or end == start:
            end = start // self.resolution + 1
        else:
            end = int(end // self.resolution + 1)

        return self.chrom_data[chrom][start:end]

    def from_bed(self,bed_file):
        f = open(bed_file,'r')
        typecast = False

        for line in f:
            data = line.split()
            chrom = data[0]

            try:
                start = int(data[1]) // self.resolution
                end = int(data[2]) // self.resolution
            except ValueError:
                continue

            try:
                value = float(data[3])
            except:
                value = data[3]
                if not typecast:
                    for c in self.sizes:
                        self.chrom_data[c] = self.chrom_data[c].astype(str)
                    typecast = True

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
    def __init__(self,name,tracks=[]):
        self.name = name
        self.tracks = {}

        for track in tracks:
            self.tracks[track.name] = track

    def __getitem__(self,track_name):
        return self.tracks[track_name]

    def add_track(self, track):
        self.tracks[track.name] = track

    def get(self, chrom, start, end=None):
        data = {}

        for track_name in self.tracks:
            data[track_name] = self.tracks[track_name].get(chrom,start,end)

        return data

class Genome:
    def __init__(self, resolution, sizes_file):
        self.resolution = int(resolution)
        self.tracks = {}
        self.sizes = chrom_sizes(sizes_file)

    def __getitem__(self, track_name):
        return self.get_track(track_name)

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
    <end> will be rounded up if <start> and <end> are equal or end is not specified,
        <end> is set to <start> + self.resolution

    Outputs a dict where each key is the track name and values are the data
    returned by each track.
    """
    def get(self, chrom, start, end=None, tracks=None, to_numpy=False):
        
        data = {}

        if tracks is None:
            tracks = [track for track in self.tracks]

        for track_path in tracks:
            track_path_array = track_path.split('/')

            track = self.tracks[track_path_array[0]]
            
            for i in range(1,len(track_path_array)):
                track = track[track_path_array[i]]

            data[track_path] = track.get(chrom,start,end)

        if to_numpy:
            return np.array([data[track] for track in data])

        return data

    """
    Queries data in one specific track in chromosome <chrom> within interval
    [<start>,<end>]

    Outputs data returned by one track.
    """
    def get_track(self, track_name):

        return self.tracks[track_name]

    """Saves genome into a pickle file format"""
    def save(self, fname):

        pickle.dump(self, open(fname,'wb'))
        
        """f = h5py.File(fname,'w')
        f['resolution'] = self.resolution
        recursively_save_dict_contents_to_group(f,'',self.tracks)

        f.close()"""

    """Saves genome from a hdf5 file format"""
    def load(self, fname):

        G = pickle.load(open(fname,'rb'))
        self.tracks = G.tracks

        """f = h5py.File(fname)
        self.resolution = f['resolution'][()]
        self.tracks = recursively_load_dict_contents_from_group(f,'/',self)"""

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, Sequential):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item.tracks)
        elif isinstance(item, Track):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item.chrom_data)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def recursively_load_dict_contents_from_group(h5file, path, genome):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            print(item)
            track_name = path.split('/')[-2]
            track = Track(track_name, genome.resolution, genome.sizes)
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            if 'value' not in dir(item):
                ans[key] = Sequential(key)
            
            ans[key].tracks = recursively_load_dict_contents_from_group(h5file, path + key + '/', genome)

    return ans
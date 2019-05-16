from Genome import Genome, Track, Sequential
import pprint

g = Genome(1e5, '../hg19.chrom.sizes')
tracks = (
    Sequential('Replication_Timing', tracks=(
        Track('G1',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_G1.bed'),
        Track('S1',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S1.bed'),
        #Track('S2',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S2.bed'),
        #Track('S3',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S3.bed'),
        #Track('S4',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S4.bed'),
        #Track('G2',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_G2.bed'),
    )),
    Track('Subcompartments',1e5,'../hg19.chrom.sizes').from_bed('data/GM12878_subcompartments.bed'),
)

""""""

for track in tracks:
    g.add_track(track)

g.save('genomes/GM12878.mg')

g = Genome(1e5,'../hg19.chrom.sizes')
g.load('genomes/GM12878.mg')
pprint.pprint(g.get('chr2',1e5,5e5,tracks=['Replication_Timing/G1','Replication_Timing/S1','Subcompartments'])['Subcompartments'])
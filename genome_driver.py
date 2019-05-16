from Genome import Genome, Track, Sequential
import pprint

g = Genome(1e5)
tracks = (
    Sequential('Replication_Timing', tracks=(
        Track('G1',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_G1.bed'),
        Track('S1',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S1.bed'),
        Track('S2',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S2.bed'),
        Track('S3',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S3.bed'),
        Track('S4',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_S4.bed'),
        Track('G2',1e5,'../hg19.chrom.sizes').from_bed('/mnt/d/data/sniper/repli_seq/GM12878/beds/GM12878_G2.bed'),
    )),
    Sequential('Histone_Marks', tracks=(
        Track('H3K4me1',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K4me1.bigWig'
        ),
        Track('H3K4me2',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K4me2.bigWig'
        ),
        Track('H3K4me3',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K4me3.bigWig'
        ),
        Track('H3K9ac',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K9ac.bigWig'
        ),
        Track('H3K9me3',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K9me3.bigWig'
        ),
        Track('H3K27ac',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K27ac.bigWig'
        ),
        Track('H3K36me3',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K36me3.bigWig'
        ),
        Track('H3K79me2',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H3K79me2.bigWig'
        ),
        Track('H4K20me1',1e5,'../hg19.chrom.sizes').from_bigWig(
            '/mnt/d/data/sniper/histone_marks/bigWigs_GM12878/GM12878_H4K20me1.bigWig'
        ),
    )),
    Track('Subcompartments',1e5,'../hg19.chrom.sizes').from_bed('data/GM12878_subcompartments.bed'),
)

""""""

for track in tracks:
    g.add_track(track)

g.save('genomes/genome_GM12878.mg')
pprint.pprint(g.get('chr2',1e5,5e5,to_numpy=True,tracks=['Replication_Timing/G1','Replication_Timing/S2','Replication_Timing/G2']))
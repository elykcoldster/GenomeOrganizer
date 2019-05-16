from Genome import Genome, Track, Sequential
import pprint

g = Genome(1e5, '../hg19.chrom.sizes')
tracks = (
    Sequential('histone_marks', tracks=(
        Sequential('cohesin_depleted', tracks=(
            Track('H3K27ac',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP010-H3K27Ac-treated.bw'
            ),
            Track('H3K4me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP012-H3K4me3-treated.bw'
            ),
            Track('H3K4me1',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP014-H3K4me1-treated.bw'
            ),
            Track('H3K36me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP016-H3K36me3-treated.bw'
            ),
            Track('H3K27me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP018-H3K27me3-treated.bw'
            ),
            Track('H3K9me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP020-H3K9me3-treated.bw'
            ),
            Track('H3K79me2',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP024-H3K79me2-treated.bw'
            ),
            Track('H4K20me1',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP026-H4K20me3-treated.bw'
            ),
        )),
        Sequential('wild_type', tracks=(
            Track('H3K27ac',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP009-H3K27Ac-untreated.bw'
            ),
            Track('H3K4me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP011-H3K4me3-untreated.bw'
            ),
            Track('H3K4me1',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP013-H3K4me1-untreated.bw'
            ),
            Track('H3K36me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP015-H3K36me3-untreated.bw'
            ),
            Track('H3K27me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP017-H3K27me3-untreated.bw'
            ),
            Track('H3K9me3',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP019-H3K9me3-untreated.bw'
            ),
            Track('H3K79me2',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP023-H3K79me2-untreated.bw'
            ),
            Track('H4K20me1',1e5,g.sizes).from_bigWig(
                '/mnt/d/data/sniper/histone_marks/bigWigs_HCT116/Rao-2017-CHIP025-H4K20me3-untreated.bw'
            ),
        ))
    )),
    Sequential('subcompartments', tracks= (
        Track('cohesin_depleted',1e5,g.sizes).from_bed('data/cd_cd_HCT116.bed'),
        Track('wild_type',1e5,g.sizes).from_bed('data/cd_wt_HCT116.bed'),
    ))
)

""""""

for track in tracks:
    g.add_track(track)

g.save('genomes/HCT116.mg')

pprint.pprint(g.get('chr2',1e5,5e5,to_numpy=True,tracks=[
    'histone_marks/cohesin_depleted/H3K27ac',
    'histone_marks/wild_type/H3K27ac',
    'histone_marks/cohesin_depleted/H3K36me3',
    'histone_marks/wild_type/H3K36me3',
]))
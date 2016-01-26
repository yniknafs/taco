'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from taco.lib.transfrag import Transfrag
from taco.lib.base import Exon, Strand


def test_introns():
    t = Transfrag(chrom='chrTest', strand=Strand.POS,
                  exons=[Exon(0, 10), Exon(20, 30), Exon(40, 50)])
    introns = list(t.iterintrons())
    assert len(introns) == 2
    assert introns[0] == (10, 20)
    assert introns[1] == (30, 40)


def test_splices():
    t = Transfrag(chrom='chrTest', strand=Strand.POS,
                  exons=[Exon(0, 10), Exon(20, 30), Exon(40, 50)])
    splice_sites = frozenset(t.itersplices())
    assert len(splice_sites) == 4
    assert splice_sites == frozenset((10, 20, 30, 40))

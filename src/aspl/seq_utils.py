import pysam
from .events import SpliceJunction, SpliceSiteType, Strand, SpliceSite

COMPLEMENT = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}


def get_site_sequence(ss: SpliceSite, fa: pysam.FastaFile):
    to_complement = ss.strand == Strand.MINUS

    if (ss.type == SpliceSiteType.DONOR) != (ss.strand == Strand.MINUS):
        seq = fa.fetch(ss.seqname, ss.coord, ss.coord + 2).upper()
    else:
        seq = fa.fetch(ss.seqname, ss.coord - 3, ss.coord - 1).upper()

    if to_complement:
        seq = "".join(COMPLEMENT[n] for n in seq)[::-1]

    return seq


def get_junction_sequence(ss: SpliceJunction, fa: pysam.FastaFile):
    return get_site_sequence(ss.donor_site, fa) + get_site_sequence(
        ss.acceptor_site, fa
    )

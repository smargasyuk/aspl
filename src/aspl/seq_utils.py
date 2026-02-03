import pysam
from .events import (
    SpliceJunction,
    SpliceSiteType,
    Strand,
    SpliceSite,
    IntervalCoordinate,
)

COMPLEMENT = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}


def reverse_complement(seq: str):
    return "".join(COMPLEMENT[n] for n in seq)[::-1]


def get_site_sequence(ss: SpliceSite, fa: pysam.FastaFile):
    to_complement = ss.strand == Strand.MINUS

    if (ss.type == SpliceSiteType.DONOR) != (ss.strand == Strand.MINUS):
        seq = fa.fetch(ss.seqname, ss.coord, ss.coord + 2).upper()
    else:
        seq = fa.fetch(ss.seqname, ss.coord - 3, ss.coord - 1).upper()

    if to_complement:
        seq = reverse_complement(seq)

    return seq


def get_junction_sequence(ss: SpliceJunction, fa: pysam.FastaFile):
    return get_site_sequence(ss.donor_site, fa) + get_site_sequence(
        ss.acceptor_site, fa
    )


def get_interval_sequence(ic: IntervalCoordinate, fa: pysam.FastaFile):
    seq = fa.fetch(ic.seqname, ic.seqStart, ic.seqEnd).upper()
    if ic.strand == Strand.MINUS:
        seq = reverse_complement(seq)
    return seq


def get_intron_sequence(sj: SpliceJunction, fa: pysam.FastaFile):
    return get_interval_sequence(sj.get_interval(), fa)

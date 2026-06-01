import pysam
from Bio.Seq import translate
from .events import (
    SpliceJunction,
    SpliceSiteType,
    Strand,
    SpliceSite,
    IntervalCoordinate,
    Exon,
)
from .transcripts import Transcript

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


def get_sequence(el: SpliceJunction | Exon, fa: pysam.FastaFile):
    return get_interval_sequence(el.get_interval(), fa)


def _upper_w_star(char: str):
    if char == "*":
        return "#"
    return char.upper()


def _lower_w_star(char: str):
    if char == "#":
        return "*"
    return char.lower()


def _translate_exon_with_splice_marks(seq: str, leftover: str) -> tuple[str, str]:
    if leftover:
        seq = leftover + seq
    last_coding = len(seq) // 3 * 3
    p1 = translate(seq[:last_coding]).lower()
    new_leftover = seq[last_coding:]
    if not p1:
        return p1, new_leftover
    p1 = _upper_w_star(p1[0]) + p1[1:]
    if last_coding == len(seq):
        return p1[:-1] + _upper_w_star(p1[-1:]), ""
    return p1, new_leftover


def _lower_first_and_last(s):
    if not s:
        return s
    if len(s) == 1:
        return _lower_w_star(s)
    return _lower_w_star(s[0]) + s[1:-1] + _lower_w_star(s[-1])


def translate_with_splice_marks(t: Transcript, fa: pysam.FastaFile) -> str:
    """
    Translates transcript to protein, converting amino acids
    at splice junctions to uppercase.
    """
    translated = ""
    nt_leftover = ""
    for e in t.iter_exons():
        exon_seq = get_sequence(e, fa)
        exon_translated, nt_leftover = _translate_exon_with_splice_marks(
            exon_seq, nt_leftover
        )
        translated += exon_translated
    # the first and the last exon boundaries are not splice sites
    translated = _lower_first_and_last(translated)
    return translated


def _last_uppercase(s):
    for i in range(len(s) - 1, -1, -1):
        if s[i].isupper() | (s[i] == "#"):
            return i
    return None


def is_poison_by_50nt(translation_with_splice_marks: str):
    end_position = translation_with_splice_marks.find("*")
    if end_position == -1:
        return False
    last_splice_position = _last_uppercase(translation_with_splice_marks)
    if last_splice_position is None:
        return False
    return last_splice_position - end_position - 1 >= 17

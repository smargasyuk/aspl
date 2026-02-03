from typing import Callable
from polars.datatypes.classes import R
from .events import IntervalCoordinate, Strand
import numpy as np

RelativeCoordinates = tuple[int, int]


# Get coordinates of an interval relative to another interval.
# This is for sequence retrieval and manipulation, and the coordinates are counted from the 5' end, not from the smaller coordinate
def get_relative_coordinates(
    inner: IntervalCoordinate, outer: IntervalCoordinate
) -> RelativeCoordinates:
    assert inner.seqname == outer.seqname
    assert inner.strand == outer.strand

    strand = inner.strand

    length = inner.seqEnd - inner.seqStart
    if strand == Strand.MINUS:
        start = outer.seqEnd - inner.seqEnd
    else:
        start = inner.seqStart - outer.seqStart

    return start, start + length


def capitalize_by_relative_coordinates(seq: str, rel_coords: RelativeCoordinates):
    return (
        seq[: rel_coords[0]].lower()
        + seq[rel_coords[0] : rel_coords[1]].upper()
        + seq[rel_coords[1] :].lower()
    )


# def highlight_by_relative_coordinates(seq: str, rel_coords: RelativeCoordinates, highlight_method = lambda s: s.upper()):
#     return (
#         seq[: rel_coords[0]]
#         + highlight_method(seq[rel_coords[0] : rel_coords[1]])
#         + seq[rel_coords[1] :]
#     )

DecorationList = list[tuple[RelativeCoordinates, Callable[[str], str]]]


def highlight_by_multiple_relative_coordinates(
    seq: str, decoration_list: DecorationList
):
    identity = lambda x: x

    decoration_list = (
        [((0, 0), identity)] + decoration_list + [((len(seq), len(seq)), identity)]
    )

    output_seq = ""
    for dc, dn in zip(decoration_list, decoration_list[1:]):
        (c1, c2), method = dc[0], dc[1]
        c3 = dn[0][0]
        output_seq += method(seq[c1:c2]) + seq[c2:c3]

    return output_seq


def cut_by_relative_coordinates(seq: str, rel_coords: RelativeCoordinates):
    return seq[rel_coords[0] : rel_coords[1]]


def relative_to(
    inner: RelativeCoordinates, outer: RelativeCoordinates
) -> RelativeCoordinates:
    # for A and B coordinates relative to S, returns coordinates of A relative to S[B].
    return inner[0] - outer[0], inner[1] - outer[0]


def lift_coordinates_to_gapped_alignment(
    ali_seq: str, rel_coords: RelativeCoordinates
) -> RelativeCoordinates:
    # non-zero positions in the seq1 gapped alignment; ali_seq[seq_mapping[i]] equals seq[i]
    seq_mapping = np.where(np.array(list(seq1_ali)) != "-")[0]
    return seq_mapping[rel_coords[0]], seq_mapping[rel_coords[1]]

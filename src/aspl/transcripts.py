from dataclasses import dataclass
from typing import Self, TypeVar
from . import events
import itertools

T = TypeVar("T")


def get_windows(data: list[T], size: int, step: int):
    for i in range(0, len(data) - size + 1, step):
        yield data[i : i + size]


# treats terminal sites as splice sites, is it correct?
@dataclass(frozen=True)
class Transcript:
    splice_sites: list[events.SpliceSite]

    @property
    def splice_site_types(self):
        return itertools.cycle(
            [events.SpliceSiteType.ACCEPTOR, events.SpliceSiteType.DONOR]
        )

    @classmethod
    def from_exons(cls, exons: list[events.Exon]) -> Self:
        return cls(splice_sites=[ss for e in exons for ss in e.splice_sites])

    @classmethod
    def from_splice_sites(cls, splice_sites: list[events.SpliceSite]) -> Self:
        return cls(splice_sites)

    def get_sites_sorted_by_coordinate(self):
        sites = self.splice_sites
        if sites[0].strand == events.Strand.MINUS:
            sites = sites[::-1]
        return sites

    def is_valid(self) -> bool:
        result = True
        result &= len(self.splice_sites) % 2 == 0
        for ss, sst in zip(self.splice_sites, self.splice_site_types):
            result &= ss.type == sst
        result &= len(set(ss.seqname for ss in self.splice_sites)) == 1
        result &= len(set(ss.strand for ss in self.splice_sites)) == 1

        if not result:
            return False

        for s1, s2 in zip(self.splice_sites, self.splice_sites[1:]):
            result &= ((s1.coord <= s2.coord) & (s1.strand == events.Strand.PLUS)) | (
                (s1.coord >= s2.coord) & (s1.strand == events.Strand.MINUS)
            )
        return result

    def iter_exons_w_flanking(self):
        for site_group in get_windows(self.splice_sites[1:-1], 4, 2):
            yield events.CassetteExon.from_splice_sites(site_group)

    def iter_exon_skipping_introns(self):
        for site_group in get_windows(self.splice_sites[1:-1], 4, 2):
            yield events.SpliceJunction(site_group[0], site_group[3])

    def iter_exons(self):
        ss_iter = iter(self.splice_sites)
        ss_pairs = zip(ss_iter, ss_iter)

        for s1, s2 in ss_pairs:
            yield events.Exon(s1, s2)

    def iter_introns(self):
        ss_iter = iter(self.splice_sites[1:-1])
        ss_pairs = zip(ss_iter, ss_iter)

        for s1, s2 in ss_pairs:
            yield events.SpliceJunction(s1, s2)

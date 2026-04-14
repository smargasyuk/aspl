from dataclasses import dataclass
from enum import Enum, StrEnum
from typing import ClassVar, override, Self
from abc import ABC, abstractmethod


class Strand(StrEnum):
    PLUS = "+"
    MINUS = "-"


class SpliceSiteType(StrEnum):
    DONOR = "D"
    ACCEPTOR = "A"
    UNKNOWN = "U"


class CoordinateSystem(Enum):
    ZERO_BASED = 0
    ONE_BASED = 1


# 1-based coordinates internally
@dataclass(frozen=True)
class PointCoordinate:
    seqname: str
    coord: int
    strand: Strand


# 0-based, right-half-open BED-like coordinates internally
@dataclass(frozen=True)
class IntervalCoordinate:
    seqname: str
    seqStart: int
    seqEnd: int
    strand: Strand | None


@dataclass(frozen=True)
class SpliceSite(PointCoordinate):
    type: SpliceSiteType

    @property
    def splice_sites(self):
        return [self]

    def get_sites_sorted_by_coordinate(self):
        return [self]

    def is_valid(self):
        return True


class SplicingSubPath(ABC):
    splice_site_types: ClassVar[list[SpliceSiteType]] = []

    @property
    @abstractmethod
    def splice_sites(self) -> list[SpliceSite]:
        pass

    @classmethod
    def from_splice_sites(cls, splice_sites: list[SpliceSite]) -> Self:
        assert len(splice_sites) == len(cls.splice_site_types)
        return cls(*splice_sites)

    def is_valid(self) -> bool:
        result = True
        result &= len(self.splice_sites) == len(self.splice_site_types)
        for ss, sst in zip(self.splice_sites, self.splice_site_types):
            result &= ss.type == sst
        result &= len(set(ss.seqname for ss in self.splice_sites)) == 1
        result &= len(set(ss.strand for ss in self.splice_sites)) == 1

        if not result:
            return False

        for s1, s2 in zip(self.splice_sites, self.splice_sites[1:]):
            result &= ((s1.coord <= s2.coord) & (s1.strand == Strand.PLUS)) | (
                (s1.coord >= s2.coord) & (s1.strand == Strand.MINUS)
            )
        return result

    def get_sites_sorted_by_coordinate(self):
        sites = self.splice_sites
        if sites[0].strand == Strand.MINUS:
            sites = sites[::-1]
        return sites


# always 1-based as of now
@dataclass(frozen=True)
class SpliceJunction(SplicingSubPath):
    donor_site: SpliceSite
    acceptor_site: SpliceSite
    splice_site_types: ClassVar[list[SpliceSiteType]] = [
        SpliceSiteType.DONOR,
        SpliceSiteType.ACCEPTOR,
    ]

    @property
    @override
    def splice_sites(self) -> list[SpliceSite]:
        return [self.donor_site, self.acceptor_site]

    def get_interval(self):
        assert self.is_valid()
        coords = [s.coord for s in self.get_sites_sorted_by_coordinate()]

        return IntervalCoordinate(
            self.donor_site.seqname, coords[0], coords[1] - 1, self.donor_site.strand
        )


@dataclass(frozen=True)
class Exon(SplicingSubPath):
    siteB: SpliceSite
    siteC: SpliceSite
    splice_site_types: ClassVar[list[SpliceSiteType]] = [
        SpliceSiteType.ACCEPTOR,
        SpliceSiteType.DONOR,
    ]

    @property
    @override
    def splice_sites(self) -> list[SpliceSite]:
        return [self.siteB, self.siteC]

    def get_interval(self):
        assert self.is_valid()
        coords = [s.coord for s in self.get_sites_sorted_by_coordinate()]

        return IntervalCoordinate(
            self.siteB.seqname, coords[0] - 1, coords[1], self.siteB.strand
        )


@dataclass(frozen=True)
class CassetteExon(SplicingSubPath):
    siteA: SpliceSite
    siteB: SpliceSite
    siteC: SpliceSite
    siteD: SpliceSite
    splice_site_types: ClassVar[list[SpliceSiteType]] = [
        SpliceSiteType.DONOR,
        SpliceSiteType.ACCEPTOR,
        SpliceSiteType.DONOR,
        SpliceSiteType.ACCEPTOR,
    ]

    @property
    @override
    def splice_sites(self) -> list[SpliceSite]:
        return [self.siteA, self.siteB, self.siteC, self.siteD]

    def get_flanking_junctions(self):
        return SpliceJunction(self.siteA, self.siteB), SpliceJunction(
            self.siteC, self.siteD
        )

    def get_exon(self):
        return Exon(self.siteB, self.siteC)

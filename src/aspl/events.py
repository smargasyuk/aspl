from dataclasses import dataclass
from enum import Enum, StrEnum
from typing import Any, Literal


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


# may be 0-based
@dataclass(frozen=True)
class SpliceSite(PointCoordinate):
    type: SpliceSiteType

    def format(
        self,
        cs: CoordinateSystem = CoordinateSystem.ONE_BASED,
        include_type: bool = False,
    ):
        coord = self.coord - 1 if cs == CoordinateSystem.ZERO_BASED else self.coord
        result = self.seqname + "_" + str(coord) + "_" + self.strand

        if include_type:
            result += "_" + self.type

        return result

    @staticmethod
    def parse(
        s: str,
        cs: CoordinateSystem = CoordinateSystem.ONE_BASED,
        type_included: bool = False,
    ):
        if type_included:
            seqname, coord, strand, site_type = s.rsplit("_", maxsplit=3)
            site_type = SpliceSiteType(site_type)
        else:
            seqname, coord, strand = s.rsplit("_", maxsplit=2)
            site_type = SpliceSiteType.UNKNOWN

        coord = int(coord)
        if cs == CoordinateSystem.ZERO_BASED:
            coord += 1

        strand = Strand(strand)

        return SpliceSite(seqname, coord, strand, site_type)


# always 1-based as of now
@dataclass(frozen=True)
class SpliceJunction:
    donor_site: SpliceSite
    acceptor_site: SpliceSite

    def is_valid(self):
        result = True
        result &= self.donor_site.type == SpliceSiteType.DONOR
        result &= self.acceptor_site.type == SpliceSiteType.ACCEPTOR
        result &= self.donor_site.seqname == self.acceptor_site.seqname
        result &= self.donor_site.strand == self.acceptor_site.strand

        if not result:
            return False

        strand = self.donor_site.strand

        if strand == Strand.PLUS:
            result &= self.acceptor_site.coord > self.donor_site.coord
        else:
            result &= self.acceptor_site.coord < self.donor_site.coord
        return result

    @staticmethod
    def parse(s: str):
        seqname, coord1, coord2, strand = s.rsplit("_", maxsplit=3)
        strand = Strand(strand)
        coord1, coord2 = int(coord1), int(coord2)
        if strand == Strand.MINUS:
            coord2, coord1 = coord1, coord2

        return SpliceJunction(
            SpliceSite(seqname, coord1, strand, SpliceSiteType.DONOR),
            SpliceSite(seqname, coord2, strand, SpliceSiteType.ACCEPTOR),
        )

    def format(self):
        assert self.is_valid()
        coord1, coord2 = self.donor_site.coord, self.acceptor_site.coord
        if self.donor_site.strand == Strand.MINUS:
            coord2, coord1 = coord1, coord2
        return (
            self.donor_site.seqname
            + "_"
            + str(coord1)
            + "_"
            + str(coord2)
            + "_"
            + self.donor_site.strand
        )


@dataclass(frozen=True)
class Exon:
    siteB: SpliceSite
    siteC: SpliceSite

    def get_splice_sites(self):
        return [self.siteB, self.siteC]

    def is_valid(self) -> bool:
        is_valid = True
        is_valid &= self.siteC.type == SpliceSiteType.DONOR
        is_valid &= self.siteB.type == SpliceSiteType.ACCEPTOR
        is_valid &= self.siteB.seqname == self.siteC.seqname
        is_valid &= self.siteB.strand == self.siteC.strand

        if not is_valid:
            return False

        strand = self.siteB.strand

        if strand == Strand.PLUS:
            is_valid &= self.siteB.coord < self.siteC.coord
        else:
            is_valid &= self.siteB.coord > self.siteC.coord
        return is_valid

    def format(self):
        assert self.is_valid()
        coords = [str(s.coord) for s in self.get_splice_sites()]
        if self.siteB.strand == Strand.MINUS:
            coords = coords[::-1]

        return self.siteB.seqname + "_" + "_".join(coords) + "_" + self.siteB.strand


@dataclass(frozen=True)
class CassetteExon:
    siteA: SpliceSite
    siteB: SpliceSite
    siteC: SpliceSite
    siteD: SpliceSite

    def get_splice_sites(self):
        return [self.siteA, self.siteB, self.siteC, self.siteD]

    def get_flanking_junctions(self):
        return SpliceJunction(self.siteA, self.siteB), SpliceJunction(
            self.siteC, self.siteD
        )

    def get_exon(self):
        return Exon(self.siteB, self.siteC)

    def is_valid(self):
        is_valid = True
        is_valid &= self.siteA.type == self.siteC.type == SpliceSiteType.DONOR
        is_valid &= self.siteB.type == self.siteD.type == SpliceSiteType.ACCEPTOR
        is_valid &= (
            self.siteA.seqname
            == self.siteB.seqname
            == self.siteC.seqname
            == self.siteD.seqname
        )
        is_valid &= (
            self.siteA.strand
            == self.siteB.strand
            == self.siteC.strand
            == self.siteD.strand
        )

        if not is_valid:
            return False

        strand = self.siteA.strand

        if strand == Strand.PLUS:
            is_valid &= (
                self.siteA.coord
                < self.siteB.coord
                < self.siteC.coord
                < self.siteD.coord
            )
        else:
            is_valid &= (
                self.siteA.coord
                > self.siteB.coord
                > self.siteC.coord
                > self.siteD.coord
            )
        return is_valid

    @staticmethod
    def parse(s: str):
        seqname, *coords, strand = s.rsplit("_", maxsplit=5)
        strand = Strand(strand)
        coords = [int(c) for c in coords]
        if strand == Strand.MINUS:
            coords = coords[::-1]

        return CassetteExon(
            SpliceSite(seqname, coords[0], strand, SpliceSiteType.DONOR),
            SpliceSite(seqname, coords[1], strand, SpliceSiteType.ACCEPTOR),
            SpliceSite(seqname, coords[2], strand, SpliceSiteType.DONOR),
            SpliceSite(seqname, coords[3], strand, SpliceSiteType.ACCEPTOR),
        )

    def format(self):
        assert self.is_valid()
        coords = [str(s.coord) for s in self.get_splice_sites()]
        if self.siteA.strand == Strand.MINUS:
            coords = coords[::-1]

        return self.siteA.seqname + "_" + "_".join(coords) + "_" + self.siteA.strand

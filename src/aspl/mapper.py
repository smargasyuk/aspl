from calendar import c
import dataclasses
from abc import ABC, abstractmethod
from typing import TypeVar, override
import polars as pl
from .events import (
    SpliceSite,
    SpliceJunction,
    CassetteExon,
    SpliceSiteType,
    CoordinateSystem,
    Exon,
    Strand,
)
from .transcripts import Transcript
from .formatters import UnderscoreSeparated
from liftover import ChainFile, get_lifter

T = TypeVar("T", CassetteExon, SpliceJunction, Exon, Transcript)


class Mapper(ABC):
    @abstractmethod
    def map_splice_site(self, ss: SpliceSite) -> SpliceSite | None:
        pass

    def map_event(self, event: T) -> T | None:
        mapped_ss = [self.map_splice_site(ss) for ss in event.splice_sites]
        if any(mss is None for mss in mapped_ss):
            return None
        cls = type(event)
        mapped_event = cls.from_splice_sites(mapped_ss)
        if not mapped_event.is_valid():
            return None
        return mapped_event


@dataclasses.dataclass
class SiteTableMapper(Mapper):
    mapping: dict[
        SpliceSite, SpliceSite
    ]  # all splice sites in the dictionary have to be untyped

    @override
    def map_splice_site(self, ss: SpliceSite):
        return dataclasses.replace(
            self.mapping[dataclasses.replace(ss, type=SpliceSiteType.UNKNOWN)],
            type=ss.type,
        )

    def map_junction(self, sj: SpliceJunction):
        return SpliceJunction(
            self.map_splice_site(sj.donor_site), self.map_splice_site(sj.acceptor_site)
        )

    def map_cassette_exon(self, ce: CassetteExon):
        return CassetteExon(
            self.map_splice_site(ce.siteA),
            self.map_splice_site(ce.siteB),
            self.map_splice_site(ce.siteC),
            self.map_splice_site(ce.siteD),
        )

    def map_exon(self, e: Exon):
        return Exon(
            self.map_splice_site(e.siteB),
            self.map_splice_site(e.siteC),
        )

    # parse the zero-based mapping table built by maptools
    @staticmethod
    def parse_site_file(dfm1: pl.DataFrame):
        formatter = UnderscoreSeparated(CoordinateSystem.ZERO_BASED)
        m1_dict = {}
        for r in dfm1.iter_rows(named=True):
            s1 = formatter.parse(r["S1"], SpliceSite)
            s2 = formatter.parse(r["S2"], SpliceSite)
            m1_dict[s1] = s2
        return SiteTableMapper(m1_dict)


def _map_single_site_with_pyliftover(ss: SpliceSite, chain: ChainFile):
    strand_inverted = {"+": "-", "-": "+"}
    targets_list = chain[ss.seqname][ss.coord]
    if len(targets_list) != 1:
        return None
    new_chrom, new_coord, new_strand = targets_list[0]
    # always converts the plus strand; if the original strand is minus, invert the target strand
    if ss.strand == Strand.MINUS:
        new_strand = strand_inverted[new_strand]
    new_strand = Strand(new_strand)
    return SpliceSite(new_chrom, new_coord, new_strand, ss.type)


@dataclasses.dataclass
class PyLiftOverMapper(Mapper):
    converter: ChainFile

    @override
    def map_splice_site(self, ss: SpliceSite) -> SpliceSite | None:
        return _map_single_site_with_pyliftover(ss, self.converter)

    @staticmethod
    def parse_chain_file(chain_file_path: str):
        chain_file = ChainFile(chain_file_path, one_based=True)
        return PyLiftOverMapper(chain_file)


@dataclasses.dataclass
class ReciprocalBestHitMapper(Mapper):
    c1: ChainFile
    c2: ChainFile

    @override
    def map_splice_site(self, ss: SpliceSite) -> SpliceSite | None:
        target_site = _map_single_site_with_pyliftover(ss, self.c1)
        if target_site is None:
            return None
        target_backmapping = _map_single_site_with_pyliftover(target_site, self.c2)
        if target_backmapping != ss:
            return None
        return target_site

    @staticmethod
    def parse_chain_files(chain_file1_path: str, chain_file2_path: str):
        c1 = ChainFile(chain_file1_path, one_based=True)
        c2 = ChainFile(chain_file2_path, one_based=True)
        return ReciprocalBestHitMapper(c1, c2)

    @staticmethod
    def from_genomes_pair(genome1, genome2):
        c1 = get_lifter(genome1, genome2, one_based=True)
        c2 = get_lifter(genome2, genome1, one_based=True)
        return ReciprocalBestHitMapper(c1, c2)

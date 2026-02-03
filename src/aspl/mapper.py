import dataclasses
import polars as pl
from .events import (
    SpliceSite,
    SpliceJunction,
    CassetteExon,
    SpliceSiteType,
    CoordinateSystem,
    Exon,
)
from .formatters import UnderscoreSeparated


@dataclasses.dataclass
class SiteMapper:
    mapping: dict[
        SpliceSite, SpliceSite
    ]  # all splice sites in the dictionary have to be untyped

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
        return SiteMapper(m1_dict)

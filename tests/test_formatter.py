from aspl.events import CassetteExon, SpliceJunction, SpliceSite, CoordinateSystem
from aspl.formatters import UnderscoreSeparated

test_ce = "chr1_10_20_30_40_-"
formatter = UnderscoreSeparated()


def test_cassette_exon():
    assert test_ce == formatter.format(formatter.parse(test_ce, CassetteExon))


def test_junction():
    sj1 = "chr1_10_40_-"
    assert sj1 == formatter.format(formatter.parse(sj1, SpliceJunction))


def test_site():
    ss1 = "chr1_10_-"
    assert ss1 == formatter.format(formatter.parse(ss1, SpliceSite))


def test_site_ipsa_format():
    ss1 = "chr1_10_-_D"
    formatter = UnderscoreSeparated(include_site_type=True)
    assert ss1 == formatter.format(formatter.parse(ss1, SpliceSite))


def test_site_coordinate_systems():
    ss1 = "chr1_9_-_D"
    formatter = UnderscoreSeparated(
        coordinate_system=CoordinateSystem.ZERO_BASED, include_site_type=True
    )
    assert formatter.parse(ss1, SpliceSite).coord == 10


def test_ce_site_relation():
    ss1 = "chr1_40_-_D"
    formatter = UnderscoreSeparated(include_site_type=True)
    assert formatter.parse(test_ce, CassetteExon).siteA == formatter.parse(
        ss1, SpliceSite
    )

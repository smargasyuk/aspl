from aspl.events import CassetteExon, SpliceJunction, SpliceSite, CoordinateSystem

test_ce = "chr1_10_20_30_40_-"


def test_cassette_exon():
    assert test_ce == CassetteExon.parse(test_ce).format()


def test_junction():
    sj1 = "chr1_10_40_-"
    assert sj1 == SpliceJunction.parse(sj1).format()


def test_site():
    ss1 = "chr1_10_-"
    assert ss1 == SpliceSite.parse(ss1).format()


def test_site_ipsa_format():
    ss1 = "chr1_10_-_D"
    assert ss1 == SpliceSite.parse(ss1, type_included=True).format(include_type=True)


def test_site_coordinate_systems():
    ss1 = "chr1_9_-_D"
    assert (
        SpliceSite.parse(ss1, type_included=True, cs=CoordinateSystem.ZERO_BASED).coord
        == 10
    )


def test_ce_site_relation():
    ss1 = "chr1_40_-_D"
    assert CassetteExon.parse(test_ce).siteA == SpliceSite.parse(
        ss1, type_included=True
    )

from aspl.events import CassetteExon, SpliceJunction, SpliceSite


def test_cassette_exon():
    ce1 = "chr1_10_20_30_40_-"
    assert ce1 == CassetteExon.parse(ce1).format()


def test_junction():
    sj1 = "chr1_10_40_-"
    assert sj1 == SpliceJunction.parse(sj1).format()


def test_site():
    ss1 = "chr1_10_-"
    assert ss1 == SpliceSite.parse(ss1).format()


def test_site_ipsa_format():
    ss1 = "chr1_10_-_D"
    assert ss1 == SpliceSite.parse(ss1, type_included=True).format(include_type=True)

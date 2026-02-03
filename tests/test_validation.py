import dataclasses as dt
from aspl.events import (
    SpliceSite,
    SpliceJunction,
    CassetteExon,
    Exon,
    Strand,
    SpliceSiteType,
)

import pytest

ce1 = CassetteExon.parse("chr1_10_20_30_40_-")
jxn1 = SpliceJunction(ce1.siteA, ce1.siteB)


ce_test_data = [
    (ce1, True),
    (
        dt.replace(ce1, siteA=dt.replace(ce1.siteA, strand=Strand.PLUS)),
        False,
    ),  # invalid strand
    (dt.replace(ce1, siteB=dt.replace(ce1.siteB, coord=5)), False),  # invalid order
    (
        dt.replace(ce1, siteB=dt.replace(ce1.siteB, seqname="chr2")),
        False,
    ),  # different chromosomes
    (
        dt.replace(ce1, siteB=dt.replace(ce1.siteC, type=SpliceSiteType.ACCEPTOR)),
        False,
    ),  # wrong site type
]

jxn_test_data = [
    (jxn1, True),
    (
        dt.replace(jxn1, donor_site=dt.replace(jxn1.donor_site, strand=Strand.PLUS)),
        False,
    ),  # invalid strand
    (
        dt.replace(jxn1, donor_site=dt.replace(jxn1.donor_site, coord=5)),
        False,
    ),  # invalid order
    (
        dt.replace(jxn1, donor_site=dt.replace(jxn1.donor_site, seqname="chr2")),
        False,
    ),  # different chromosomes
    (
        dt.replace(
            jxn1, donor_site=dt.replace(jxn1.donor_site, type=SpliceSiteType.ACCEPTOR)
        ),
        False,
    ),  # wrong site type
]


@pytest.mark.parametrize("instance, expected", ce_test_data)
def test_ce_validation(instance: CassetteExon, expected: bool):
    assert instance.is_valid() == expected


@pytest.mark.parametrize("instance, expected", jxn_test_data)
def test_jxn_validation(instance: SpliceJunction, expected: bool):
    assert instance.is_valid() == expected

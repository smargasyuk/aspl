from dataclasses import replace
from . import events, transcripts


def _stranded_sites_leq(a: events.SpliceSite, b: events.SpliceSite):
    if a.strand == events.Strand.PLUS:
        return a.coord <= b.coord
    else:
        return a.coord >= b.coord


def _cut_exon(e: events.Exon, s: events.SpliceSite | None, side: str):
    if s is None:
        return e
    if side == "left":
        if _stranded_sites_leq(e.siteC , s):
            return e
        if not _stranded_sites_leq(e.siteB , s): 
            return None
        return replace(e, siteC = replace(s, type = events.SpliceSiteType.DONOR))

    if side == "right":
        if _stranded_sites_leq(s, e.siteB):
            return e
        if not _stranded_sites_leq(s, e.siteC): 
            return None
        return replace(e, siteB = replace(s, type = events.SpliceSiteType.ACCEPTOR))
    
    raise NotImplementedError


def _cut_exon_twosites(e: events.Exon, s1: events.SpliceSite | None, s2: events.SpliceSite | None):
    ec1 = _cut_exon(e, s1, "right")
    if ec1 is None:
        return None
    ec2 = _cut_exon(ec1, s2, "left")
    return ec2
    

def cut_transcript(t: transcripts.Transcript, s1: events.SpliceSite | None, s2: events.SpliceSite | None):
    exons_new = [_cut_exon_twosites(e, s1, s2) for e in t.iter_exons()]
    return transcripts.Transcript.from_exons([e for e in exons_new if e is not None])


def insert_exon_into_transcript(t: transcripts.Transcript, e: events.Exon):
    strand = e.siteB.strand
    sorting_key = lambda x: x.coord if strand == events.Strand.PLUS else -x.coord
    return transcripts.Transcript(sorted(t.splice_sites + [e.siteB, e.siteC], key=sorting_key))
import aspl
from aspl import seq_utils
from pathlib import Path
from pysam import FastaFile


def get_genome_fasta():
    data_dir = Path(__file__).parent / "data"
    fasta_path = data_dir / "chr21.fa.gz"
    return FastaFile(str(fasta_path))


def test_exon_sequence():
    test_exon = aspl.formatters.UnderscoreSeparated().parse(
        "chr21_10541112_10541165_+", aspl.events.Exon
    )
    fasta = get_genome_fasta()
    exon_seq = seq_utils.get_interval_sequence(test_exon.get_interval(), fasta)
    assert exon_seq == "TCCTGATCCGACTGACCTGGCGGGAGTCATCATTGAGCTCGGCCCCAATGACAG"


def test_exon_translate():
    exon_seq = "TCCTGATCCGACTGACCTGGCGGGAGTCATCATTGAGCTCGGCCCCAATGACAG"
    protein_seq, leftover = seq_utils._translate_exon_with_splice_marks(exon_seq, "AG")
    assert protein_seq == "Spdptdlagviielgpnd"
    assert leftover == "AG"


def test_short_exon_translate():
    exon_seq = "AT"
    protein_seq, leftover = seq_utils._translate_exon_with_splice_marks(exon_seq, "")
    assert protein_seq == ""
    assert leftover == "AT"


def test_star_uppercase():
    exon_seq = "AAAAAAA"
    protein_seq, leftover = seq_utils._translate_exon_with_splice_marks(exon_seq, "T")
    assert protein_seq == "#k"
    assert leftover == "AA"


def test_transcript_nice_translate():
    # POFUT2-203 CDS
    transcript_str = "chr21_45265485_45265635_45267590_45267713_45269839_45270019_45277017_45277142_45278103_45278169_45282349_45282459_45283383_45283527_45285678_45285928_45287741_45287871_-"
    fasta = get_genome_fasta()
    transcript = aspl.formatters.UnderscoreSeparated().parse(
        transcript_str, aspl.transcripts.Transcript
    )
    transcript_coding_seq = seq_utils.translate_with_splice_marks(transcript, fasta)

    print(transcript_coding_seq)
    transcript_n_exons = len(list(transcript.iter_exons()))
    count_capitals = lambda s: sum(1 for c in s if c.isupper())
    assert transcript_coding_seq.lower().startswith(
        "MATLSFVFLLLGAVSWPPASASGQEFWPGQSAADILSGAASRRRYLLYDVNPPEGFNLRR".lower()
    )
    assert count_capitals(transcript_coding_seq) >= (transcript_n_exons - 1)
    assert "*" not in transcript_coding_seq


def test_poison_negative():
    seq = "aaaaGAaaaaAaaa*" + "ae" * 20
    assert not seq_utils.is_poison_by_50nt(seq)


def test_poison_positive():
    seq = "aaaaGAaaaaAaaa*" + "ae" * 20 + "Paaa"
    assert seq_utils.is_poison_by_50nt(seq)


def test_poison_star():
    seq = "aaaaGAaaaaAaaa*" + "ae" * 20 + "#aaa"
    assert seq_utils.is_poison_by_50nt(seq)

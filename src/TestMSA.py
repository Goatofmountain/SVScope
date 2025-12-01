from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline, ClustalOmegaCommandline, MafftCommandline
from Bio import AlignIO, SeqIO
import tempfile, subprocess, io, os


def msa_from_list_stdin(seq_records, engine="muscle"):
    """Run expensive MSA process"""
    # move fasta file into memory 
    fasta_str = io.StringIO()
    SeqIO.write(seq_records, fasta_str, "fasta")
    if engine == "muscle":
        cmd = ["muscle", "-quiet"]
    elif engine == "clustalo":
        cmd = ["clustalo", "--infile=-", "--outfile=-", "--force"]
    elif engine == "mafft":
        cmd = ["mafft", "--quiet", "--thread", "1", "-"]
    else:
        raise ValueError("stdin mode only support muscle/clustalo")
    proc = subprocess.run(cmd,
                          input=fasta_str.getvalue(),
                          text=True,
                          capture_output=True,
                          check=True)
    aln = AlignIO.read(io.StringIO(proc.stdout), "fasta")
    seqList = [str(rec.seq).upper() for rec in aln]
    return seqList

# usage
seq_list = [
    SeqRecord(Seq("ATGGCATACGCTAGCTAGCTAG"), id="Alpha"),
    SeqRecord(Seq("ATGGCATACGCTAGCTAGCTAG"), id="Beta"),
    SeqRecord(Seq("ATCGTACGCTAGCTAGCTAG"), id="Gamma"),
]
aln = msa_from_list_stdin(seq_list, engine='mafft')
print(aln)

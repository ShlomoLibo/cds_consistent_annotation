from os import PathLike
import os

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbitblastxCommandline


def tblastx(
        query: PathLike,
        db: PathLike,
        out: PathLike,
        blast_start=None,
        blast_end=None,
        outfmt=5,
        evalue=0.01,
        override=True
):
    """
    check BLAST documentation: https://www.ncbi.nlm.nih.gov/books/NBK279684/
    :param override: whether to do blast even if a file with the specified name already exists
    """
    if override or not os.path.exists(out):

        # Create BLAST database

        cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=db)
        os.popen(str(cline)).read()

        # Run BLAST
        if blast_start is not None and blast_end is not None:
            cline = NcbitblastxCommandline(query=query, db=db, out=out,
                                           outfmt=outfmt, evalue=evalue, query_loc=f"{blast_start}-{blast_end}")
        else:
            cline = NcbitblastxCommandline(query=query, db=db, out=out,
                                           outfmt=outfmt, evalue=evalue)

        os.popen(str(cline)).read()


def cd_hit_est(
        input_fasta,
        output_file,
        sequence_identity,
        word_size
):
    """
    check cd-hit documentation: https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide
    """
    cline = f"cd-hit-est -i {input_fasta} -o {output_file} -c {sequence_identity} -n {word_size} -d 100"
    os.popen(cline).read()
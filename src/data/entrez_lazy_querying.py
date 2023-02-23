from collections import defaultdict
from typing import Dict, List
import logging

from Bio import Entrez, SeqIO, Seq

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)


class LazySeq:
    """
    A class that represents sequence data, that is designed to aggregate the sequences that are needed to
    retrieve (e.g. from a server) and retrieve them all in once, to minimize requests to the server.
    uses e-utilities.
    """
    def get_sequence(self) -> Seq:
        """
        When one sequence is needed evaluation, all sequences are requested from the server.
        :return: The actual sequence.
        """


class LazyGenomeSeq(LazySeq):
    completed_queries = dict()
    mrna_start: int  # the starting position of the mrna (according to the gff file)
    mrna_end: int  # the starting position of the mrna (according to the gff file)
    genome_id: str  # the chromosome id (in NCBI) on which the mrna located
    genome_to_pending_queries: Dict[str, List['LazyGenomeSeq']] = defaultdict(list)

    def __init__(self, email, start, end, genome):
        self.start = start
        self.end = end
        Entrez.email = email
        self.genome_id = genome
        self.query_id = (genome, start, end)
        LazyGenomeSeq.genome_to_pending_queries[genome].append(self)

    def get_sequence(self):
        if LazyGenomeSeq.genome_to_pending_queries:
            for genome, pending_queries in LazyGenomeSeq.genome_to_pending_queries.items():
                start = min(pending_queries, key=lambda x: x.start).start
                end = max(pending_queries, key=lambda x: x.end).end
                handle = Entrez.efetch(db="nucleotide", id=genome, rettype="fasta", retmode="text", seq_start=start,
                                       seq_stop=end)

                logger = logging.getLogger(__name__)
                logger.info(f"made entrez query with id: {genome}, start: {start}, end: {end}")
                record = next(SeqIO.parse(handle, "fasta"))
                for query in pending_queries:
                    LazyGenomeSeq.completed_queries[query.query_id] = record.seq[
                                                                    query.start - start: query.end - start]
        LazyGenomeSeq.genome_to_pending_queries = defaultdict(list)
        return LazyGenomeSeq.completed_queries[self.query_id]


class LazyMRNASeq(LazySeq):
    completed_queries = dict()
    mrna_start: int  # the starting position of the mrna (according to the gff file)
    mrna_end: int  # the starting position of the mrna (according to the gff file)
    query_id: str
    genome_seq: LazyGenomeSeq = None

    def __init__(self, query_id, email, mrna_start, mrna_end, genome):
        self.query_id = query_id
        self.mrna_start = mrna_start
        self.mrna_end = mrna_end
        Entrez.email = email
        self.genome_seq = LazyGenomeSeq(email, mrna_start, mrna_end, genome)

    def get_sequence(self):
        """
        # The following code is for querying ncbi for mrnas directly, unfortunately doesn't work, since they don't 
        # contain introns (although, it might be possible to do this approach to recover CDS, by inferring which residues belong
        # to which cds piece.
        query_ids = [lazy_mrna.query_id for lazy_mrna in LazyMRNASeq.pending_queries]
        if query_ids:
            print(f"made query with ids: {','.join(query_ids)}")
            handle = Entrez.efetch(db="nucleotide", id=",".join(query_ids), rettype="fasta", retmode="text")
            for record in SeqIO.parse(handle, "fasta"):
                LazyMRNASeq.completed_queries[record.id] = record.seq
                LazyMRNASeq.pending_queries = set()
        return LazyMRNASeq.completed_queries[self.query_id]
        """

        # second try: query per genome
        return self.genome_seq.get_sequence()


class LazyExonSeq(LazySeq):
    """
    Lazy exon depends on a lazyMRNA, since exons / cds pieces don't have unique ids in genbank.
    """
    lazy_mrna: LazyMRNASeq
    exon_start: int  # the starting position of the exon (according to the gff file)
    exon_end: int  # the starting position of the exon (according to the gff file)

    def __init__(self, lazy_mrna, exon_start, exon_end):
        self.lazy_mrna = lazy_mrna
        self.exon_start = exon_start
        self.exon_end = exon_end

    def get_sequence(self):
        return self.lazy_mrna.get_sequence()[self.exon_start - self.lazy_mrna.mrna_start: self.exon_end - self.lazy_mrna.mrna_start]

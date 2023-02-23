from Bio.Data import CodonTable
from typing import List
from Bio.Seq import Seq
import warnings
from pyfaidx import Fasta
import gffutils
from gffutils import FeatureDB
from visualization.visualize import VisualizableRange
from data.entrez_lazy_querying import LazyMRNASeq, LazyExonSeq
from visualization.range_tools import Range


class Indexable:
    instances = {}

    def __init__(self, id_, dont_index=False):
        if self.__class__.get(id_) is not None:
            warnings.warn(f"class {self.__class__.__name__}: object of id: {id_} already exists")
        self.id = id_
        if not dont_index:
            self.__class__.instances[id_] = self

    @classmethod
    def get(cls, id_):
        return cls.instances.get(id_)

    @classmethod
    def delete(cls, id_):
        if id_ in cls.instances:
            del cls.instances[id_]

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id


class ExonAnnotation(Indexable, VisualizableRange):
    start: int
    end: int
    _cluster: 'ExonCluster' = None
    _cds: 'CDS' = None
    id_: str
    lazy_seq: LazyExonSeq = None

    def __init__(self, start, end, cds, exon_name, cluster=None, **kwargs):
        Range.__init__(self, start, end)
        self.cds = cds
        self.exon_name = exon_name
        self.cluster = cluster
        self.subject = self.create_exon_subject(cds, exon_name)
        Indexable.__init__(self, id_=self.subject, **kwargs)
        self.lazy_seq = LazyExonSeq(self.mrna.lazy_seq, start, end)
        self.group = self.genome.name

    @property
    def mrna(self):
        return self.cds.mrna

    @property
    def gene(self):
        return self.cds.gene

    @property
    def genome(self):
        return self.cds.genome

    @property
    def cluster(self):
        return self._cluster

    @cluster.setter
    def cluster(self, cluster):
        if self._cluster:
            self._cluster.exons.remove(self)
        self._cluster = cluster
        if cluster is not None:
            cluster.exons.append(self)

    @property
    def cds(self):
        return self._cds

    @cds.setter
    def cds(self, cds):
        if self._cds:
            self._cds.exons.remove(self)
        self._cds = cds
        if cds is not None:
            cds.exons.append(self)

    @classmethod
    def create_exon_subject(cls, cds: 'CDS', exon_name):
        """
        Used to work with commandline applications such as blast and cd-hit. Used in conjunction with parse_exon_text
        """
        return f"{cds.genome.name}__{cds.gene.gene_name}__{exon_name}"

    @classmethod
    def parse_exon_subject(cls, subject):
        """
            Used to work with commandline applications such as blast and cd-hit. Used in conjunction with crete_exon_subject
        """
        genome, gene, exon = tuple(subject.split("__"))
        return AnnotatedGenome.get(genome), Gene.get(gene), ExonAnnotation.get(subject)

    @property
    def color(self):
        return self.cluster.color

    @property
    def marker(self):
        return self.gene.marker

    @property
    def vertical_offset(self):
        return self.gene.display_position * 0.15

    @property
    def label(self):
        return self.genome.name.split("_")[0]


class ExonCluster(Indexable):
    """
    Represents a cluster of similar exons (from a sequence perspective)
    """
    color: str  # color in which the cluster will be visualized
    id: int # identifier in the .clstr file

    def __init__(self, id_):
        Indexable.__init__(self, id_)
        self._exons: List[ExonAnnotation] = []

    @property
    def exons(self) -> List[ExonAnnotation]:
        return self._exons


class Gene(Indexable):
    gene_name: str
    marker: str  # matplotlib marker for visualization
    display_position: int  # the vertical display position of the gene in the visualization
    _mrnas: List['MRNA'] = None

    def __init__(self, gene_name, marker=None, display_position=None):
        Indexable.__init__(self, id_=gene_name)
        self.gene_name=gene_name
        self.marker=marker
        self.display_position=display_position
        self._mrnas = list()

    @property
    def mrnas(self):
        return self._mrnas


class ExonAlignment(VisualizableRange):
    def __init__(self, start, end, source_exon: ExonAnnotation, target_genome=None):
        Range.__init__(self, start, end)
        self.source_exon = source_exon
        self.target_genome=target_genome
        self.group = self.source_exon.genome.name

    @property
    def color(self):
        return self.source_exon.color

    @property
    def marker(self):
        return self.source_exon.marker

    @property
    def vertical_offset(self):
        return self.source_exon.vertical_offset

    @property
    def label(self):
        return self.source_exon.label


class CDS(Indexable):
    _mrna: 'MRNA' = None

    def __init__(self, id_):
        Indexable.__init__(self, id_)
        self._exons: List[ExonAnnotation] = []

    @property
    def exons(self) -> List[ExonAnnotation]:
        return self._exons

    @classmethod
    def get_cds_from_exon_name(cls, exon_name):
        cds_id = "_".join(exon_name.split("_")[:-1])  # decoding cds_key function output inside get_or_create_db
        cds = cls.get(cds_id)
        if cds is None:
            return cls(cds_id)
        else:
            return cds

    def nucl_seq(self, reverse_complement=False):
        sequence = list()
        for exon in sorted(self._exons, key=lambda e: e.start):
            sequence.append(exon.lazy_seq.get_sequence())
        sequence = Seq(''.join(str(s) for s in sequence))
        if reverse_complement:
            sequence = sequence.reverse_complement()
        return sequence

    def aa_seq(self, reverse_complement=False):
        dna_seq = self.nucl_seq(reverse_complement=reverse_complement)
        # use the seq2seq function to translate the DNA sequence into its corresponding amino acid sequence
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        aa_seq = dna_seq.translate(table=standard_table)
        return aa_seq

    @property
    def gene(self):
        return self.mrna.gene

    @property
    def mrna(self):
        return self._mrna

    @mrna.setter
    def mrna(self, mrna):
        if self.mrna:
            self.mrna.cdss.remove(self)
        self._mrna = mrna
        if mrna is not None:
            mrna.cdss.append(self)


    @property
    def genome(self):
        return self.mrna.genome


class MRNA(Indexable, Range):
    _gene: Gene = None
    _cdss: List['CDS'] = None
    genome: 'AnnotatedGenome' = None
    lazy_seq: LazyMRNASeq

    def __init__(self, id_, start, end, email, genome):
        """
        :param email:  an email is needed to retrieve the sequence from genbank
        """
        Indexable.__init__(self, id_)
        Range.__init__(self, start, end)
        self._cdss = list()
        self.genome = genome
        if id_.startswith("rna-"):
            self.lazy_seq = LazyMRNASeq(id_[4:], email, mrna_start=start, mrna_end=end, genome=genome.chromosome_name)
        else:
            self.lazy_seq = LazyMRNASeq(id_, email, mrna_start=start, mrna_end=end, genome=genome.chromosome_name)

    @property
    def cdss(self):
        return self._cdss

    @property
    def gene(self):
        return self._gene

    @gene.setter
    def gene(self, gene):
        if self._gene:
            self._gene.mrnas.remove(self)
        self._gene = gene
        if gene is not None:
            gene.mrnas.append(self)


class AnnotatedGenome(Indexable):
    name: str
    gff_file: str
    fasta_file: Fasta
    gff_db: FeatureDB
    gff_db_file: str
    chromosome_name: str

    def __init__(self, name, gff_db_file, chromosome_name=""):
        Indexable.__init__(self, id_=name)
        self.name = name
        self.gff_db_file = gff_db_file
        self.gff_db = gffutils.FeatureDB(self.gff_db_file)
        self.chromosome_name = chromosome_name

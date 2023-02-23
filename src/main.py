#!/usr/bin/env python
# coding: utf-8
import argparse
from collections import defaultdict
from itertools import groupby
from os import PathLike
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import FeatureNotFoundError
from Bio.Blast import NCBIXML
from typing import Iterable
from visualization.visualize import visualize_blast_results, visualize_sets, align_sets
from cluster_utils import colorize_clusters, get_cluster_map
from external_tools import tblastx
from data.entrez_lazy_querying import LazyGenomeSeq
from annotations import *
import logging


# has to be one of the following: exon or CDS
# CDS would run the algorithm for annotated CDS pieces
# exon would run the algorithm for exons
FEATURE_TYPE = 'CDS'


log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)


def get_exon_positions(blast_file: PathLike) -> List[List[ExonAlignment]]:
    """Get the positions of the exons for the given gene in the given BLAST output file.

    Args:
        blast_file: The name of the BLAST output file.

    Returns:
         A list of tuples, where each tuple contains the start and end positions of an exon.
    """
    # Open the file
    with open(blast_file, "r") as fh:
        # Parse the file
        blast_records = NCBIXML.parse(fh)

        # Create a list to store the positions of the exons
        exon_position_dict = defaultdict(list)

        # Iterate over the Blast records
        for blast_record in blast_records:
            # Iterate over the alignments
            for alignment in blast_record.alignments:
                # Get the name of the subject sequence
                subject = alignment.hit_def

                genome, gene, exon = ExonAnnotation.parse_exon_subject(subject)
                if gene is not None and exon is not None:
                    for hsp in alignment.hsps:
                        start = hsp.query_start
                        end = hsp.query_end
                        exon_alignment = ExonAlignment(start, end, exon)
                        exon_position_dict[genome].append(exon_alignment)

        return list(exon_position_dict.values())


def cluster_exons(exon_annotations: List[ExonAnnotation], cluster_map, default_cluster):
    for annotation in exon_annotations:
        if annotation.subject not in cluster_map:
            # the exon exist, but didn't get clustered because it's too small/ faulty -> set to default cluster:
            annotation.cluster = default_cluster
        else:
            annotation.cluster = cluster_map[annotation.subject]


def get_exons_from_gene(gene: Gene, genome: AnnotatedGenome, email:str) -> List[ExonAnnotation]:

    try:
        exon_annotations = list()
        for mrna_record in genome.gff_db.children(genome.gff_db[gene.gene_name], featuretype='mRNA'):
            mrna = MRNA.get(mrna_record.id)
            if mrna is None:
                mrna = MRNA(mrna_record.id, mrna_record.start, mrna_record.end, email, genome)
                mrna.gene = gene
                mrna.genome = genome
            for exon_record in genome.gff_db.children(mrna_record, featuretype=FEATURE_TYPE):
                cds = CDS.get_cds_from_exon_name(exon_record.id)
                cds.mrna = mrna
                annotation = ExonAnnotation.get(ExonAnnotation.create_exon_subject(cds, exon_record.id))
                if annotation is None:
                    exon_annotations.append(ExonAnnotation(exon_record.start, exon_record.end, exon_name=exon_record.id, cds=cds))
                else:
                    exon_annotations.append(annotation)
    except FeatureNotFoundError:
        exon_annotations = list()
        logger = logging.getLogger(__name__)
        logger.info(f"gene {gene.gene_name} not found in {genome.name}")
    return exon_annotations


def get_exons_from_gene_group(gene_group: List[Gene], genome, email) -> List[ExonAnnotation]:
    exons = []
    for gene in gene_group:
        exons.extend(get_exons_from_gene(gene, genome, email))
    return exons


def refine_exon_boundaries(chromosome, chromosome_name, exon_boundaries):
    """
    todo: ~~~CURRENTLY WIP~~ (code works in principle, but the project was refactored since, so it needs to be adjusted)
    Refines the boundaries of an exon to comply with the AGGT rule.

    Parameters:
        chromosome (str): The name of the chromosome.
        exon_boundaries (tuple): A tuple with the original start and end coordinates of the exon.
    Returns:
        tuple: A tuple with the refined start and end coordinates of the exon.
    """
    # Load the FASTA file into a Fasta object
    dna_sequence = Fasta(chromosome)[chromosome_name]

    # Extract the DNA sequence for the region represented by the start and end coordinates
    start, end = exon_boundaries

    # Find the starting and ending points for the exon inside the boundary
    # Find the starting and ending points for the exon inside the boundary
    ag_start_inside = -1
    gt_end_inside = -1
    for i in range(start, end):
        if dna_sequence[i:i + 2] == 'AG':
            ag_start_inside = i + 2
            break
    for i in range(end - 1, start - 1, -1):
        if dna_sequence[i:i + 2] == 'GT':
            gt_end_inside = i
            break

    # Find the starting and ending points for the exon outside the boundary
    ag_start_outside = -1
    gt_end_outside = -1
    for i in range(start - 1, -1, -1):
        if dna_sequence[i:i + 2] == 'AG':
            ag_start_outside = i + 2
            break
    for i in range(end, len(dna_sequence)):
        if dna_sequence[i:i + 2] == 'GT':
            gt_end_outside = i
            break

    # Compare the starting and ending points for the exon inside and outside the boundary,
    # and choose the ones that are closest to the original boundary
    if abs(ag_start_inside - start) < abs(ag_start_outside - start):
        best_start = ag_start_inside
    else:
        best_start = ag_start_outside

    if abs(gt_end_inside - end) < abs(gt_end_outside - end):
        best_end = gt_end_inside
    else:
        best_end = gt_end_outside

    # Return the refined start and end coordinates of the exon
    return (best_start, best_end)


def visualize_genes(gene_group: Iterable[Gene]):
    for gene in gene_group:
        exons_in_gene: List[List[ExonAnnotation]] = list()
        for mrna in gene.mrnas:
            for cds in mrna.cdss:
                exons_in_gene.append(cds.exons)

        visualize_sets(align_sets(exons_in_gene))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", help="Email address to make queries to to ncbi's eutilities", required=True)
    args = parser.parse_args()
    email = args.email
    data_path = Path(__file__).parents[1] / Path("data")
    interim_path = data_path / Path("interim")

    human = AnnotatedGenome(name="human", gff_db_file=data_path / "processed/gff_dbs/human.db", chromosome_name='NC_000001.11')
    gorilla = AnnotatedGenome(name="gorilla", gff_db_file=data_path / "processed/gff_dbs/gorilla.db", chromosome_name='NC_044602.1')
    chimpanzee = AnnotatedGenome(name="chimpanzee", gff_db_file=data_path / "processed/gff_dbs/chimpanzee.db", chromosome_name='NC_036879.1')
    dog = AnnotatedGenome(name="dog", gff_db_file=data_path / "processed/gff_dbs/dog.db", chromosome_name='NC_051806.1')
    cat = AnnotatedGenome(name="cat", gff_db_file=data_path / "processed/gff_dbs/cat.db", chromosome_name='NC_058375.1')
    hedgehog = AnnotatedGenome(name="hedgehog", gff_db_file=data_path / "processed/gff_dbs/hedgehog.db", chromosome_name='NW_006804681.1')

    gene_group = [Gene("gene-PLA2G2E", ">", 0),
                  Gene("gene-PLA2G2D", "+", 1),
                  Gene("gene-PLA2G2C", "*", 2),
                  Gene("gene-PLA2G2A", "s", 3),
                  Gene("gene-PLA2G2F", "x", 4)]

    ref_genomes = [human, gorilla, chimpanzee, dog]  # genomes whose exons will be aligned to target genome
    target_genome = cat  # the genome that we are trying to annotate

    # process all relevant exons (that are part of the gene group) from the reference genomes, based on the gff files
    exons: List[ExonAnnotation] = list()
    for genome in ref_genomes:
        exons.extend(get_exons_from_gene_group(gene_group, genome, email))

    # Create exon fasta, to prepare for making a blast database
    exon_records = [SeqRecord(exon.lazy_seq.get_sequence(), id=exon.subject, description="") for exon in exons]
    exon_fasta_path = interim_path / Path(f"exon_fasta/exons.fasta")
    with open(exon_fasta_path, "w") as f:
        SeqIO.write(exon_records, f, "fasta")

    # cluster exons
    default_cluster = ExonCluster(id_=-1)  # some exons are not clustered, associate a default cluster to them
    default_cluster.color = 'grey'
    cluster_map = get_cluster_map(exon_fasta_path, interim_path / Path("exon_clusters/exon_clusters"))
    colorize_clusters(list(set(cluster_map.values())))
    cluster_exons(exons, cluster_map, default_cluster)

    # visualize how the exons are annotated in reference genomes:
    visualize_genes(gene_group)

    # annotated regions in target:
    exons_target = get_exons_from_gene_group(gene_group, target_genome, email)
    cluster_exons(exons_target, cluster_map=cluster_map, default_cluster=default_cluster)

    # translate cds into amino acid encodings (currently doesn't take inverse complement into consideration)
    for gene in gene_group:
        with open(interim_path / Path(f"exon_fasta/cds_aa_{gene.gene_name}.fasta"), "w") as f:
            records = []
            for mrna in gene.mrnas:
                for cds in mrna.cdss:
                    records.append(SeqRecord(cds.aa_seq(), id=cds.id, description=f"{cds.genome.name}"))
            SeqIO.write(records, f, "fasta")

    # define tblastx search region - tblastx is expensive so we want to limit the search to a specific region
    # currently we are "cheating" by looking at the region suggested by the existing annotation. In the future we
    # might look at flanking genes, blastn results (to estimate the region - and then run tblsatx), or both.
    BLAST_SEARCH_MARGIN = 10000
    blast_start = min(exons_target, key=lambda x: x.start).start - BLAST_SEARCH_MARGIN
    blast_end = max(exons_target, key=lambda x: x.end).end + BLAST_SEARCH_MARGIN

    # prepare target sequence (load from ncbi)
    blast_seq = LazyGenomeSeq(email, start=blast_start, end=blast_end, genome=target_genome.chromosome_name)
    target_fasta_file = interim_path / Path(f"fasta_files/{target_genome.chromosome_name}.fasta")
    with open(target_fasta_file, 'w') as f:
        SeqIO.write([SeqRecord(blast_seq.get_sequence(), id=f"{target_genome.chromosome_name}")], f, "fasta")

    # BLAST!
    blast_results = interim_path / Path(f"blast_results/{target_genome.name}_tblastx.xml")
    tblastx(query=target_fasta_file, db=exon_fasta_path, out=blast_results, evalue=0.02)

    # process blast results
    exon_positions = get_exon_positions(blast_results)

    # group exons in target sequence by cds
    grouped_exons_target = [list(group) for key, group in groupby(exons_target, key=lambda exon: exon.cds)]

    # visualize blast results
    # first, visualize each region separately, then together
    visualize_blast_results(exon_positions, grouped_exons_target, separate_plots=True, compare_offset=blast_start)
    visualize_blast_results(exon_positions, grouped_exons_target, separate_plots=False, compare_offset=blast_start)

    # todo: WORK IN PROGRESS
    """
    print("Region Proposals")

    for region_proposal in matched_gene_regions:
        print(f"Region proposal: {region_proposal}")
        refined_region = refine_exon_boundaries(target_genome.fasta_file, target_genome.chromosome_name,
                                                region_proposal)
        print(f"Refined proposal according to AGGT rule: {refined_region}")
    
    """


if __name__ == "__main__":
    main()

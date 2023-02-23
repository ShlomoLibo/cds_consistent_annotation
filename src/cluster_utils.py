import warnings
from os import PathLike
from typing import List, Dict
import matplotlib
import networkx as nx
from Bio.Blast import NCBIXML
from networkx.algorithms.community import greedy_modularity_communities
import os
from annotations import ExonCluster
from external_tools import cd_hit_est, tblastx


def parse_cdhit_clusters_to_dict(cdhit_clusters_file: PathLike) -> tuple[dict[str, int], dict[int, str]]:
    """
    make a dict[subject]->cluster number, out of cd-hit .clstr file and a dict from each cluster to it's centroid
    """
    subjects_to_clusters = {}
    centroids = {}
    with open(cdhit_clusters_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                cluster = int(line.split()[1])
            else:
                subject = line.split(', >')[1].split('...')[0]
                subjects_to_clusters[subject] = cluster
                if line.strip().endswith("*"):
                    centroids[cluster] = subject
    return subjects_to_clusters, centroids


def get_cluster_map(input_fasta: PathLike, output_file: PathLike, sequence_identity=0.95, word_size=10) -> Dict[str, ExonCluster]:
    """
    Uses CD-HIT to cluster exon fasta, see docs: https://github.com/weizhongli/cdhit/wiki/
    Then uses all against all tblastx to further cluster by finding the connected componenets
    Using only CD-HIT doesn't work, too many cluster are created (even for lower sequence_identity)
    Using only tblastx is not efficient
    :param word_size: word_size used for cd-hit
    :param sequence_identity: sequence_identity used for cd-hit
    :param input_fasta: exon fasta file path
    :param output_file: output file path (used file will be {output_file}.clstr)
    :return: dictionary mapping exon annotation headings (matching the blast file) to the their cluster.
    """

    # Cluster CD-HIT
    cd_hit_est(input_fasta, output_file, sequence_identity, word_size)

    # Cluster representatives with BLAST:
    blast_file = os.path.join(os.path.dirname(output_file), "clusters_tblastx.xml")
    tblastx(query=output_file, db=output_file, out=blast_file, evalue=0.02)
    blast_graph = nx.Graph()
    subjects_to_clusters = {}

    parsed_cd_hit, centroids = parse_cdhit_clusters_to_dict(f"{output_file}.clstr")

    for _, centroid in centroids.items():
        blast_graph.add_node(centroid)

    with open(blast_file, "r") as fh:
        # Parse the file
        blast_records = NCBIXML.parse(fh)

        # Iterate over the Blast records
        for blast_record in blast_records:
            # Iterate over the alignments
            for alignment in blast_record.alignments:

                blast_graph.add_edge(blast_record.query, alignment.hit_def)
    communities = greedy_modularity_communities(blast_graph)
    for i, cc in enumerate(communities):
        cluster = ExonCluster(i)
        for subject in cc:
            subjects_to_clusters[subject] = cluster

    for subject, cd_hit_cluster in parsed_cd_hit.items():
        centroid_cluster = subjects_to_clusters[centroids[cd_hit_cluster]]
        subjects_to_clusters[subject] = centroid_cluster

    return subjects_to_clusters


def colorize_clusters(clusters: List[ExonCluster]):
    # set the color attribute of ExonCluster that will be used for visualization:
    if len(clusters) > 60:
        warnings.warn("Colorizing more than 60 clusters is currently not supported")
    cmaps=[matplotlib.cm.get_cmap('tab20'), matplotlib.cm.get_cmap('tab20b'), matplotlib.cm.get_cmap('tab20c')]
    for i, cluster in enumerate(clusters):
        cluster.color = cmaps[i//20](i / 20.0)
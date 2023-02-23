# -*- coding: utf-8 -*-
import logging
from pathlib import Path
import gffutils
# from dotenv import find_dotenv, load_dotenv


def create_db(gff_file: Path, db_name: Path):
    """
    makes a sqlite3 database out of the gff_file
    """

    def cds_key(feature: gffutils.Feature):
        return f"autoincrement:{feature.attributes['ID'][0]}"

    id_spec = {"exon": "ID", "gene": "ID", "mRNA": "ID", "CDS": [cds_key]}
    if not db_name.exists():
        gffutils.create_db(str(gff_file), str(db_name), id_spec=id_spec)


def main():
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making annotation databases from .gff files')
    project_dir = Path(__file__).resolve().parents[2]
    gff_files_dir = project_dir / Path("data/raw/gff_files")
    db_dir = project_dir / Path("data/processed/gff_dbs")

    for gff_path in gff_files_dir.iterdir():
        if gff_path.suffix == ".gff":
            db_name = gff_path.with_suffix(".db").name
            create_db(gff_path, db_dir / db_name)


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    # load_dotenv(find_dotenv())

    main()

import argparse
import logging
import sys
from os import path

# Add path Scripts to sys.path
sys.path.append(path.abspath(path.join(path.dirname(__file__), 'Scripts')))

from Scripts.trimmomatic_run import trimmomatic_run
from Scripts.SPAdes_run import SPAdes_run
from Scripts.bowtie_run import bowtie_run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('REFERENCE', type=str, help='Reference file for alignment')
    
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    reference = args.REFERENCE
    
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s - %(name)s - %(pathname)s:%(lineno)d',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[
                            logging.FileHandler(
                                f'{PROJECT_NAME}_pipeline_multiline.log',
                                mode="w"),
                            logging.StreamHandler()
                        ])
    logger = logging.getLogger(__name__)

    # Execute trimmomatic process
    trimmomatic_run(PROJECT_NAME)

    # # Execute SPAdes analisis
    # SPAdes_run(PROJECT_NAME)

    # # Execute bowtie analisis
    # bowtie_run(PROJECT_NAME, reference)

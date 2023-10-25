import argparse
import logging
import sys
from os import path, makedirs

# Add path Scripts to sys.path
sys.path.append(path.abspath(path.join(path.dirname(__file__), 'Scripts')))

from Scripts.modules.general_functions import read_config
from Scripts.trimmomatic_run import trimmomatic_run
from Scripts.SPAdes_run import SPAdes_run
from Scripts.bowtie_run import bowtie_run
from Scripts.resfinder_run import resfinder_run
from Scripts.oprD_run import oprD_run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('REFERENCE', type=str, help='Reference file for alignment')
    args = parser.parse_args()

    general_config = path.join("configs","general.json")
    config_general = read_config(general_config)
    PROJECTS_PATH = config_general["PROJECTS_PATH"]
    
    PROJECT_NAME = args.PROJECT_NAME
    reference = args.REFERENCE
    
    makedirs(path.join(PROJECTS_PATH, "Logs"), exist_ok=True)

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s - %(name)s - %(pathname)s:%(lineno)d',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[
                            logging.FileHandler(
                                path.join(PROJECTS_PATH, "logs", f'{PROJECT_NAME}_pipeline_multiline.log'),
                                mode="w"),
                            logging.StreamHandler()
                        ])
    logger = logging.getLogger(__name__)

    # Execute trimmomatic process
    trimmomatic_run(PROJECT_NAME)

    # # Execute SPAdes analisis
    SPAdes_run(PROJECT_NAME)

    # # Execute bowtie analisis
    bowtie_run(PROJECT_NAME, reference)

    # # Execute resfinder analisis
    resfinder_run(PROJECT_NAME)

    # # Execute oprD analisis
    oprD_run(PROJECT_NAME)

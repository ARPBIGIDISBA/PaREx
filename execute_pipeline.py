import argparse
import logging
import sys
import os
import glob

# Add path Scripts to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'programs_scripts')))

from programs_scripts.modules.general_functions import read_config, check_project
from programs_scripts.trimmomatic_run import trimmomatic_run
from programs_scripts.SPAdes_run import SPAdes_run
from programs_scripts.bowtie_run import bowtie_run
from programs_scripts.resfinder_run import resfinder_run
from programs_scripts.oprD_run import oprD_run



logger = logging.getLogger(__name__)

if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('operation', type=str, help='Different operations')
    args = parser.parse_args()

    PROJECT_NAME = args.PROJECT_NAME
    OPERATION = args.operation
    
    general_config = os.path.join("configs","general.json")
    config_general = read_config(general_config)
    PROJECTS_PATH = config_general["PROJECTS_PATH"]
    project_path = os.path.join(PROJECTS_PATH, PROJECT_NAME)
    
    check_project(project_path)
    
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s - %(name)s - %(pathname)s:%(lineno)d',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(
                            os.path.join(project_path, "Logs", f'execute_pipeline.log'),
                            mode="w"),
                        logging.StreamHandler()
                    ])

    
    if OPERATION == "create_project":
        # comand line question if you are sure to continue
        logger.info(f"Creating project {PROJECT_NAME} estructure")
        project_path = os.path.join(PROJECTS_PATH, PROJECT_NAME)
        os.makedirs(project_path, exist_ok=True)
        os.makedirs(os.path.join(project_path, f"FASTQ_{PROJECT_NAME}"), exist_ok=True)
        os.makedirs(os.path.join(project_path, f"ANALYSIS_{PROJECT_NAME}"), exist_ok=True)
        with open(os.path.join(project_path, f"SAMPLES_LIST_{PROJECT_NAME}"), 'w') as file:
            pass 
    if OPERATION == "create_sample_list":
        file_name = f"SAMPLES_LIST_{PROJECT_NAME}"
        path = os.path.join(project_path, file_name)
        logger.info(f"Creating {file_name} file")
        logger.debug(f"File path: {path}")
        with open(path, 'w') as file:
            path_fasta = os.path.join(project_path, f"FASTQ_{PROJECT_NAME}")
            fastaq_files = glob.glob(f'{path_fasta}/*.fastq')
            fastaq_gz_files = glob.glob(f'{path_fasta}/*.fastq.gz')
            fastaq_files = fastaq_files + fastaq_gz_files
            if len(fastaq_files) == 0:
                logger.warning(f"No files found in {path_fasta} Using denovo SPAdes path now")
                path_fasta = os.path.join(project_path, f"ANALYSIS_{PROJECT_NAME}", "denovo_assemblies_SPAdes")
                logger.debug(f'{path_fasta}/*.fasta')
                fastaq_files = glob.glob(f'{path_fasta}/*.fasta')
                fastaq_files = [os.path.basename(f).replace(".SPAdes.denovoassembly.fasta", "") for f in fastaq_files]

            else:
                fastaq_files = [os.path.basename(f) for f in fastaq_files]
                fastaq_files = [f.split("_R1")[0] for f in fastaq_files]
                fastaq_files = [f.split("_R2")[0] for f in fastaq_files]
                fastaq_files = list(set(fastaq_files))
            # Print out all filenames
            for filename in fastaq_files:
                print(filename)
                file.write(filename + '\n')
    else:
        logger.warning("Operation not found %s", OPERATION)
        logger.info("Operations available: create_project")
    
    
    
    # reference = args.REFERENCE
    
    # makedirs(path.join(PROJECTS_PATH, "Logs"), exist_ok=True)


    # # Execute trimmomatic process
    # trimmomatic_run(PROJECT_NAME)

    # # # Execute SPAdes analisis
    # SPAdes_run(PROJECT_NAME)

    # # # Execute bowtie analisis
    # bowtie_run(PROJECT_NAME, reference)

    # # # Execute resfinder analisis
    # resfinder_run(PROJECT_NAME)

    # # # Execute oprD analisis
    # oprD_run(PROJECT_NAME)

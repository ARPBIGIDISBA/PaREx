'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
    Here are the command options for spades https://github.com/ablab/spades#sec3.2

'''
import os
import sys
import shutil
import argparse
import logging
from modules.general_functions import read_args, execute_command


logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
default_config_json = os.path.join(script_directory, "SPAdes_config.json")


def SPAdes_run(project_name, config_file=default_config_json):
    # Read command line arguments, sample list and config file
    samples, config = read_args(project_name, config_file)

    PROJECTS_PATH = config["PROJECTS_PATH"]
    # list of coma separated options https://github.com/ablab/spades#sec3.2
    SPADES_PROGRAM_PATH = config['SPADES_PATH']
    SPADES_OPTIONS = config['SPADES_OPTIONS']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    TRIMMOMATIC_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "FASTQ_Trimmomatic")
    if not os.path.exists(TRIMMOMATIC_FILES_PATH):
        logger.error("You have to run first the trimmomatic process")
        logger.error("This forlder does not exist:%s", TRIMMOMATIC_FILES_PATH)
        sys.exit(1)

    OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        sample_name = sample_name.strip()
        logger.info("Processing %s", sample_name)

        # Definir los ficheros de entrada 1 y 2 
        input_r1_path = os.path.join(TRIMMOMATIC_FILES_PATH, f"{sample_name}_trim_R1.fastq")
        input_r2_path = os.path.join(TRIMMOMATIC_FILES_PATH, f"{sample_name}_trim_R2.fastq")
        execute = True
        if not os.path.exists(input_r1_path):
            execute = False
            logger.error("You have to run first the trimmomatic process")
            logger.error("This file does not exist: %s", input_r1_path)
        
        if not os.path.exists(input_r2_path):
            execute = False
            logger.error("You have to run first the trimmomatic process")
            logger.error("This file does not exist: %s", input_r2_path)
        
        if execute:
            command = ["python3", SPADES_PROGRAM_PATH, "-o", OUTPUT_PATH, "-1", input_r1_path, "-2", input_r2_path] + SPADES_OPTIONS

            result = execute_command(command)

            if result:
                logger.info("SPAdes assembly finished")
                # Renombrar los archivos de salida
                old_file_path = os.path.join(OUTPUT_PATH, "contigs.fasta")
                new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

                shutil.move(old_file_path, new_file_path)
                logger.info(f"Rename files for other analysis {new_file_path}")
            else:
                logger.error("SPAdes assembly failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME

    # Start the python logging variable to generate a file
    LOG_MODE = "w"  # "a" to append or "w" to overwrite
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[
                            logging.FileHandler(f'{PROJECT_NAME}_trimmomatic.log', mode=LOG_MODE),
                            logging.StreamHandler()
                        ])
    logger = logging.getLogger(__name__)
    SPAdes_run(PROJECT_NAME, logger)

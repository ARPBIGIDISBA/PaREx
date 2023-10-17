'''
    Este script aplica el progama bowtie para buscar una referencia de alineamiento
    utilizar como entrada dos parametros PROJECT NAME  y fichero de referencia
    https://bowtie-bio.sourceforge.net/tutorial.shtml
'''
import os
import sys
import argparse
import logging
from modules.general_functions import read_args, execute_command, init_configs
from modules.general_functions import configure_logs


logger = logging.getLogger(__name__)

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "bowtie.json")


def bowtie_run(project_name, reference, config=config):
    ''' 
        this function is used to apply the bowtie program to the fastq.gz files
        with a reference file

        parameters:
            project_name (str): Name of the project
            reference (str): Name of the reference file
            config dict with the configuration parameters merge of general.json and bowtie.json

        results:
            this generates a sam file in OUTPUT_PATH
            PROJECTS_PATH/project_name/ANALYSIS_{project_name}/sam_files
    '''

    # Read command line arguments, sample list and config file
    samples = read_args(project_name, config)

    
    # list of coma separated options https://bowtie-bio.sourceforge.net/tutorial.shtml
    BOWTIE_PROGRAM = config['BOWTIE_PATH']
    BOWTIE_OPTIONS = config['BOWTIE_OPTIONS']
    REFERENCE_PATH = config['REFERENCE_PATH']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(config["PROJECTS_PATH"], project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    base_ref = os.path.join(REFERENCE_PATH, reference)
    reference_file = os.path.join(base_ref, reference)
    
    if not os.path.exists(base_ref):
        logger.error("The reference path is not correct")
        logger.error(f"This forlder does not exist:{reference_file}")
        sys.exit(1)

    OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "sam_files")
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    
    TRIMMOMATIC_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "FASTQ_Trimmomatic")
    if not os.path.exists(TRIMMOMATIC_FILES_PATH):
        logger.error("You have to run first the trimmomatic process")
        logger.error("This forlder does not exist:%s", TRIMMOMATIC_FILES_PATH)
        sys.exit(1)

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        sample_name = sample_name.strip()
        logger.info("Processing %s", sample_name)

        # Definir los ficheros de entrada 1 y 2 
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
            output_file = os.path.join(OUTPUT_PATH, f"{sample_name}.sam")
            command = [BOWTIE_PROGRAM, "-x", reference_file, "-q",
                       "-1", input_r1_path,
                       "-2", input_r2_path,
                       "-S", output_file] + BOWTIE_OPTIONS

            result = execute_command(command)

            if result:
                logger.info("Bowtie alignment finished")
            else:
                logger.error("Bowtie alignment failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('REFERENCE', type=str, help='Reference for alignment')
    
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    REFERENCE = args.REFERENCE

    configure_logs(PROJECT_NAME, f"bowtie_{REFERENCE}", config)
    logger = logging.getLogger(__name__)

    bowtie_run(args.PROJECT_NAME, args.REFERENCE, config)

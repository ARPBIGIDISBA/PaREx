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
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "SPAdes.json")


def SPAdes_run(project_name, config=config):
    ''' 
        this function is used to apply the SPAdes program to the fastq.gz files

        parameters:
            project_name (str): Name of the project
            config_json (str): Path to the config file by default is SPAdes_config.json

        results:
            One fasta file with the assembly SPAdes.denovoassembly.fasta written in
            PROJECTS_PATH/project_name/ANALYSIS_{project_name}/denovo_assemblies_SPAdes

    '''

    # Read command line arguments, sample list and config file
    samples = read_args(project_name, config)

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
            old_file_path = os.path.join(OUTPUT_PATH, "contigs.fasta")
            new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")
            result = False

            if not os.path.exists(new_file_path):
                command = ["python3", SPADES_PROGRAM_PATH, "-o", OUTPUT_PATH, "-1", input_r1_path, "-2", input_r2_path] + SPADES_OPTIONS
                result = execute_command(command)

                if result:
                    logger.info("SPAdes assembly finished")
                    # Renombrar los archivos de salida
<<<<<<< HEAD
                    old_file_path = os.path.join(OUTPUT_PATH, "contigs.fasta")
=======
>>>>>>> c1e40a62bc10d3a705aad5ff96bbed3286469bfb
                    new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

                    shutil.move(old_file_path, new_file_path)
                    logger.info(f"Rename files for other analysis {new_file_path}")
                else:
                    logger.error("SPAdes assembly failed")
<<<<<<< HEAD

=======
            else:
                logger.warning(f"File {new_file_path} already exists, skipping")
>>>>>>> c1e40a62bc10d3a705aad5ff96bbed3286469bfb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, args.json_config)

    # Start the python logging variable to generate a file
    configure_logs(project_name, "SPAdes", config)

    logger = logging.getLogger(__name__)

    SPAdes_run(project_name, config)

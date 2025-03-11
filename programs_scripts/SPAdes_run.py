"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the SPAdes program to assemble the reads of the samples
"""
import os
import shutil
import argparse
import logging
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "SPAdes.json", required_keys=["SPADES_PATH", "SPADES_OPTIONS"])


def SPAdes_run(project_name, config=config, extra_config={"force": False, "keep_output": False}):
    ''' 
        this function is used to apply the SPAdes program to the fastq.gz files

        parameters:
            project_name (str): Name of the project
            config_json (str): Path to the config file by default is SPAdes_config.json

        results:
            One fasta file with the assembly SPAdes.denovoassembly.fasta written in
            PROJECTS_PATH/project_name/ANALYSIS_{project_name}/denovo_assemblies_SPAdes

    '''

    direct_file = None
    if extra_config["file"] is not None:
        direct_file = extra_config["file"]
    
    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]

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
        if not os.path.exists(input_r1_path) or not os.path.exists(input_r2_path):
            logger.warning("You have to run first the trimmomatic process")
            if not os.path.exists(input_r1_path):
                execute = False
                logger.warning("This file does not exist: %s", input_r1_path)
            
            if not os.path.exists(input_r2_path):
                execute = False
                logger.warning("This file does not exist: %s", input_r2_path)
        
        if not execute:
            logger.warning("ussing the untrimmed files")
            # Crear los paths de entrada y salida
            input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq.gz")
            input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq.gz")
            logger.warning("Files used: %s %s", input_r1_path, input_r2_path)
            if os.path.exists(input_r1_path) and os.path.exists(input_r2_path):
                execute = True      
            else:
                logger.error("The original FASTQ files does not exist either")
        
        if execute:
            old_file_path = os.path.join(OUTPUT_PATH, "contigs.fasta")
            new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")
            result = False

            if not os.path.exists(new_file_path) or extra_config["force"]:
                command = ["python3", SPADES_PROGRAM_PATH, "-o", OUTPUT_PATH, "-1", input_r1_path, "-2", input_r2_path] + SPADES_OPTIONS
                result = execute_command(command)

                if result:
                    logger.info("SPAdes assembly finished")
                    # Renombrar los archivos de salida
                    new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

                    shutil.move(old_file_path, new_file_path)
                    logger.info(f"Rename files for other analysis {new_file_path}")
                else:
                    logger.error("SPAdes assembly failed")
            else:
                logger.warning(f"File {new_file_path} already exists, skipping")
            
    # Clear output files
    if not extra_config["keep_output"]:
        for file_ in os.listdir(OUTPUT_PATH):
            # remove all the files and directory minus .log and SPAdes.denovoassembly.fasta
            if not file_.endswith("SPAdes.denovoassembly.fasta") and not file_.endswith(".log"):
                file_path = os.path.join(OUTPUT_PATH, file_)
                if os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                else:
                    os.remove(file_path)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default="INFO")
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Keep the output files')

    
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, args.json_config, required_keys=["SPADES_PATH", "SPADES_OPTIONS"])

    # Start the python logging variable to generate a file
    configure_logs(project_name, "SPAdes", config)

    logger = logging.getLogger(__name__)

    SPAdes_run(project_name, config, extra_config={"force": args.force, "keep_output": args.keep_output})

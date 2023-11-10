'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
    Here are the command options for spades https://github.com/ablab/spades#sec3.2

'''
import os
import sys
import argparse
import logging
import json
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "mlst.json")

def mlst_run(project_name, config=config, direct_file = None):
    ''' 
        this function is used to apply the resfinder program to the denovo files output of SPAdes

        parameters:
            project_name (str): Name of the project
            config dict (dict): readed from Path  is mlst.json
            direct_file (str): path to the file instead of list in sample list

        results:

    '''

    # Read command line arguments, sample list and config file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]

    PROJECTS_PATH = config["PROJECTS_PATH"]
    # list of coma separated options https://github.com/ablab/spades#sec3.2
    
    MLST_OPTIONS = config['MLST_OPTIONS']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "mlst_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    
    results_data = [
        ["sample_name", "scheme", "sequence_type", "alleles"]
    ]
    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        if direct_file:
            SPADES_FILE = sample_name
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
        else:
            sample_name = sample_name.strip()
            logger.info("Processing sample %s", sample_name)
            # Definir los ficheros de entrada 1 y 2 
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

        logger.info("Using SPAdes file: %s", SPADES_FILE)
        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the trimmomatic / SPades process or use a difect file --file path_to_file")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            os.makedirs(os.path.join(OUTPUT_PATH, "outputs"), exist_ok=True)
            json_file = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}.json")
            command_mslt = [config["MLST_PATH"], SPADES_FILE] + MLST_OPTIONS + ["--json", json_file, "--label", sample_name]
            result = execute_command(command_mslt)
            if result:
                if os.path.exists(json_file):
                    with open(json_file) as json_file:
                        data = json.load(json_file)
                        record = data[0]
                        id_ = record["id"]
                        scheme = record["scheme"]
                        sequence_type = record["sequence_type"]
                        alleles = record["alleles"]
                        output = f"{id_} {scheme} {sequence_type} "
                        logger.info("MLST output: %s", output)
                        # Sort alleles alphabetically by gene name
                        sorted_alleles = sorted(alleles.items(), key=lambda x: x[0])
                        output += " ".join([f"{gene}({allele})" for gene, allele in sorted_alleles])
                        alleles = " ".join([f"{gene}({allele})" for gene, allele in sorted_alleles])
                        results_data.append([sample_name, scheme, sequence_type, alleles])
                        logger.info(output)

    # Crear y escribir en el archivo CSV usando punto y coma como delimitadorç
    filename = os.path.join(OUTPUT_PATH, f"{project_name}_mlst_results.csv")
    with open(filename, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file, delimiter=';')
        writer.writerows(results_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "oprD", config)

    logger = logging.getLogger(__name__)

    mlst_run(project_name, config, args.file)

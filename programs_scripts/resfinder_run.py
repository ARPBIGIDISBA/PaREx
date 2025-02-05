'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
    Here are the command options for spades https://github.com/ablab/spades#sec3.2

'''
import os
import argparse
import logging
import json
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "resfinder.json", required_keys=["RESFINDER_PATH", "RESFINDER_OPTIONS", "INTRINSIC_PAER_GENES"])


def print_metadata(data):
    logger.info("**********************")
    logger.info("***** Metadata ******")
    logger.info("**********************")
    for db_key, db_info in data["databases"].items():
        logger.info(f"Database used: {db_key}")
        string = ""
        for info_key, info_value in db_info.items():
            string += f"{info_key}: {info_value} "
        logger.info(string)

    logger.info("Software %s version %s", data["software_name"], data["software_version"])
    logger.info("**********************")


def filter_output(data, ignore_list):
    '''
        This function is used to filter the output of the resfinder program 
        and generate csv output clean
    '''
    csv_fullcoverage = "name;identity;coverage;ref_start_pos;ref_end_pos;ref_seq_lenght;coverage;ref_id;query_id;query_start_pos;query_end_pos;ref_acc;grade;phenotypes\n"
    csv_partialcoverage = csv_fullcoverage
    posible = False
    full = False
    for seq_key, seq_info in data["seq_regions"].items():
        name = seq_info["name"]
        if name not in ignore_list:
            alignment = seq_info["alignment_length"]
            seq_length = seq_info["ref_seq_length"]
            identity = seq_info["identity"]
            start_pos = seq_info["ref_start_pos"]
            end_pos = seq_info["ref_end_pos"]
            ref_query = seq_info["ref_seq_length"]
            coverage = seq_info["coverage"]
            phenotypes = ', '.join(seq_info['phenotypes'])
            #logh keys seq_info
            logger.info("Gene: %s identity %2.f. (%s, %s) %s %s", name, identity, start_pos, end_pos, identity, coverage)
            line = f"{name};{identity};{coverage};{start_pos};{end_pos};{ref_query};{coverage};{seq_info['ref_id']};{seq_info['query_id']};"
            line += f"{seq_info['query_start_pos']};{seq_info['query_end_pos']};{seq_info['ref_acc']};{seq_info['grade']};{phenotypes}\n"

            # Change to postive first
            if (alignment == seq_length and coverage >= 100 and identity>=100) or name=="crpP":
                full = True
                csv_fullcoverage = csv_fullcoverage + line
            else:
                if alignment != seq_length:
                    logger.info("   Distint lenght: %s %s", alignment, seq_length)
                if coverage < 100:
                    logger.info("   Coverage minus 100%%: %s", coverage)
                posible = True
                csv_partialcoverage = csv_partialcoverage + line
            

    # To not generate the csv file if there is no result
    if not full:
        full = False
    else:
        csv_fullcoverage = csv_fullcoverage.replace(".", ",")
    if not posible:
        posible = False
    else:
        csv_partialcoverage = csv_partialcoverage.replace(".", ",")
    return csv_fullcoverage, csv_partialcoverage


def resfinder_run(project_name, config=config, only_output=False, direct_file = None, extra_config={"force": False, "keep_output": False}):
    ''' 
        this function is used to apply the resfinder program to the denovo files output of SPAdes

        parameters:
            project_name (str): Name of the project
            config dict (dict): readed from Path  is resfinder.json

        results:
            

    '''

    if direct_file is None and extra_config["file"] is not None:
        direct_file = extra_config["file"]
        
    # Read command line arguments, sample list and config file
    if direct_file:
        samples = [direct_file]
    else:
        samples = read_args(project_name, config)

    PROJECTS_PATH = config["PROJECTS_PATH"]
    # list of coma separated options https://github.com/ablab/spades#sec3.2
    RESFINDER_PROGRAM_PATH = config['RESFINDER_PATH']
    RESFINDER_OPTIONS = config['RESFINDER_OPTIONS']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "resfinder_results" )
    OUTPUT_PATH_SCRIPT = os.path.join(OUTPUT_PATH, "output" )
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    previous_dir = os.getcwd()
    os.chdir(RESFINDER_PROGRAM_PATH)

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        sample_name = sample_name.strip()
        
        if direct_file:
            SPADES_FILE = sample_name
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
        else:
            sample_name = sample_name.strip()
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

        logger.info("Processing sample %s", sample_name)

        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the SPAdes process")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            command = ["python3", os.path.join(RESFINDER_PROGRAM_PATH, "run_resfinder.py"), "-o", OUTPUT_PATH_SCRIPT, "-s", "OTHER", "-ifa", SPADES_FILE] + RESFINDER_OPTIONS
            output_json = os.path.join(OUTPUT_PATH_SCRIPT, f"{sample_name}.json")
            
            if only_output:
                if os.path.exists(output_json):
                    result = True
                else:
                    logger.error("You have to run first the resfinder process")
                    logger.error("File not found: %s", output_json)
            else:
                result = execute_command(command)

            if result:
                # Read the json file and get the results
                with open(output_json) as json_file:
                    data = json.load(json_file)
                    logger.info("Resfinder results for sample %s", sample_name)
                    csv_fullcoverage, csv_partialcoverage = filter_output(data, config["INTRINSIC_PAER_GENES"])
                    os.makedirs(os.path.join(OUTPUT_PATH,"csv_samples"), exist_ok=True)
                    if csv_fullcoverage:
                        with open(os.path.join(OUTPUT_PATH, "csv_samples", f"{sample_name}.fullcoverage.csv"), "w") as file:
                            file.write(csv_fullcoverage)
                    if csv_partialcoverage:
                        with open(os.path.join(OUTPUT_PATH, "csv_samples", f"{sample_name}.partialcoverage.csv"), "w") as file:
                            file.write(csv_partialcoverage)
            else:
                logger.error("Resfinder failed assembly failed on sample %s", sample_name)

    os.chdir(previous_dir)
    if not extra_config["keep_output"]:
        os.system(f"rm -rf {OUTPUT_PATH_SCRIPT}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Keep the output files')

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, args.json_config, required_keys=["RESFINDER_PATH", "RESFINDER_OPTIONS", "INTRINSIC_PAER_GENES"])
        
    # Start the python logging variable to generate a file
    configure_logs(project_name, "resfinder", config)

    resfinder_run(project_name, config, args.parse_output, args.file, extra_config={"force": args.force, "keep_output": args.keep_output})

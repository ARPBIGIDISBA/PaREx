"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the oprD program to align the reads with a reference file
"""
import os
import argparse
import logging
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
from modules.blast_functions import get_differences, analize_sample

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "oprD.json", required_keys=["BLAST_OPTIONS", "BLASTN_OPTIONS"])

def oprD_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False, extra_config= {"force": False, "keep-output": False}):
    ''' 
        this function is used to apply the oprD program to the denovo files output of SPAdes

        parameters:
            project_name (str): Name of the project
            config dict (dict): readed from Path  is resfinder.json
            only_output (bool): Set the flag to not execute but only process json file
            direct_file (str): path to the file instead of list in sample list
            normal_output (bool): Set the flag to not execute but only process normal blast outputs
        results:

    '''

    if direct_file is None and extra_config["file"] is not None:
        direct_file = extra_config["file"]
    
    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]

    PROJECTS_PATH = config["PROJECTS_PATH"]
    BLAST_OPTIONS = config['BLAST_OPTIONS']
    BLASTN_OPTIONS = config['BLASTN_OPTIONS']
    NUCLEOTIDE_PATH = os.path.join(config['DATABASE_PATH'], 'oprD', 'nucleotide')
    PROTEIN_PATH = os.path.join(config['DATABASE_PATH'], 'oprD', 'protein')

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    # Filter only files finishing in .fasta
    files_nucleotide = [f for f in os.listdir(NUCLEOTIDE_PATH) if f.endswith(".fasta")]

    
    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "oprD_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    results_data = [
        ["sample_name", "oprD", "oprD_REFERENCE", "bit_score", "gaps", "identity"]
    ]

    for sample_name in samples:
        logger.info("***********************************************")
        if direct_file:
            SPADES_FILE = sample_name
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
        else:
            sample_name = sample_name.strip()
            logger.debug("Processing sample %s", sample_name)
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

        logger.debug("Using SPAdes file: %s", SPADES_FILE)
        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the SPades process or use a difect file --file path_to_file")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            max_bitscore = {
                "name": "deleted",
                "value": -1,
                "gaps":-1,
                "identity": -1,
                "hsps": -1,
                "differences": ""
            }

            for nucleotide_file in files_nucleotide:
                suffix_name = nucleotide_file[0:-len("_nucleotide.fasta")]
                name = suffix_name.replace("oprD_", "")
                logger.debug("Processing sample %s", name)

                ## nucleotide analysis
                nucleotide_file = os.path.join(NUCLEOTIDE_PATH, nucleotide_file)
                logger.debug("Using nucleotide file: %s", nucleotide_file)
                output_file_nucleotide = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{name}_nucleotide.json")
                os.makedirs(os.path.join(OUTPUT_PATH, "output"), exist_ok=True)

                logger.debug("Output file %s", output_file_nucleotide)
                if normal_output:
                    command_nucleotide = ["blastn", "-query", nucleotide_file, "-subject", SPADES_FILE, "-out", output_file_nucleotide.replace(".json", "")] + BLAST_OPTIONS
                else:
                    command_nucleotide = ["blastn", "-query", nucleotide_file, "-subject", SPADES_FILE, "-out", output_file_nucleotide, "-outfmt", "15"] + BLAST_OPTIONS
                
                if only_output and not normal_output:
                    if os.path.exists(output_file_nucleotide):
                        result = True
                    else:
                        logger.error("File not found: %s", output_file_nucleotide)
                        logger.error("You have to run first the oprD process")
                        result = False
                else:
                    result_nuc = execute_command(command_nucleotide)
                    result = result_nuc
                if result and not normal_output:
                    # Read the json file and get the results
                    logger.debug("***********************************************")
                    logger.debug("oprD Analysis sample %s againts %s", sample_name, name)
                    logger.debug("***********************************************")
                    
                    logger.debug("-----------------------------------------------")
                    logger.debug("--- Nucleotide analysis %s ----------------", sample_name)
                    logger.debug("-----------------------------------------------")
                    results = analize_sample(output_file_nucleotide, name, "nucleotide")
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]

                    logger.debug("Gaps %s", gaps)
                    logger.debug("Bit score %s", bit_score)                   
                    
                    if bit_score > max_bitscore["value"]:
                        logger.debug("------------New max bit score %s", bit_score)
                        max_bitscore["name"] = name
                        max_bitscore["value"] = bit_score
                        max_bitscore["gaps"] = gaps
                        max_bitscore["identity"] = identity
                        max_bitscore["hsps"] = results["hsps"]
                        max_bitscore["differences"] = results["differences"]
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("oprD failed assembly failed on sample %s", sample_name)
            
            logger.info("Max bit score for %s value %s gaps %s Identity %s", max_bitscore["name"], max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"])
            
            if max_bitscore["gaps"] == -1:
                logger.warning("oprD failed assembly failed on sample %s", sample_name)
                logger.warning(results)
                results_data.append([sample_name, results["differences"], max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                logger.info("This is a Wild Type (WT). No gaps in nucleotide")
                logger.info("***********************************************")  
                results_data.append([sample_name, "WT", max_bitscore["name"],  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] < 100:
                logger.info("Analyse the differences in protein")
                logger.debug("***********************************************")        
            
                protein_file = os.path.join(PROTEIN_PATH, f"oprD_{ max_bitscore['name']}_protein.fasta")
                logger.debug("Use protein file %s", protein_file)
                if os.path.exists(protein_file):
                    logger.debug("Using protein file: %s", protein_file)
                    output_file_protein = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{ max_bitscore['name']}_protein.json")
                    os.makedirs(os.path.join(OUTPUT_PATH, "output"), exist_ok=True)
                
                    logger.debug("Output file %s", output_file_protein)
                    if normal_output:
                        command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein.replace(".json", "")] + BLASTN_OPTIONS
                    else:
                        command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein, "-outfmt", "15"] + BLASTN_OPTIONS
                    
                    if only_output and not normal_output:
                        if os.path.exists(output_file_protein):
                            result = True
                        else:
                            logger.error("File not found: %s", output_file_protein)
                            logger.error("You have to run first the oprD process")
                            result = False
                    else:
                        result = execute_command(command_protein)

                    if result and not normal_output:
                        # Read the json file and get the results
                        logger.debug("***********************************************")
                        logger.info("Protein Analysis sample %s againts %s", sample_name, max_bitscore['name'])
                        logger.debug("***********************************************")
                        
                        results = analize_sample(output_file_protein, max_bitscore['name'], "protein")
                        if results["gaps"] == 0 and results["identity"] == 100:
                            logger.info("This is a Wild Type (WT). No gaps in protein")
                            results_data.append([sample_name, "WT", results["name"],  results["bit_score"], results["gaps"], results["identity"]])
                        else:
                            result = [sample_name, ",".join(results["differences"]), max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]]
                            results_data.append(result)
                    else:
                        if not normal_output and not only_output:
                            logger.error("oprD failed assembly failed on sample %s", sample_name)
            else:
                differences = get_differences(max_bitscore["hsps"],  max_bitscore["name"], max_bitscore["gaps"], "nucleotide")
                results_data.append([sample_name, ",".join(differences), max_bitscore["name"],  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])    
                
    # Crear y escribir en el archivo CSV usando punto y coma como delimitador
    filename = os.path.join(OUTPUT_PATH, f"{project_name}_oprD_results.csv")
    logger.debug(filename)
    with open(filename, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file, delimiter=';')
        writer.writerows(results_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--normal-output', action='store_true', help='Produce only screen process normal blast outputs') 
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep-output', action='store_true', help='Force the execution of the program')

    
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "oprD", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    oprD_run(project_name, config, args.parse_output, args.file, args.normal_output, extra_config= {"force": args.force, "keep-output": args.keep_output})

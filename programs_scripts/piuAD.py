"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the oprD program to align the reads with a reference file
"""
import os
import argparse
import logging
import csv
from modules.general_functions import read_args
from modules.general_functions import configure_logs, init_configs
from modules.blast_functions import get_differences, analize_sample, run_blast

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "piuAD.json", required_keys=["BLAST_OPTIONS", "BLASTN_OPTIONS"])

def piuAD_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False, extra_config= {"force": False, "keep-output": False}):
    ''' 
        This function is used to apply a blast with piuA and piuD genes to the denovo files output of SPAdes

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
    NUCLEOTIDE_PATH = os.path.join(config['DATABASE_PATH'], 'piuAD', 'nucleotide')
    PROTEIN_PATH = os.path.join(config['DATABASE_PATH'], 'piuAD', 'protein')

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    file_list = [{
        "fasta_filename": "PA4514_piuA.fasta",
        "name": "piuA",
        "locus_tag": "PA4514",
        "result": "A"
    }, {
        "fasta_filename": "PALES_48941_piuD.fasta",
        "name": "piuD",
        "locus_tag": "PALES_48941",
        "result": "D"
    }]

    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "piuAD_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    results_data = [
        ["sample_name", "piuA/D", "piuA/D_REFERENCE", "gaps", "identity"]
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

            for fast_file_compare in file_list:
                
                logger.debug("Processing sample %s", fast_file_compare)

                ## nucleotide analysis
                file = os.path.join(NUCLEOTIDE_PATH, fast_file_compare["fasta_filename"])
                name = f"{fast_file_compare['name']}_nucleotide"
                logger.debug("Using nucleotide file: %s", file)
                output_file_nucleotide = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}_{name}.json")                
                os.makedirs(os.path.join(OUTPUT_PATH, "outputs"), exist_ok=True)
                
                result = run_blast(sample_name, name, file, OUTPUT_PATH, SPADES_FILE, BLAST_OPTIONS, normal_output, only_output, tblastn=False)
                
                if result and not normal_output:
                    # Read the json file and get the results
                    logger.debug("***********************************************")
                    logger.debug("piuAD Analysis sample %s againts %s", sample_name, name)
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
                        max_bitscore["result"] = fast_file_compare["result"]
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("piuAD failed assembly failed on sample %s", sample_name)
            
            logger.debug("***********************************************")
            logger.info("Final piuAD Analysis sample %s", sample_name)
            logger.debug("***********************************************")
            logger.info("Max bit score %s against %s", max_bitscore["value"], max_bitscore["name"])
            piuAD = "deleted"
            if max_bitscore["name"].startswith("piuD"):
                piuAD = "D"
            elif max_bitscore["name"].startswith("piuA"):
                piuAD = "A"
            
            if max_bitscore["gaps"] == -1:
                logger.warning("piuAD failed assembly failed on sample %s", sample_name)
                logger.warning(results)
                results_data.append([sample_name, piuAD, results["differences"], results["gaps"], results["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                logger.info("This is a Wild Type (WT). No gaps in nucleotide")
                logger.info("***********************************************")  
                results_data.append([sample_name, piuAD, "WT", max_bitscore["gaps"], max_bitscore["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] < 100:
                logger.info("Analyse the differences in protein")
                logger.debug("***********************************************")

                file = os.path.join(PROTEIN_PATH, fast_file_compare["fasta_filename"])
                name = f"{fast_file_compare['name']}_protein"
                logger.debug("Using protein file: %s", file)
                output_file_protein = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}_{name}.json")
                os.makedirs(os.path.join(OUTPUT_PATH, "outputs"), exist_ok=True)
                result = run_blast(sample_name, name, file, OUTPUT_PATH, SPADES_FILE, BLASTN_OPTIONS, normal_output, only_output, tblastn=True) 
          
                if result and not normal_output:
                    # Read the json file and get the results
                    logger.debug("***********************************************")
                    logger.info("Protein Analysis sample %s againts %s", sample_name, max_bitscore['name'])
                    logger.debug("***********************************************")
                    
                    results = analize_sample(output_file_protein, max_bitscore['name'], "protein")

                    if results["gaps"] == 0 and results["identity"] == 100:
                        logger.info("This is a Wild Type (WT). No gaps in protein")
                        results_data.append([sample_name, piuAD, "WT", results["gaps"], results["identity"]])
                    else:
                        result = [sample_name, piuAD, ",".join(results["differences"]), results["gaps"], results["identity"]]
                        results_data.append(result)
                else:
                    if not normal_output and not only_output:
                        logger.error("piuAD failed assembly failed on sample %s", sample_name)
            else:
                differences = get_differences(max_bitscore["hsps"],  max_bitscore["name"], max_bitscore["gaps"], "nucleotide")
                results_data.append([sample_name, piuAD, ",".join(differences), max_bitscore["gaps"], max_bitscore["identity"]])    
                
    # Crear y escribir en el archivo CSV usando punto y coma como delimitador
    filename = os.path.join(OUTPUT_PATH, f"{project_name}_piuAD_results.csv")
    logger.info("Generated filename: %s", filename)
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
    configure_logs(project_name, "piuAD", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    piuAD_run(project_name, config, args.parse_output, args.file, args.normal_output, extra_config= {"force": args.force, "keep-output": args.keep_output})

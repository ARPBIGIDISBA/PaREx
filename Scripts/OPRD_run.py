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
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "OPRD.json")

def print_metadata(json_file):
    protein_data = json.load(open(json_file))['BlastOutput2'][0]['report']
    logger.info("Program: %s version %s",protein_data['program'], protein_data['version'])
    params_str = ', '.join(f"{k}: {v}" for k, v in protein_data["params"].items())
    logger.info("Params used: %s", params_str)

    for key, value in protein_data["results"].items():
        for result in value:
            logger.info("-----------------------------------------------")
            logger.info("Result of %s", result["query_title"])
            logger.info("Result of %s", key)
            logger.info("Number of hits: %s", len(result["hits"]))
            stats_str = ', '.join(f"{k}: {v}" for k, v in result["stat"].items())
            logger.info("Stats: %s", stats_str)
            for hit in result["hits"]:
                logger.info("Hit: %s", hit["description"][0]["title"])

def OPRD_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False):
    ''' 
        this function is used to apply the resfinder program to the denovo files output of SPAdes

        parameters:
            project_name (str): Name of the project
            config dict (dict): readed from Path  is resfinder.json

        results:
            

    '''

    # Read command line arguments, sample list and config file
    if not direct_file:
        samples = read_args(project_name, config)

    PROJECTS_PATH = config["PROJECTS_PATH"]
    # list of coma separated options https://github.com/ablab/spades#sec3.2
    
    BLAST_OPTIONS = config['BLAST_OPTIONS']
    BLASTN_OPTIONS = config['BLASTN_OPTIONS']
    
    NUCLEOTIDE_PATH = config["NUCLEOTIDE_PATH"]
    PROTEIN_PATH = config["PROTEIN_PATH"]

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    files_nucleotide = os.listdir(NUCLEOTIDE_PATH)
    
    
    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "OPRD_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    # Create a CSV file and write the header
    output_csv = os.path.join(OUTPUT_PATH, 'summary.csv')
    logger.info("Generate generic CSV output file: %s", output_csv)
    with open(output_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        headers = ["Sample", "Gene", "Identity", "Coverage", "Start", "End", "Strand", "Evalue", "Accession", "Product"]
        csv_writer.writerow(headers)

    if direct_file:
        samples = [direct_file]

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        if direct_file:
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
            SPADES_FILE = direct_file
        else:
            sample_name = sample_name.strip()
            logger.info("Processing sample %s", sample_name)
            # Definir los ficheros de entrada 1 y 2 
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")
        logger.info("Using SPAdes file: %s", SPADES_FILE)
        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the trimmomatic process")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            samples_tested = []
            for nucleotide_file in files_nucleotide:
                name = nucleotide_file[0:-len("_nucleotide.fasta")]
                samples_tested.append(name)
                logger.info("Processing sample %s", name)

                ## nucleotide analysis
                nucleotide_file = os.path.join(NUCLEOTIDE_PATH, nucleotide_file)
                logger.info("Using nucleotide file: %s", nucleotide_file)
                output_file_nucleotide = os.path.join(OUTPUT_PATH, f"{sample_name}_{name}_nucleotide.json")
                logger.info("Output file %s", output_file_nucleotide)
                if normal_output:
                    command_nucleotide = ["blastn", "-query", nucleotide_file, "-subject", SPADES_FILE, "-out", output_file_nucleotide.replace(".json", "")] + BLAST_OPTIONS
                else:
                    command_nucleotide = ["blastn", "-query", nucleotide_file, "-subject", SPADES_FILE, "-out", output_file_nucleotide, "-outfmt", "15"] + BLAST_OPTIONS
                
                # Protein analysis
                protein_file = os.path.join(PROTEIN_PATH, f"{name}_protein.fasta")
                if os.path.exists(protein_file):
                    logger.info("Using protein file: %s", protein_file)
                    output_file_protein = os.path.join(OUTPUT_PATH, f"{sample_name}_{name}_protein.json")
                    logger.info("Output file %s", output_file_protein)
                    if normal_output:
                        command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein.replace(".json", "")] + BLASTN_OPTIONS
                    else:
                        command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein, "-outfmt", "15"] + BLASTN_OPTIONS
                else:
                    logger.error("File not found: %s", protein_file)
                    logger.error("Every nucleotide file must have a protein file")
                    result = False
                if only_output and not normal_output:
                    if os.path.exists(output_file_nucleotide) and os.path.exists(output_file_protein):
                        result = True
                    else:
                        logger.error("You have to run first the OPRD process")
                        logger.error("File not found: %s", output_file_nucleotide)
                        logger.error("File not found: %s", output_file_protein)
                else:
                    result_nuc = execute_command(command_nucleotide)
                    result_pro = execute_command(command_protein)
                    result = result_nuc and result_pro

                if result and not normal_output:
                    # Read the json file and get the results
                    logger.info("***********************************************")
                    logger.info("OPRD Analysis sample %s againts %s", sample_name, name)
                    logger.info("***********************************************")
                    
                    logger.info("**** Nucleotide analysis ****")
                    print_metadata(output_file_nucleotide)
                    # logger.info("**** Protein analysis ****")
                    # print_metadata(output_file_protein)
                    


                else:
                    logger.error("Resfinder failed assembly failed on sample %s", sample_name)
                break

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--normal-output', action='store_true', help='Produce only screen process normal blast outputs')   

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    # Start the python logging variable to generate a file
    configure_logs(project_name, "OPRD", config)

    logger = logging.getLogger(__name__)

    OPRD_run(project_name, config, args.parse_output, args.file, args.normal_output)

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
from modules.general_functions import get_spades_file
from modules.blast_functions import analize_sample, run_blast 

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "geneabsence.json", required_keys=["BLAST_OPTIONS"])

def gene_absence_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False, extra_config= {"force": False, "keep-output": False}):
    ''' 
        this function is used to apply the resfinder program to the denovo files output of SPAdes

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

    FASTA_PATH = os.path.join(config['DATABASE_PATH'], 'gene_absence', 'fasta_sequences')

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    # Filter only files finishing in .fasta
    
    # Read database
    GENES_STUDY = []
    for f in os.listdir(FASTA_PATH):
        test = f.endswith(".fasta")
        if test:
            name = f[0:-len(".fasta")]
            EXTRA_PATH = os.path.join(FASTA_PATH, name+"_extra")
            output = {
                "name": name,
                "extra": []
            }
            if os.path.exists(EXTRA_PATH):
                for ef in os.listdir(EXTRA_PATH):
                    if ef.endswith(".fasta"):
                        name_extra = ef[0:-len(".fasta")]
                        output["extra"].append(name_extra)
            GENES_STUDY.append(output)
    
    
    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "gene_absence_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    header = ["sample_name"]
    for gene in GENES_STUDY:
        header.append(gene['name'])
        for extra in gene['extra']:
            header.append(extra)
    # add genes of study and
    results_data = [header]

    logger.debug("Samples to process: %s", header)

    for sample_name in samples:
        logger.info("***********************************************")
        logger.info("Processing sample %s", sample_name)
        logger.info("***********************************************")

        SPADES_FILE, sample_name, execute = get_spades_file(sample_name, direct_file=direct_file, SPADES_FILES_PATH=SPADES_FILES_PATH)

        row = [sample_name]
        if execute:
            for gene in GENES_STUDY:
                gene_name = gene["name"]

                logger.info("Analyzing gene %s", gene_name)
                compare_file = os.path.join(FASTA_PATH, gene_name+".fasta")
                output_file = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}_{gene_name}.json")
                result = run_blast(sample_name, gene_name, compare_file, OUTPUT_PATH, SPADES_FILE, BLAST_OPTIONS, normal_output, only_output, tblastn=False)
                                
                if result and not normal_output:
                    # Read the json file and get the results
                    logger.debug("***********************************************")
                    logger.debug(" Gene Absence Analysis sample %s againts %s", sample_name, gene_name)
                    logger.debug("***********************************************")
                    results = analize_sample(output_file, gene_name, "nucleotide")
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]
                    differences = results["differences"]
                    if differences == "deleted":
                        row.append("Deleted")
                        #analize extra files
                        for extra in gene["extra"]:
                            logger.info("Analyzing extra gene %s", extra)
                            compare_file_extra = os.path.join(FASTA_PATH, gene_name+"_extra", extra+".fasta")
                            output_file_extra = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}_{extra}.json")
                            result_extra = run_blast(sample_name, extra, compare_file_extra, OUTPUT_PATH, SPADES_FILE, BLAST_OPTIONS, normal_output, only_output, tblastn=False)
                            if result_extra and not normal_output:
                                logger.debug("***********************************************")
                                logger.debug(" Gene Absence Analysis sample %s againts extra %s", sample_name, extra)
                                logger.debug("***********************************************")
                                results_extra = analize_sample(output_file_extra, extra, "nucleotide")

                                differences = results_extra["differences"]
                                if differences == "deleted":
                                    row.append("Deleted")
                                else:
                                    row.append("")
                            else:
                                row.append("Gene absence analysis failed on sample {}".format(sample_name))
                    else:
                        logger.debug("Gene present Gaps: {}, Bit score: {}, Identity: {:.2f}%".format(gaps, bit_score, identity))
                        row.append("")
                        for extra in gene["extra"]:
                            row.append("")
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("Gene absence analysis failed on sample %s", sample_name)
            
            results_data.append(row)
                
    # Crear y escribir en el archivo CSV usando punto y coma como delimitador
    filename = os.path.join(OUTPUT_PATH, f"{project_name}_gene_absence_results.csv")
    logger.info(f"Generated {filename}")
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
    configure_logs(project_name, "gene_absense", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    gene_absence_run(project_name, config, args.parse_output, args.file, args.normal_output, extra_config= {"force": args.force, "keep-output": args.keep_output})

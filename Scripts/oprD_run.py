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
config = init_configs(script_directory, "oprD.json")

def get_differences_protein(hsps, name, gaps = 0):
        qseq = hsps["qseq"]
        midline = hsps["midline"]
        hseq = hsps["hseq"]
        differences = []
        qstate = False
        hstate = False
        mstate = False
        for i, (q, m, h) in enumerate(zip(qseq, midline, hseq)):
            if q =="-":
                if not qstate:
                    qstate = True
                    if gaps > 0:
                        differences.append(f"nt{name}ins{gaps}")
                    else:
                        differences.append(f"X{i+1}{h}")
            else:
                qstate = False
            if h =="-":
                if not hstate:
                    hstate = True
                    if gaps > 0:
                        differences.append(f"nt{name}del{gaps}")
                    else:
                        differences.append(f"{q}{i+1}X")
            else:
                hstate = False
            
            if m ==" ":
                if not mstate:
                    hstate = True
                    if gaps > 0:
                        differences.append(f"nt{name}del{gaps}")
                    else:
                        differences.append(f"{q}{i+1}X")
            else:
                hstate = False

        for difference in differences:
            logger.info("%s %s",name,difference)
        return differences

def read_output(json_file):
    data = json.load(open(json_file))['BlastOutput2'][0]['report']
    return data

def analize_sample(json_file, name, nucleotide_protein = "nucleotide"):
    protein_data = read_output(json_file)
    
    logger.info("Program: %s version %s",protein_data['program'], protein_data['version'])
    params_str = ', '.join(f"{k}: {v}" for k, v in protein_data["params"].items())
    logger.info("Params used: %s", params_str)

    for key, value in protein_data["results"].items():
        for result in value:
            logger.info("Result of %s key %s", result["query_title"], key)
            if len(result["hits"]) > 1:
                logger.info("Number of hits: %s", len(result["hits"]))
            stats_str = ', '.join(f"{k}: {v}" for k, v in result["stat"].items())
            logger.info("Stats: %s", stats_str)
                
            query_len = result["query_len"]
            
            if nucleotide_protein == "nucleotide":
                if query_len >=1296 and query_len <= 1352:
                    for hit in result["hits"]:
                        title = hit["description"][0]["title"]
                        logger.info("%s Hit: %s", name, title)
                        bit_score = 0
                        for hsps in hit["hsps"]:
                            gaps = hsps["gaps"]
                            identity = hsps["identity"]/hsps["align_len"]*100
                            if hsps["bit_score"] > bit_score:
                                bit_score = hsps["bit_score"]
                    return {"gaps": gaps, "bit_score": bit_score, "identity": identity}
                
            elif nucleotide_protein == "protein":
                if query_len >=441 and query_len <= 443:
                    best_match = {
                        "bit_score": 0,
                        "hsps": None
                    }
                    for hit in result["hits"]:
                        title = hit["description"][0]["title"]
                        for hsps in hit["hsps"]:
                            bit_score = hsps["bit_score"]
                            if bit_score >= best_match["bit_score"]:
                                best_match["bit_score"] = bit_score
                                best_match["hsps"] = hsps
                    
                    
                    differences = get_differences_protein(best_match["hsps"], name, best_match["hsps"]["gaps"])
                    return {"name": name, "differences": differences, "bit_score": best_match["bit_score"]}
            else:
                logger.error("type must be nucleotide or protein")
                sys.exit(1)

def oprD_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False):
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

    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]

    PROJECTS_PATH = config["PROJECTS_PATH"]
    BLAST_OPTIONS = config['BLAST_OPTIONS']
    BLASTN_OPTIONS = config['BLASTN_OPTIONS']
    NUCLEOTIDE_PATH = config["NUCLEOTIDE_PATH"]
    PROTEIN_PATH = config["PROTEIN_PATH"]

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    files_nucleotide = os.listdir(NUCLEOTIDE_PATH)
    

    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "oprD_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    for sample_name in samples:
        
        if direct_file:
            SPADES_FILE = sample_name
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
        else:
            sample_name = sample_name.strip()
            logger.info("Processing sample %s", sample_name)
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

        logger.info("Using SPAdes file: %s", SPADES_FILE)
        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the trimmomatic process")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            max_bitscore = {
                "name": None,
                "value": 0
            }

            for nucleotide_file in files_nucleotide:
                suffix_name = nucleotide_file[0:-len("_nucleotide.fasta")]
                name = suffix_name.replace("oprD_", "")
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
                    logger.info("***********************************************")
                    logger.info("oprD Analysis sample %s againts %s", sample_name, name)
                    logger.info("***********************************************")
                    
                    logger.info("-----------------------------------------------")
                    logger.info("--- Nucleotide analysis %s ----------------", sample_name)
                    logger.info("-----------------------------------------------")
                    results = analize_sample(output_file_nucleotide, name, "nucleotide")
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]

                    logger.info("Gaps %s", gaps)
                    logger.info("Bit score %s", bit_score)
                    
                    if bit_score > max_bitscore["value"]:
                        logger.info("------------New max bit score %s", bit_score)
                        max_bitscore["name"] = name
                        max_bitscore["value"] = bit_score
                        max_bitscore["gaps"] = gaps
                        max_bitscore["identity"] = identity
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("oprD failed assembly failed on sample %s", sample_name)
            
            logger.info("***********************************************")
            logger.info("Max bit score for %s value %s gaps %s", max_bitscore["name"], max_bitscore["value"], max_bitscore["gaps"])
            if gaps == 0 and max_bitscore["identity"] == 100:
                logger.info("This is a Wild Type (WT) sample. No gaps Rock and Roll!!!")
                logger.info("***********************************************")        
            else:
                logger.info("Analyse the differences in protein")
                logger.info("***********************************************")        
            
                protein_file = os.path.join(PROTEIN_PATH, f"oprD_{ max_bitscore['name']}_protein.fasta")
                logger.info("Use protein file %s", protein_file)
                if os.path.exists(protein_file):
                    logger.info("Using protein file: %s", protein_file)
                    output_file_protein = os.path.join(OUTPUT_PATH, f"{sample_name}_{ max_bitscore['name']}_protein.json")
                    logger.info("Output file %s", output_file_protein)
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
                        result_pro = execute_command(command_protein)
                        result = result_pro

                    if result and not normal_output:
                        # Read the json file and get the results
                        logger.info("***********************************************")
                        logger.info(" Protein Analysis sample %s againts %s", sample_name, max_bitscore['name'])
                        logger.info("***********************************************")
                        
                        results = analize_sample(output_file_protein, max_bitscore['name'], "protein")
                        logger.debug(results)
                        
                    else:
                        if not normal_output and not only_output:
                            logger.error("oprD failed assembly failed on sample %s", sample_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--normal-output', action='store_true', help='Produce only screen process normal blast outputs') 
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "oprD", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    oprD_run(project_name, config, args.parse_output, args.file, args.normal_output)

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
config = init_configs(script_directory, "PDC.json")

def get_differences(hsps, name, gaps = 0, nucleotide_protein = "nucleotide"):
        if hsps == -1:
            return []
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
                    if nucleotide_protein == "nucleotide":
                        differences.append(f"nt{i}ins{gaps}")
                    else:
                        differences.append(f"X{i+1}{h}")
            else:
                qstate = False
            if h =="-":
                if not hstate:
                    hstate = True
                    if nucleotide_protein == "nucleotide":
                        differences.append(f"nt{i}del{gaps}")
                    else:
                        differences.append(f"{q}{i+1}X")
            else:
                hstate = False
                
            # solo miramos para protein
            if nucleotide_protein == "protein":
                if m ==" ":
                    if not mstate:
                        hstate = True
                        differences.append(f"{q}{i+1}X")
                else:
                    hstate = False

        logger.debug("Differences %s %s",name,",".join(differences))
        return differences

def read_output(json_file):
    data = json.load(open(json_file))['BlastOutput2'][0]['report']
    return data

def analize_sample(json_file, name, nucleotide_protein = "nucleotide"):
    protein_data = read_output(json_file)
    
    logger.debug("Program: %s version %s",protein_data['program'], protein_data['version'])
    params_str = ', '.join(f"{k}: {v}" for k, v in protein_data["params"].items())
    logger.debug("Params used: %s", params_str)

    for key, value in protein_data["results"].items():
        for result in value:
            logger.debug("Result of %s key %s", result["query_title"], key)
            if len(result["hits"]) > 1:
                logger.debug("Number of hits: %s", len(result["hits"]))
            stats_str = ', '.join(f"{k}: {v}" for k, v in result["stat"].items())
            logger.debug("Stats: %s", stats_str)
                
            query_len = result["query_len"]
            
            if nucleotide_protein == "nucleotide":
                # sequence in N contigs check manually  
                if len(result["hits"]) > 0:
                    if len(result["hits"]) == 1:
                        hit = result["hits"][0]
                        title = hit["description"][0]["title"]
                        logger.debug("%s Hit: %s", name, title)
                        bit_score = 0
                        gaps = 100000
                        best_hsps = None
                        for hsps in hit["hsps"]:
                            gaps = hsps["gaps"]
                            identity = hsps["identity"]/hsps["align_len"]*100
                            query_from = hsps["query_from"]
                            query_to = hsps["query_to"]
                            
                            if hsps["bit_score"] > bit_score:
                                bit_score = hsps["bit_score"]
                                best_hsps = hsps
                            if query_from > 1 or query_to < query_len:
                                return {"gaps": -1, "bit_score": bit_score, "identity": -1, "hsps": [], 
                                        "differences": f"Not complete ({query_from}-{query_to})" }
                        return {"gaps": gaps, "bit_score": bit_score, "identity": identity, "hsps": best_hsps, "differences":"" }
                    else:
                        return {"gaps": -1, "bit_score": 10, "identity": -1, "hsps": [], "differences": f"sequence in {len(result['hits'])} contigs check manually" }
                else:
                    logger.debug("No hits found for %s", name)
                    return {"gaps": -1, "bit_score": -1, "identity": -1, "hsps": [], "differences": "deleted"}
            elif nucleotide_protein == "protein":
                best_match = {
                    "bit_score": 0,
                    "hsps": None
                }
                if len(result["hits"]) > 0:
                    for hit in result["hits"]:
                        title = hit["description"][0]["title"]
                        for hsps in hit["hsps"]:
                            bit_score = hsps["bit_score"]
                            if bit_score >= best_match["bit_score"]:
                                best_match["bit_score"] = bit_score
                                best_match["hsps"] = hsps
                                best_match["identity"] = hsps["identity"]/hsps["align_len"]*100
                    differences = get_differences(best_match["hsps"], name, best_match["hsps"]["gaps"], "protein")
                    return {"name": name, "differences": differences, "bit_score": best_match["bit_score"], 
                            "gaps": best_match["hsps"]["gaps"], "identity": best_match["identity"]}
                else:
                    return {"name": name, "differences": ["error"], "bit_score": 0, 
                            "gaps": 100, "identity": 0}
                
            else:
                logger.error("type must be nucleotide or protein")
                sys.exit(1)

def PDC_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False):
    ''' 
        This function is the main function to run the PDC anylisis

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
    TBLASTN_OPTIONS = config['TBLASTN_OPTIONS']
    BLASTP_OPTIONS = config['BLASTP_OPTIONS']
    PROTEIN_PATH = config["PROTEIN_PATH"]

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    files_protein = os.listdir(PROTEIN_PATH)
    

    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "PDC_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    results_data = [
        ["sample_name","PDC", "PDC_REFERENCE", "bit_score", "gaps", "identity"]
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
            logger.error("You have to run first the trimmomatic process")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            max_bitscore = {
                "name": "deleted",
                "value": -1,
                "gaps":-1,
                "identity": -1,
                "hsps": -1,
                "differences": "",
                "path": ""
            }
            
            for index, protein_file in enumerate(files_protein):
                if protein_file.find(".fasta") == -1:
                    continue
                name = protein_file[0:-len(".fasta")]
                logger.debug("Processing sample %s", name)
                protein_file = os.path.join(PROTEIN_PATH, protein_file)
                logger.debug("Using protein file: %s", protein_file)
                output_file_protein = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{name}_protein.json")
                os.makedirs(os.path.join(OUTPUT_PATH, "output"), exist_ok=True)

                logger.debug("Output file %s", output_file_protein)
                
                if normal_output:
                    command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein.replace(".json", "")] + TBLASTN_OPTIONS
                else:
                    command_protein = ["tblastn", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein, "-outfmt", "15"] + TBLASTN_OPTIONS
                
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
                    logger.debug("***********************************************")
                    logger.info("(%s/%s)) sample %s againts %s",index, len(files_protein), sample_name, name)
                    logger.debug("***********************************************")
                    
                    results = analize_sample(output_file_protein, name, "protein")
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]

                    logger.debug("Gaps %s", gaps)
                    logger.debug("Bit score %s", bit_score)
                    
                    if bit_score > max_bitscore["value"]:
                        logger.debug("------------New max bit score %s", bit_score)
                        max_bitscore["name"] = name
                        max_bitscore["path"] = protein_file
                        max_bitscore["value"] = bit_score
                        max_bitscore["gaps"] = gaps
                        max_bitscore["identity"] = identity
                        max_bitscore["differences"] = results["differences"]
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("PDC failed assembly failed on sample %s", sample_name)

            if max_bitscore["gaps"] == -1:
                logger.warning("oprD failed assembly failed on sample %s", sample_name)
                logger.warning(results)
                results_data.append([sample_name, results["differences"], max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                logger.info("This is a Wild Type (WT) sample. No gaps Rock and Roll!!!")
                logger.info("***********************************************")  
                results_data.append([sample_name, "WT", max_bitscore["name"],  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] < 100:
                logger.info("Analyse the differences in protein %s", max_bitscore['name'])
                logger.debug("***********************************************") 
                protein_file = max_bitscore["path"]
                if os.path.exists(protein_file):
                    logger.debug("Using protein file: %s", protein_file)
                    output_file_protein = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{ max_bitscore['name']}_protein.json")
                    os.makedirs(os.path.join(OUTPUT_PATH, "output"), exist_ok=True)
                
                    logger.debug("Output file %s", output_file_protein)
                    if normal_output:
                        command_protein = ["blastp", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein.replace(".json", "")] + BLASTP_OPTIONS
                    else:
                        command_protein = ["blastp", "-query", protein_file, "-subject", SPADES_FILE, "-out", output_file_protein, "-outfmt", "15"] + BLASTP_OPTIONS
                    
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
                            logger.info("This is a Wild Type (WT) sample. No gaps Rock and Roll!!!")
                            results_data.append([sample_name, "WT", results["name"],  results["bit_score"], results["gaps"], results["identity"]])
                        else:
                            result = [sample_name, ",".join(results["differences"]), max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]]
                            results_data.append(result)
                    else:
                        if not normal_output and not only_output:
                            logger.error("PDC failed assembly failed on sample %s", sample_name)

    # Crear y escribir en el archivo CSV usando punto y coma como delimitador
    filename = os.path.join(OUTPUT_PATH, f"{project_name}_PDC_results.csv")
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
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "PDC", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    PDC_run(project_name, config, args.parse_output, args.file, args.normal_output)

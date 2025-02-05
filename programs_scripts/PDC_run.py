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
import re
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "PDC.json", required_keys=["TBLASTN_OPTIONS", "PROTEIN_PATH"])

def get_differences(hsps, name, gaps = 0):
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
                    differences.append(f"X{i+1}{h}")
            else:
                qstate = False
            if h == "-":
                if not hstate:
                    hstate = True
                    differences.append(f"{q}{i+1}X")
            else:
                hstate = False
                
            if m ==" ":
                if not mstate:
                    hstate = True
                    differences.append(f"{q}{i+1}{h}")
            elif m =="+":
                if not mstate:
                    hstate = True
                    differences.append(f"{q}{i+1}{h}")
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
                    differences = get_differences(best_match["hsps"], name, best_match["hsps"]["gaps"])
                    return {"name": name, "differences": differences, "bit_score": best_match["bit_score"], 
                            "gaps": best_match["hsps"]["gaps"], "identity": best_match["identity"]}
                else:
                    return {"name": name, "differences": ["error"], "bit_score": 0, 
                            "gaps": 100, "identity": 0}
                
            else:
                logger.error("type must be nucleotide or protein")
                sys.exit(1)

def PDC_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False, extra_config={"force": False, "keep_output": False}):
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
    
    if direct_file is None and extra_config["file"] is not None:
        direct_file = extra_config["file"]

    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]
    
    PROJECTS_PATH = config["PROJECTS_PATH"]
    TBLASTN_OPTIONS = config['TBLASTN_OPTIONS']
    PROTEIN_PATH = config["PROTEIN_PATH"]

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    files_protein = os.listdir(PROTEIN_PATH)
    # Get only files with extension .fasta
    files_protein = [file for file in files_protein if file.endswith(".fasta")]
    
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
        
        if direct_file:
            SPADES_FILE = sample_name.strip()
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
        else:
            sample_name = sample_name.strip()
            logger.debug("Processing sample %s", sample_name)
            SPADES_FILE = os.path.join(SPADES_FILES_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")
            
        logger.info("********* %s ***********************", sample_name)   
        logger.debug("Using SPAdes file: %s", SPADES_FILE)
        execute = True
        if not os.path.exists(SPADES_FILE):
            execute = False
            logger.error("You have to run first the spades process")
            logger.error("This file does not exist: %s", SPADES_FILE)
        
        if execute:
            PDC1 = {}
            max_bitscore = {
                "name": "deleted",
                "value": -1,
                "gaps":-1,
                "identity": -1,
                "hsps": -1,
                "differences": "",
                "path": ""
            }
            
            pattern = r"(PDC-\d+)"
            # Usar re.search para encontrar el patrón
            files = []    
            for name in files_protein:
                match = re.search(pattern, name)                        
                pdc_name = match.group(0) if match else None
                files.append([name, pdc_name])
                
            # sort files by pdc_name
            files = sorted(files, key=lambda x: int(x[1].split("-")[1]))
            
            for index, protein_file in enumerate(files):
                protein_file, name = files[index]

                if protein_file.find(".fasta") == -1:
                    continue
                
                logger.debug("Processing sample %s", name)
                protein_file = os.path.join(PROTEIN_PATH, protein_file)
                logger.debug("Using protein file: %s", protein_file)
                output_file_protein = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{name}.json")
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
                        logger.error("You have to run first the PDC process")
                        result = False
                else:
                    result_pro = execute_command(command_protein)
                    result = result_pro
                
            
                if result and not normal_output:
                    # Read the json file and get the results
                    results = analize_sample(output_file_protein, name, "protein")
                    
                    logger.debug("***********************************************")
                    logger.info("(%s/%s)) gaps: %s identity:%.2f  pdc %s againts %s ",index, len(files_protein), results["gaps"], results["identity"], name, sample_name)
                    logger.debug("***********************************************")
                    
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]

                    logger.debug("Gaps %s", gaps)
                    logger.debug("Bit score %s", bit_score)
                    logger.debug("Identity %s", identity)
                    
                    if name == "PDC-1":
                        logger.debug("PDC-1 found")
                        PDC1["path"] = protein_file
                        PDC1["value"] = bit_score
                        PDC1["gaps"] = gaps
                        PDC1["identity"] = identity
                        PDC1["differences"] = results["differences"]
                        
                    if bit_score >= max_bitscore["value"]:
                        logger.debug("------------New max bit score %s", bit_score)
                        max_bitscore["name"] = name
                        max_bitscore["path"] = protein_file
                        max_bitscore["value"] = bit_score
                        max_bitscore["gaps"] = gaps
                        max_bitscore["identity"] = identity
                        max_bitscore["differences"] = results["differences"]
                    
                    if max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                        logger.info("***********************************************")
                        logger.info("PDC found '%s' on sample %s", name, sample_name)
                        logger.info("***********************************************")
                        break
                else:
                    if not normal_output and not only_output:
                        logger.error("PDC failed assembly failed on sample %s", sample_name)

            if max_bitscore["gaps"] == -1:
                logger.warning("PDC failed assembly failed on sample %s", sample_name)
                logger.warning(results)
                results_data.append([sample_name, ",".join(results["differences"]), max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]])
            
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                logger.info("***********************************************")  
                results_data.append([sample_name, ",".join(PDC1["differences"]), max_bitscore["name"],  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] < 100:
                results_data.append([sample_name, ",".join(PDC1["differences"]), "new type",  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            
            # Crear y escribir en el archivo CSV usando punto y coma como delimitador
            filename = os.path.join(OUTPUT_PATH, f"{project_name}_PDC_results.csv")
            with open(filename, mode='w', newline='', encoding='utf-8') as file:
                writer = csv.writer(file, delimiter=';')
                writer.writerows(results_data)
    
    logger.info("PDC analysis finished")
    if not extra_config["keep_output"]:
        os.system(f"rm -r {os.path.join(OUTPUT_PATH, 'output')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--normal-output', action='store_true', help='Produce only screen process normal blast outputs') 
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Keep the output')
                        
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json", required_keys=["TBLASTN_OPTIONS", "PROTEIN_PATH"])

    # Start the python logging variable to generate a file
    configure_logs(project_name, "PDC", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    PDC_run(project_name, config, args.parse_output, args.file, args.normal_output, extra_config={"force": args.force, "keep_output": args.keep_output})    

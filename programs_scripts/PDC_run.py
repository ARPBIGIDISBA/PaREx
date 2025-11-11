"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the PDC program to compare the nucleotide sequences with the protein sequences
"""

import os
import sys
import argparse
import logging
import json
import csv
import re
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
from typing import List, Tuple, Optional
from modules.blast_functions import analize_sample, get_differences


logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory, "PDC.json", required_keys=["TBLASTN_OPTIONS"])

def _parse_del(token: str) -> Optional[Tuple[int, int, str, str]]:
    """
    Devuelve (ini_pos, fin_pos, ini_AA, fin_AA) si es borrado, si no None.
    Acepta: "G229-", "G229-Y230del"
    """
    s = token.strip()

    m1 = re.fullmatch(r'([A-Z])(\d+)-', s)
    if m1:
        aa, p = m1.group(1), int(m1.group(2))
        return (p, p, aa, aa)

    m2 = re.fullmatch(r'([A-Z])(\d+)-([A-Z])(\d+)del', s)
    if m2:
        aa1, p1, aa2, p2 = m2.group(1), int(m2.group(2)), m2.group(3), int(m2.group(4))
        if p2 < p1:
            p1, p2, aa1, aa2 = p2, p1, aa2, aa1
        return (p1, p2, aa1, aa2)

    return None

def _fmt_del(p1: int, p2: int, aa1: str, aa2: str) -> str:
    return f"{aa1}{p1}del" if p1 == p2 else f"{aa1}{p1}-{aa2}{p2}del"

def merge_deletions_preserve(items: List[str]) -> List[str]: 
    """
    Funde borrados consecutivos/solapados en un solo rango.
    Mantiene en la lista todo lo demás sin tocar.
    Reemplaza el primer elemento del run por el rango mergeado y elimina el resto del run.
    """
    out = items[:]  # copia para modificar por índice
    run_active = False
    run_start_idx = None
    run_p1 = run_p2 = None
    run_aa1 = run_aa2 = None

    for i, tok in enumerate(items):
        parsed = _parse_del(tok)
        if parsed:
            p1, p2, aa1, aa2 = parsed

            if not run_active:
                # iniciar run
                run_active = True
                run_start_idx = i
                run_p1, run_p2 = p1, p2
                run_aa1, run_aa2 = aa1, aa2
                # el primer token se reemplazará al cerrar el run
            else:
                # ¿continúa/solapa?
                if p1 <= run_p2 + 1:
                    # extender
                    if p2 > run_p2:
                        run_p2, run_aa2 = p2, aa2
                    out[i] = None  # eliminar este token, quedará absorbido
                else:
                    # cerrar run anterior y abrir otro
                    out[run_start_idx] = _fmt_del(run_p1, run_p2, run_aa1, run_aa2)
                    run_start_idx = i
                    run_p1, run_p2 = p1, p2
                    run_aa1, run_aa2 = aa1, aa2
        else:
            # token no-borrado: cerrar run si estaba abierto
            if run_active:
                out[run_start_idx] = _fmt_del(run_p1, run_p2, run_aa1, run_aa2)
                run_active = False
                run_start_idx = None
            # dejar tok tal cual

    # fin de lista: cerrar run si quedó abierto
    if run_active:
        out[run_start_idx] = _fmt_del(run_p1, run_p2, run_aa1, run_aa2)

    # limpiar eliminados
    return [t for t in out if t not in (None, "")]

def PDC_run(project_name, config=config, direct_file = None, extra_config={"force": False, "keep_output": True}):
    ''' 
        This function is the main function to run the PDC anylisis

        parameters:
            project_name (str): Name of the project
            config dict (dict): readed from Path  is resfinder.json
            direct_file (str): path to the file instead of list in sample list
            extra_config (dict): Extra config to set the parameters
                force (bool): Force
                keep_output (bool): Keep the output
                protein (bool): Set the flag to work with protein
                nucleotide (bool): Set the flag to work with nucleotide
        results:

    '''
    
    if direct_file is None and extra_config["file"] is not None:
        direct_file = extra_config["file"]

    if "nucleotide" not in extra_config.keys() and "protein" not in extra_config.keys():
        extra_config["nucleotide"] = True
        extra_config["protein"] = False
    
    if extra_config["protein"] == False:
        extra_config["nucleotide"] = True
    elif extra_config["nucleotide"] == False:
        extra_config["protein"] = True

    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]
    
    PROJECTS_PATH = config["PROJECTS_PATH"]
    TBLASTN_OPTIONS = config['TBLASTN_OPTIONS']

    
    # Create project directory in case it is not created, read files and create output directory
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)
    SPADES_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "denovo_assemblies_SPAdes")
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "PDC_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)


    # Read PDC Database
    PDC_DATABASE_PATH = os.path.join(config['DATABASE_PATH'], 'PDC', 'PDCs_seq')

    # we iterate over the files in the nucleotide path and the search the associate protein file in the protein path
    files_protein = os.listdir(PDC_DATABASE_PATH)
    # Get only files with extension .fasta
    files_protein = [file for file in files_protein if file.endswith(".fasta")]
    pattern = r"(PDC-\d+)"
    # Usar re.search para encontrar el patrón
    files = []    
    for name in files_protein:
        match = re.search(pattern, name)                        
        pdc_name = match.group(0) if match else None
        files.append([name, pdc_name])
        
    # sort files by pdc_name
    files = sorted(files, key=lambda x: int(x[1].split("-")[1]))

    results_data = [
        ["sample_name","PDC", "PDC_REFERENCE", "bit_score", "gaps", "identity"]
    ]
    filename = None
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
            
            
            
            for index, file in enumerate(files):
                file, name = files[index]

                if index !=0 and index >0:
                    # logger.warning("Skipping file %s with index %s", file, index)
                    continue
                
                if file.find(".fasta") == -1:
                    continue
                
                logger.debug("Processing sample %s", name)
                file = os.path.join(PDC_DATABASE_PATH, file)
                logger.debug("Using protein file: %s", file)
                
                output_file = os.path.join(OUTPUT_PATH, "output", f"{sample_name}_{name}.json")
                os.makedirs(os.path.join(OUTPUT_PATH, "output"), exist_ok=True)
                logger.debug("Output file %s", output_file)
                
                # This is for nucleotide
                if extra_config["nucleotide"]:
                    command_protein = ["tblastn", "-query", file, "-subject", SPADES_FILE, "-out", output_file, "-outfmt", "15"] + TBLASTN_OPTIONS
                # This is for protein
                elif extra_config["protein"]:
                    command_protein = ["blastp", "-query", file, "-subject", SPADES_FILE, "-out", output_file, "-outfmt", "15"] + TBLASTN_OPTIONS
                
                if extra_config["force"] or not os.path.exists(output_file):
                    logger.debug("Executing command: %s", " ".join(command_protein))
                    result_pro = execute_command(command_protein)
                    result = True
                else:
                    if os.path.exists(output_file):
                        logger.debug("Output file %s already exists", output_file)
                        result = True
                    else:
                        logger.error("Output file %s does not exist", output_file)
                        result = False
            
                if result:
                    # Read the json file and get the results
                    results = analize_sample(output_file, name, "protein")
                    
                    logger.debug("***********************************************")
                    logger.debug("(%s/%s)) gaps: %s identity:%.2f  pdc %s againts %s ",index, len(files_protein), results["gaps"], results["identity"], name, sample_name)
                    logger.debug("***********************************************")
                    
                    gaps = results["gaps"]
                    bit_score = results["bit_score"]
                    identity = results["identity"]

                    logger.debug("Gaps %s", gaps)
                    logger.debug("Bit score %s", bit_score)
                    logger.debug("Identity %s", identity)
                    
                    if name == "PDC-1":
                        logger.debug("PDC-1 found")
                        PDC1["path"] = file
                        PDC1["value"] = bit_score
                        PDC1["gaps"] = gaps
                        PDC1["identity"] = identity
                        PDC1["differences"] = results["differences"]
                        
                    if bit_score >= max_bitscore["value"]:
                        # read text from protein file
                        with open(file, "r") as f:
                            # From first line >WP_063864573.1 extended-spectrum class C beta-lactamase PDC-2 [Pseudomonas aeruginosa] extract WP_063864573.1 
                            protein_text = f.readline()
                            # Example >WP_063864573.1 extended-spectrum class C beta-lactamase PDC-2 [Pseudomonas aeruginosa] extract WP_063864573.1 firts oart if split(" ")[0]
                            full_name = "{} ({})".format(name, protein_text.split(" ")[0][1:])

                        # Save the max bitscore
                        max_bitscore["name"] = full_name
                        max_bitscore["path"] = file
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
                    logger.error("PDC failed assembly failed on sample %s", sample_name)

            differences = PDC1.get("differences", [])
            logger.debug("PDC-1 differences %s", ",".join(differences))
            merged_differences = merge_deletions_preserve(differences)
            logger.debug("PDC-1 merged differences %s", ",".join(merged_differences))
            PDC1["differences"] = merged_differences
            
            if max_bitscore["gaps"] == -1:
                logger.warning("PDC failed assembly failed on sample %s", sample_name)
                results_data.append([sample_name, ",".join(results["differences"]), max_bitscore["name"], results["bit_score"], results["gaps"], results["identity"]])
            elif max_bitscore["gaps"] == 0 and max_bitscore["identity"] == 100:
                logger.info("***********************************************")  
                results_data.append([sample_name, ",".join(PDC1["differences"]), max_bitscore["name"],  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            else:
                results_data.append([sample_name, ",".join(PDC1["differences"]), "new type",  max_bitscore["value"], max_bitscore["gaps"], max_bitscore["identity"]])
            
            # Crear y escribir en el archivo CSV usando punto y coma como delimitador
            filename = os.path.join(OUTPUT_PATH, f"{project_name}_PDC_results.csv")
            with open(filename, mode='w', newline='', encoding='utf-8') as file:
                writer = csv.writer(file, delimiter=';')
                writer.writerows(results_data)
    
    if filename is not None:
        logger.info("PDC analysis finished: result in %s", filename)
        if not extra_config["keep_output"]:
            os.system(f"rm -r {os.path.join(OUTPUT_PATH, 'output')}")

    else:
        logger.error("PDC analysis finished but no results found")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Keep the output')
    # Define argumnents if we want to work in protein or nucleotide default is protein
    
    parser.add_argument('--protein', action='store_true', default=False, help='Set the flag to work with protein')
    parser.add_argument('--nucleotide', action='store_true', default=True, help='Set the flag to work with nucleotide')
                            
    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json", required_keys=["TBLASTN_OPTIONS", "PROTEIN_PATH"])

    # Start the python logging variable to generate a file
    configure_logs(project_name, "PDC", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    PDC_run(project_name, config, args.file, extra_config={"force": args.force, "keep_output": args.keep_output, "protein": args.protein, "nucleotide": args.nucleotide, "file": args.file})    




# # To delete, only to check differences function with new blast parser
# def get_differences2(hsps, name, gaps = 0):
#         if hsps == -1:
#             return []
#         qseq = hsps["qseq"]
#         midline = hsps["midline"]
#         hseq = hsps["hseq"]
#         differences = []
#         qstate = False
#         hstate = False
#         mstate = False
#         for i, (q, m, h) in enumerate(zip(qseq, midline, hseq)):
#             if q =="-":
#                 if not qstate:
#                     qstate = True
#                     # differences.append(f"X{i+1}{h}")
#             else:
#                 qstate = False
#             if h == "-":
#                 if not hstate:
#                     hstate = True
#                     # differences.append(f"{q}{i+1}X")
#             else:
#                 hstate = False
                
#             if m ==" ":
#                 if not mstate:
#                     hstate = True
#                     differences.append(f"{q}{i+1}{h}")
#             elif m =="+":
#                 if not mstate:
#                     hstate = True
#                     differences.append(f"{q}{i+1}{h}")
#             else:
#                 hstate = False

#         logger.debug("Differences %s %s",name,",".join(differences))
#         return differences

# # To delete, only to check differences function with new blast parser
# def read_output(json_file):
#     try:
#         data = json.load(open(json_file))['BlastOutput2'][0]['report']
#     except json.decoder.JSONDecodeError as e:
#         logger.error("Error decoding JSON: %s", e)
#         logger.error("File %s is empty or not a valid JSON file", json_file)
#         data = {}
#     return data

# # To delete, only to check differences function with new blast parser
# def analize_sample2(json_file, name, nucleotide_protein = "nucleotide"):
#     protein_data = read_output(json_file)
#     if not protein_data:
#         logger.error("Error reading json file %s", json_file)
#         # Remove the file if it is empty
#         os.remove(json_file)
#         return {"gaps": -1, "bit_score": -1, "identity": -1, "hsps": [], "differences": "deleted"}
#     logger.debug("Program: %s version %s",protein_data['program'], protein_data['version'])
#     params_str = ', '.join(f"{k}: {v}" for k, v in protein_data["params"].items())
#     logger.debug("Params used: %s", params_str)

#     for key, value in protein_data["results"].items():
#         for result in value:
#             logger.debug("Result of %s key %s", result["query_title"], key)
#             if len(result["hits"]) > 1:
#                 logger.debug("Number of hits: %s", len(result["hits"]))
#             stats_str = ', '.join(f"{k}: {v}" for k, v in result["stat"].items())
#             logger.debug("Stats: %s", stats_str)
                
#             query_len = result["query_len"]
            
#             if nucleotide_protein == "nucleotide":
#                 # sequence in N contigs check manually  
#                 if len(result["hits"]) > 0:
#                     if len(result["hits"]) == 1:
#                         hit = result["hits"][0]
#                         title = hit["description"][0]["title"]
#                         logger.debug("%s Hit: %s", name, title)
#                         bit_score = 0
#                         gaps = 100000
#                         best_hsps = None
#                         for hsps in hit["hsps"]:
#                             gaps = hsps["gaps"]
#                             identity = hsps["identity"]/hsps["align_len"]*100
#                             query_from = hsps["query_from"]
#                             query_to = hsps["query_to"]
                            
#                             if hsps["bit_score"] > bit_score:
#                                 bit_score = hsps["bit_score"]
#                                 best_hsps = hsps
#                             if query_from > 1 or query_to < query_len:
#                                 return {"gaps": -1, "bit_score": bit_score, "identity": -1, "hsps": [], 
#                                         "differences": f"Not complete ({query_from}-{query_to})" }
#                         return {"gaps": gaps, "bit_score": bit_score, "identity": identity, "hsps": best_hsps, "differences":"" }
#                     else:
#                         return {"gaps": -1, "bit_score": 10, "identity": -1, "hsps": [], "differences": f"sequence in {len(result['hits'])} contigs check manually" }
#                 else:
#                     logger.debug("No hits found for %s", name)
#                     return {"gaps": -1, "bit_score": -1, "identity": -1, "hsps": [], "differences": "deleted"}
#             elif nucleotide_protein == "protein":
#                 best_match = {
#                     "bit_score": 0,
#                     "hsps": None
#                 }
#                 if len(result["hits"]) > 0:
#                     for hit in result["hits"]:
#                         title = hit["description"][0]["title"]
#                         for hsps in hit["hsps"]:
#                             bit_score = hsps["bit_score"]
#                             if bit_score >= best_match["bit_score"]:
#                                 best_match["bit_score"] = bit_score
#                                 best_match["hsps"] = hsps
#                                 best_match["identity"] = hsps["identity"]/hsps["align_len"]*100
#                     differences = get_differences(best_match["hsps"], name, best_match["hsps"]["gaps"])
#                     return {"name": name, "differences": differences, "bit_score": best_match["bit_score"], 
#                             "gaps": best_match["hsps"]["gaps"], "identity": best_match["identity"]}
#                 else:

#                     logger.debug("No hits found for %s", name)
#                     return {"name": name, "differences": ["error"], "bit_score": 0, 
#                             "gaps": 100, "identity": 0}
                
#             else:
#                 logger.error("type must be nucleotide or protein")
#                 sys.exit(1)

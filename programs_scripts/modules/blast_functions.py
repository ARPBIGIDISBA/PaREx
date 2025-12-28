import json
import logging
import os
import subprocess
import sys


logger = logging.getLogger(__name__)
def get_differences(hsps, name, gaps=0, nucleotide_protein= "nucleotide"):
    """ 
    This function is used to get the differences between the query and the subject
    """
    if hsps == -1 or "qseq" not in hsps or "hseq" not in hsps or "midline" not in hsps:
        return []
    
    qseq = hsps["qseq"]
    midline = hsps["midline"]
    hseq = hsps["hseq"]
    differences = []
    qstate = False
    hstate = False
    mstate = False
    query_from = hsps.get("query_from", 1)
    for i, (q, m, h) in enumerate(zip(qseq, midline, hseq)):
        q = q.upper()
        m = m.upper()
        h = h.upper()
        if query_from>1:
            index = i + query_from -1
        else:
            index = i
        if q == "-":
            if not qstate:
                qstate = True
                if nucleotide_protein == "nucleotide":
                    differences.append(f"nt{index}ins{gaps}")
                else:
                    differences.append(f"X{index+1}{h}")
        else:
            qstate = False
        if h =="-":
            if not hstate:
                hstate = True
                if nucleotide_protein == "nucleotide":
                    differences.append(f"nt{index}del{gaps}")
                # else:
                #     print(f"{q}{index+1}X" )
                #     # differences.append(f"{q}{index+1}X")
        else:
            hstate = False
            
        # solo miramos para protein
        if nucleotide_protein == "protein" or nucleotide_protein == "nucleotide":
            if h=='*':
                if not mstate:
                    hstate = True
                    differences.append(f"{q}{index+1}X")
            elif m==" " or m=="+":
                if not mstate:
                    hstate = True
                    differences.append(f"{q}{index+1}{h}")
            else:
                hstate = False

    return differences


def read_output(json_file):
    """ This function is used to read the output json file from blast"""
    try:
        data = json.load(open(json_file))['BlastOutput2'][0]['report']
    except Exception as e:
        logger.error("Error reading json file %s", json_file)
        logger.error(e)
        data = None
    return data


def analize_sample(json_file, name, nucleotide_protein = "nucleotide", cover_limit=90):
    """
    This function is used to analyze the sample using the provided json file
    """
    protein_data = read_output(json_file)
    logger.debug("Analyzing sample %s", name)
    logger.debug("Protein data read from json file %s", json_file)
    
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
                best_match = {
                    "bit_score": 0,
                    "hsps": None,
                    "identity": 0,
                    "cover": 0,
                    "differences": []
                }

                if len(result["hits"]) > 0:
                    if len(result["hits"]) > 1:
                        logger.warning("%s has multiple contigs (%d hits)", name, len(result["hits"]))
                    bit_score = -1
                    
                    for hit in result["hits"]:
                        for hsps in hit["hsps"]:
                            bit_score = hsps["bit_score"]
                            if bit_score >= best_match["bit_score"]:
                                best_match["bit_score"] = bit_score
                                best_match["hsps"] = hsps
                                best_match["identity"] = hsps["identity"]/hsps["align_len"]*100
                                best_match["cover"] = hsps["align_len"]/result["query_len"]*100

                    if best_match:

                        if best_match["cover"]<cover_limit:
                            logger.warning("Coverage %.2f%% below limit %.2f%% for %s", best_match["cover"], cover_limit, name)
                            return {"gaps": -1, "bit_score": -1, "identity": -1, "hsps": [], "differences": "deleted"}

                        # if query_from > 1 or query_to < query_len:
                        #     return {
                        #         "gaps": -1,
                        #         "bit_score": bit_score,
                        #         "identity": -1,
                        #         "hsps": [],
                        #         "differences": f"Not complete ({query_from}-{query_to})"
                        #     }
                        return {
                            "gaps": best_match["hsps"]["gaps"],
                            "bit_score": best_match["bit_score"],
                            "identity": best_match["identity"],
                            "hsps": best_match["hsps"],
                            "differences": get_differences(best_match["hsps"], name, best_match["hsps"]["gaps"], "nucleotide")
                        }
                else:
                    logger.info("No hits found for %s", name)
                    return {"gaps": -1, "bit_score": -1, "identity": -1, "hsps": [], "differences": "deleted"}
                
            elif nucleotide_protein == "protein":
                best_match = {
                    "bit_score": 0,
                    "hsps": None,
                    "identity": 0,
                    "cover": 0,
                    "differences": []
                }
                if len(result["hits"]) > 0:
                    for hit in result["hits"]:
                        for hsps in hit["hsps"]:
                            bit_score = hsps["bit_score"]
                            if bit_score >= best_match["bit_score"]:
                                best_match["bit_score"] = bit_score
                                best_match["hsps"] = hsps
                                best_match["identity"] = hsps["identity"]/hsps["align_len"]*100
                                best_match["cover"] = hsps["align_len"]/result["query_len"]*100
                    if best_match["cover"]<cover_limit:
                        return {"gaps":  best_match["hsps"]["gaps"], "bit_score": best_match["bit_score"] , "identity": best_match["identity"], "hsps": [], "differences": ["deleted"]}
                    
                    differences = get_differences(best_match["hsps"], name, best_match["hsps"]["gaps"], "protein")
                    return {"name": name, "differences": differences, "bit_score": best_match["bit_score"], 
                            "gaps": best_match["hsps"]["gaps"], "identity": best_match["identity"]}
                else:
                    logger.debug("No hits found for %s", name)
                    logger.debug("Differences are %s", differences)
                    differences = ["deleted"]
                    return {"name": name, "differences": differences, "bit_score": best_match["bit_score"], 
                            "gaps": best_match["hsps"]["gaps"], "identity": best_match["identity"]}
                    
            else:
                logger.error("type must be nucleotide or protein")
                sys.exit(1)



def run_blast(sample_name, query_name, query_file, OUTPUT_PATH, SPADES_FILE, BLAST_OPTIONS=[], normal_output=False, only_output=False, tblastn=False):
        """
        Run analysis using tblastn
        
        Args:
            sample_name (str): Name of the sample
            query_name (str): Name of the query
            query_file (str): Path to the query file
            OUTPUT_PATH (str): Path for output files
            SPADES_FILE (str): Path to SPAdes assembly file
            BLASTN_OPTIONS (list): BLAST options as a list of strings  default empty list
            normal_output (bool): Flag for normal output format default True
            only_output (bool): Flag to only process existing output files
            tblastn (bool): Flag for tblastn vs blastn (default False)
        
        Returns:
            bool: True if analysis was successful, False otherwise
        """
        
        logger.debug("Use file %s", query_file)
        
        if os.path.exists(query_file):
            logger.debug("Using file: %s", query_file)
            
            os.makedirs(os.path.join(OUTPUT_PATH, "outputs"), exist_ok=True)
            outputfile = os.path.join(OUTPUT_PATH, "outputs", f"{sample_name}_{query_name}.json")
        
            logger.debug("Output file %s", outputfile)
            if tblastn:
                command = "tblastn"
            else:
                command = "blastn"
            if normal_output:
                command_full = [command, "-query", query_file, "-subject", SPADES_FILE, "-out", outputfile.replace(".json", "")] + BLAST_OPTIONS
            else:
                command_full = [command, "-query", query_file, "-subject", SPADES_FILE, "-out", outputfile, "-outfmt", "15"] + BLAST_OPTIONS

            logger.debug("Command line: %s", ' '.join(command_full))

            if only_output and not normal_output:
                if os.path.exists(outputfile):
                    result = True
                else:
                    logger.error("File not found: %s", outputfile)
                    logger.error("You have to run first the process")
                    result = False
            else:
                result = execute_command(command_full)
            
            return result
        else:
            logger.error("Protein file not found: %s", query_file)
            return False
        
# Execute a command and log the output
def execute_command(command):
    try:
        logger.debug(f"Executing command line: {' '.join(command)}")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        # Leer salida en tiempo real
        while True:
            output = process.stdout.readline()
            if output:
                logger.debug(output.strip())
            return_code = process.poll()
            if return_code is not None:
                break

        for output in process.stdout.readlines():
            logger.debug(output.strip())
        for output in process.stderr.readlines():
            logger.debug(output.strip())  # Log errors separately

        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return False
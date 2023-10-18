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

def print_differences_protein(hsps):
        # "qseq": "ATGAAAGTGATGAAGTGGAGCGCCATTGCACTGGCGGTTTCCGCAGGTAGCACTCAGTTCGCCGTGGCCGACGCATTCGTCAGCGATCAGGCCGAAGCGAAGGGGTTCATCGAAGACAGCAGCCTCGACCTGCTGCTCCGCAACTACTATTTCAACCGTGACGGCAAGAGCGGCAGCGGGGACCGCGTCGACTGGACCCAAGGCTTCCTCACCACCTATGAATCCGGCTTCACCCAAGGCACTGTGGGCTTCGGCGTCGATGCCTTCGGCTACCTCGGCCTGAAGCTCGACGGTACCTCTGACAAGAGCGGCACCGGCAACCTGCCAGTAATGAACGACGGAACGCCCCGTGACGACTACAGCCGCGCCGGTGGCGCCGTGAAGGTACGCATCTCCAAGACCATGCTGAAGTGGGGCGAGATGCAGCCGACCGCTCCGGTCTTCGCCGCTGGGGGCAGCCGCCTGTTCCCCCAGACCGCGACCGGCTTCCAGCTGCAGAGCAGCGAACTCGAAGGGCTCGACCTCGAGGCAGGCCACTTCACCGAGGGCAAGGAGCCGACCACCGTCAAATCGCGTGGCGAACTCTATGCCACCTACGCAGGCGAGACCGCCAAGAGCGCCGATTTCATTGGGGGCCGCTACGCAATCACCGATAACCTCAGCGCCTCCCTGTACGGCGCCGAACTCGAAGACATCTATCGCCAGTATTACCTGAACAGCAACTACACCATCCCACTGGCATCCGACCAATCGCTGGGCTTCGATTTCAACATCTACCGCACAAACGATGAAGGCAAGGCCAAGGCCGGCGACATCAGCAACACCACTTGGTCCCTGGCGGCAGCCTACACTCTGGATGCGCACACTTTCACCTTGGCCTACCAGAAGGTCCATGGCGATCAGCCGTTTGATTATATCGGCTTCGGCCGCAACGGCTCTGGCGCAGGTGGCGACTCGATTTTCCTCGCCAACTCTGTCCAGTACTCCGACTTCAACGGCCCTGGCGAGAAATCCTGGCAGGCTCGCTACGACCTGAACCTAGCCTCCTATGGCGTTCCCGGCCTGACTTTCATGGTCCGCTATATCAATGGCAAGGACATCGATGGCACCAAGATGTCTGACAACAACGTCGGCTATAAGAACTACGGCTACGGCGAGGACGGCAAGCACCACGAGACCAACCTCGAAGCCAAGTACGTGGTCCA-GTCCGGTCCGGCCAAGGACCTGTCGTTCCGCATCCGCCAGGCCTGGCACCGTGCCAACGCCGACCAGGGCGAAGGCGACCAGAACGAGTTCCGCCTGATCGTCGACTATCCGCTGTCGATCCTGTAA",
        # "hseq": "ATGAAAGTGATGAAGTGGAGCGCCATTGCACTGGCGGTTTCCGCAGGTAGCACTCAGTTCGCCGTGGCCGACGCATTCGTCAGCGATCAGGCCGAAGCGAAGGGGTTCATCGAAGACAGCAGCCTCGACCTGCTGCTCCGCAACTACTATTTCAACCGTGACGGCAAGAGCGGCAGCGGGGACCGCGTCGACTGGACCCAAGGCTTCCTCACCACCTATGAATCCGGCTTCACCCAAGGCACCGTGGGCTTCGGCGTCGATGCCTTCGGCTACCTCGGTCTGAAGCTCGACGGCACCTCGGACAAGAGCGGTACCGGCAACCTGCCGGTGATGAACGACGGCACGCCCCGTGACGACTACAGCCGCGCCGGTGGCGCCGTGAAGGTACGCATCTCCAAGACCATGTTGAAGTGGGGCGAGATGCAGCCGACCGCTCCGGTCTTCGCCGCCGGCGGCAGCCGCCTGTTCCCGCAGACCGCGACCGGCTTCCAACTGCAGAGCAGCGAACTCGAAGGGCTCGATCTCGAAGCGGGCCACTTCACCGAAGGCAAGCAGGGCACCACCACCAAGTCGCGCGGCGAACTCTACGCAACCTATGCAGGCGAGACCGCCAAGAGCGCCGATTTCATTGGGGGCCGCTACGCAATCACCGATAACCTCAGCGCCTCCCTGTACGGTGCTGAACTCGAAGACATCTATCGTCAGTATTACCTGAACAGCAACTACACCATCCCACTGGCATCCGACCAATCGCTGGGCTTCGATTTCAACATCTACCGCACAAACGATGAAGGCAAGGCCAAGGCCGGCGACATCAGCAACACCACTTGGTCCCTGGCGGCAGCCTACACTCTGGATGCGCACACTTTCACCTTGGCCTACCAGAAGGTCCATGGCGATCAGCCGTTTGATTATATCGGCTTCGGCGAGAACGGTTCCGGCGGCGGCGGTGACTCGATTTTCCTCGCCAACTCCGTGCAGTACTCCGACTTCAACGGCCCCGGCGAGAAATCCTGGCAGGCCCGCTACGACCTGAACCTCGCCTCCTATGGCGTTCCCGGCCTGACTTTCATGGTCCGCTATATCAATGGCAAGGACATCGATGGCACCAAGATGTCTGACAACAACGTCGGCTATAAGAACTACGGCTACGGCGAGGACGGCAAGCACCACGAGACCAACCTGGAAGCCAAGTACGTGGTCCACGTCCGGTCCGGCCAAGGACCTGTCGTTCCGCATCCGCCAGGCCTGGCACCGCGCCAACGCCGACCAGGCCGAAGGCGACCAGAACGAGTTCCGCCTGATCGTCGACTATCCGCTGTCGATCCTGTAA",
        # "midline": "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||| |||||||||||||| ||||| ||||||||||| |||||||||||||| || ||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||| || ||||||||||||||||| |||||||||||||||||||| ||||||||||||||||||||||||||||| ||||| || |||||||||||||| |||||| ||   ||||||  ||| ||||| ||||||||||| || ||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| || |||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   ||||| || ||||  || || ||||||||||||||||||||||| || ||||||||||||||||||||||| |||||||||||||||||||| ||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        qseq = hsps["qseq"]
        hseq = hsps["hseq"]
        midline = hsps["midline"]
        segments = []
        start = None

        for i in range(len(midline)):
            if midline[i] == ' ' or qseq[i] != hseq[i]:
                if start is None:
                    start = i
            else:
                if start is not None:
                    segments.append({
                        'start': start,
                        'end': i - 1,
                        'difference': qseq[start:i] + " -> " + hseq[start:i]
                    })
                    start = None

        if start is not None:
            segments.append({
                'start': start,
                'end': len(midline) - 1,
                'difference': qseq[start:] + " -> " + hseq[start:]
            })
        for segment in segments:

            logger.info(f"Difference: {segment['difference']} Start: {segment['start']}, End: {segment['end']}")
            start = segment['start']
            if start >1:
                start = start - 1
            end = segment['end']
            if end < len(midline) - 1:
                end = end + 2
            logger.info(f"Qseq: {qseq[start:end]}   ")
            logger.info(f"Midl: {midline[start:end]}   ")
            logger.info(f"Hseq: {hseq[start:end]}   ")

def print_metadata(json_file, name, nucleotide_protein = "nucleotide"):
    protein_data = json.load(open(json_file))['BlastOutput2'][0]['report']
    logger.info("Program: %s version %s",protein_data['program'], protein_data['version'])
    params_str = ', '.join(f"{k}: {v}" for k, v in protein_data["params"].items())
    logger.info("Params used: %s", params_str)

    for key, value in protein_data["results"].items():
        for result in value:
            logger.info("Result of %s", result["query_title"])
            logger.info("Result of %s", key)
            logger.info("Number of hits: %s", len(result["hits"]))
            stats_str = ', '.join(f"{k}: {v}" for k, v in result["stat"].items())
            logger.info("Stats: %s", stats_str)
            query_len = result["query_len"]
            
            if nucleotide_protein == "nucleotide":
                if query_len >=1296 and query_len <= 1352:
                    for hit in result["hits"]:
                        title = hit["description"][0]["title"]
                        logger.info("%s Hit: %s", name, title)
                        for hsps in hit["hsps"]:
                            gaps = hsps["gaps"]
                            
                            if gaps == 0:
                                logger.info("--- #%s", hsps['num'])
                                output_str = f"---     bit_score: {hsps['bit_score']}, evalue: {hsps['evalue']}, identity: {hsps['identity']}"
                                logger.info(output_str)
                            else:
                                logger.info("*** %s #%s %s gaps Hit: %s with ", name, hsps['num'], gaps,  title)
                                output_str = f"*** QS {hsps['query_strand']} HS {hsps['hit_strand']} Identity {hsps['identity']}"
                                logger.info(output_str)
                                output_str = f"*** From {hsps['query_from']} to {hsps['query_to']} and from {hsps['hit_from']} to {hsps['hit_to']}"
                                logger.info(output_str)
                                print_differences_protein(hsps)

            elif nucleotide_protein == "protein":
                if query_len >=441 and query_len <= 443:
                    for hit in result["hits"]:
                        title = hit["description"][0]["title"]
                        logger.info("--- %s Hit: %s", name, title)
                        for hsps in hit["hsps"]:
                            gaps = hsps["gaps"]
                            
                            if gaps == 0:

                                logger.info("--- #%s", hsps['num'])
                                output_str = f"--- bit_score: {hsps['bit_score']}, evalue: {hsps['evalue']}, identity: {hsps['identity']}"
                                logger.info(output_str)
                            else:
                                logger.info("*** %s #%s Hit: %s with %s gaps", name, hsps['num'], title, gaps)
                                output_str = f"*** bit_score: {hsps['bit_score']}, evalue: {hsps['evalue']}, identity: {hsps['identity']}"
                                logger.info(output_str)
                                output_str = f"*** Hit frame {hsps['hit_frame']} Positive {hsps['positive']} align_len {hsps['align_len']}"
                                logger.info(output_str)
                                output_str = f"*** from {hsps['query_from']} to {hsps['query_to']} and from {hsps['hit_from']} to {hsps['hit_to']}"
                                logger.info(output_str)

            else:
                logger.info("type must be nucleotide or protein")

def oprD_run(project_name, config=config, only_output = False, direct_file = None, normal_output = False):
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
    
    OUTPUT_PATH =  os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "oprD_results")
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
            SPADES_FILE = sample_name
            sample_name = os.path.basename(sample_name)
            sample_name = sample_name[0:-len(".fasta")]
            
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
                suffix_name = nucleotide_file[0:-len("_nucleotide.fasta")]
                name = suffix_name.replace("oprD_", "")
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
                protein_file = os.path.join(PROTEIN_PATH, f"{suffix_name}_protein.fasta")
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
                        logger.error("File not found: %s", output_file_nucleotide)
                        logger.error("File not found: %s", output_file_protein)
                        logger.error("You have to run first the oprD process")
                        result = False
                else:
                    result_nuc = execute_command(command_nucleotide)
                    result_pro = execute_command(command_protein)
                    result = result_nuc and result_pro

                if result and not normal_output:
                    # Read the json file and get the results
                    logger.info("***********************************************")
                    logger.info("oprD Analysis sample %s againts %s", sample_name, name)
                    logger.info("***********************************************")
                    
                    logger.info("-----------------------------------------------")
                    logger.info("--- Nucleotide analysis %s ----------------", sample_name)
                    logger.info("-----------------------------------------------")
                    print_metadata(output_file_nucleotide, name, "nucleotide")

                    # logger.info("-----------------------------------------------")
                    # logger.info("--- Protein analysis %s --------------", sample_name)
                    # logger.info("-----------------------------------------------")
                    # print_metadata(output_file_protein, name, "protein")
                    
                else:
                    if not normal_output and not only_output:
                        logger.error("oprD failed assembly failed on sample %s", sample_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process json file')   
    parser.add_argument('--file', type=str, help='Path to the file', default=None)
    parser.add_argument('--normal-output', action='store_true', help='Produce only screen process normal blast outputs')   

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    # Start the python logging variable to generate a file
    configure_logs(project_name, "oprD", config)

    logger = logging.getLogger(__name__)

    oprD_run(project_name, config, args.parse_output, args.file, args.normal_output)

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
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
import glob

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory)

def read_csv_results(csv_path, sample_id_row ="sample_name"):
    samples = {}
    if os.path.exists(csv_path):
        with open(csv_path, 'r') as file:
            reader = csv.DictReader(file, delimiter=';')
            for row in reader:
                strain_id = row[sample_id_row]
                samples[strain_id] = row
    else:
        logger.error(f"File not found {csv_path}")
        logger.warning("Execute first the analysis")
        return None
    return samples
    
def generate_excel_run(project_name, config=config):
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

    PROJECTS_PATH = config["PROJECTS_PATH"]
    OUTPUT_PATH = os.path.join(PROJECTS_PATH, project_name, f"ANALYSIS_{project_name}","")

    # Crear y escribir en el archivo CSV usando punto y coma como delimitador
    # Object with information about the samples key is STRAIN_ID and values from json files
    samples = {}


    oprD_csv = os.path.join(OUTPUT_PATH, "oprD_results",f"{project_name}_oprD_results.csv")
    oprD_samples = read_csv_results(oprD_csv)
    
    mlst_csv = os.path.join(OUTPUT_PATH, "mlst_results",f"{project_name}_mlst_results.csv")
    mlst_samples = read_csv_results(mlst_csv)
    
    snippy_csv = os.path.join(OUTPUT_PATH, "snippy_results","processed", f"*_snippy.csv")
    files = glob.glob(snippy_csv)
    
    snippy_data = []
    for file in files:
        sample_name = os.path.basename(file).replace("_snippy.csv", "")
        snippy_samples = read_csv_results(file)
        mutations = {}
        for name in snippy_samples:
            logger.info("name %s mutation %s", name, snippy_samples[name]["mutations"])
            mutations[name] = snippy_samples[name]["mutations"]
        
        
    resfinder_path = os.path.join(OUTPUT_PATH, "resfinder_results","csv_samples")
    files = glob.glob(os.path.join(resfinder_path, "*processed.csv"))
    resfinder_samples = None
    if len(files)>0:
        resfinder_samples = {}
        for file in files:  
            sample_name = os.path.basename(file).replace(".processed.csv", "")
            samples = read_csv_results(file, sample_id_row ="name")
            resfinder_samples[sample_name] = {
                "beta":[],
                "aminologlycoside":[],
                "fluoroquinolones":[],
                "other":[]
            }
            for name in samples:
                added_pheno=False
                phenotypes = samples[name].get("phenotypes","").split(",")
                for phenotype in phenotypes:
                    phenotype = phenotype.strip()
                    
                    if phenotype in ["tobramycin", "gentamycin", "amikacin"]:
                        resfinder_samples[sample_name]["aminologlycoside"].append(name)
                        added_pheno = True
                        break
                    
                    if phenotype in ["fluoroquinolones", "ciprofloxacin"]: # "cephalosporins", "penicillins" propuestos IA
                        resfinder_samples[sample_name]["fluoroquinolones"].append(name)
                        added_pheno = True
                        break
                if not added_pheno:
                    if name.startswith("bla"):
                        resfinder_samples[sample_name]["beta"].append(name)
                    else:
                        resfinder_samples[sample_name]["other"].append(name)      
            
    
    results_data = [] 
    samples = read_args(project_name, config)
    results_data_names = ["STRAIN ID"]
    header_added = {}
    for sample in samples:
        sample_id = sample.strip()
        row={}
        row["STRAIN ID"] = sample_id
        if mlst_samples:
            if "MLST" not in header_added:
                results_data_names.append("SEQUENCE TYPE") 
                header_added["MLST"] = True

            if mlst_samples[sample_id]["sequence_type"] == "-":
                row["SEQUENCE TYPE"] = mlst_samples[sample_id]["alleles"]
            else:
                row["SEQUENCE TYPE"] = mlst_samples[sample_id]["sequence_type"]
        
        if resfinder_samples:
            if "RESFINDER" not in header_added:
                results_data_names +=["ACQUIRED BETA-LACTAMASESE", "ACQUIRED AMINOGLYCOSIDE MODIFYING ENZYMES", 
                          "FLUOROQUINOLONES RESISTANCE DETERMINANTS","OTHER ACQUIRED RESISTANCE DETERMINANTS"]
                header_added["RESFINDER"] = True
        
            row["ACQUIRED BETA-LACTAMASESE"] = ",".join(resfinder_samples[sample_id]["beta"])
            row["ACQUIRED AMINOGLYCOSIDE MODIFYING ENZYMES"] =  ",".join(resfinder_samples[sample_id]["aminologlycoside"])
            row["FLUOROQUINOLONES RESISTANCE DETERMINANTS"] =  ",".join(resfinder_samples[sample_id]["fluoroquinolones"])
            row["OTHER ACQUIRED RESISTANCE DETERMINANTS"] =  ",".join(resfinder_samples[sample_id]["other"])
        
        if oprD_samples:
            if "OprD" not in header_added:
                results_data_names.append("oprD")
                results_data_names.append("oprD_REFERENCE")
                header_added["OprD"] = True
            
            oprD_sample = oprD_samples.get(sample_id, None)
            row["oprD"] = oprD_sample["oprD"]
            row["oprD_REFERENCE"] = oprD_sample["oprD_REFERENCE"]
        results_data.append(row)

    filename = os.path.join(OUTPUT_PATH, f"{project_name}_summary.csv")
    logger.info("Writing results in %s", filename)

    with open(filename, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, delimiter=';', fieldnames=results_data_names)
        writer.writeheader()
        writer.writerows(results_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "generate_excel", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    generate_excel_run(project_name)


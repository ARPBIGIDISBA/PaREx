'''
    Este script se utiliza para aplicar el programa Trimmomatic a los ficheros fastq.gz
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado dos ficheros R1_001.fastq.gz y R2_001.fastq.gz escritos en OUT
'''

import os
import shutil
import re
import argparse
import logging
from modules.general_functions import configure_logs, init_configs


logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "novasec.json")


def novasec_run(project_name, config=config):
    ''' 
        this function is used to apply the Trimmomatic program to the fastq.gz files

        parameters:
            project_name (str): Name of the project
            config_json (str): Path to the config file by default is trimmomatic_config.json

        results:
            Two files R1_001.fastq.gz and R2_001.fastq.gz written in 
            PROJECTS_PATH/project_name/ANALYSIS_{project_name}/FASTQ_Trimmomatic
            filename format will be  sample_name_trim_R1.fastq.gz and sample_name_trim_R2.fastq.gz
    '''

    # Parametros de configuración de este script
    PROJECTS_PATH = config['PROJECTS_PATH']
    
    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    FILES_PATH = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}")
    
    # read folder list of files
    files = os.listdir(FILES_PATH)
    # filter only folders
    files = [f for f in files if os.path.isdir(os.path.join(FILES_PATH, f))]

    if len(files) == 0:
        logger.warning(f"No files found in {FILES_PATH}")
        logger.warning(f"The novasec program needs the files in the following format: sample_L1_ds and sample_L2_ds")
        logger.warning(f"Please check the files in {FILES_PATH}")
        return
    logger.debug(files)
    
    pattern = re.compile(r'(?P<sample_name>\d+)_L(?P<lane>\d)_ds\..+')
    samples = {}
    for item in files:
        match = pattern.match(item)
        if match:
            sample_name = match.group('sample_name')
            lane = f"L{match.group('lane')}"
            if sample_name not in samples:
                samples[sample_name] = {}
            samples[sample_name]["sample_name"] = sample_name
            samples[sample_name][lane] = item    
    
    ## {'10062444': {'sample_name': '10062444', 'L1': '10062444_L1_ds.94ade02377de4179ae594bcb401f1f65', 'L2': '10062444_L2_ds.3e2e7bfc7d764b60a6f29a53eb43676f'}, '10004571': {'sample_name': '10004571', 'L2': '10004571_L2_ds.fe1442f809c848f5b1c9c95c7a5467f7', 'L1': '10004571_L1_ds.bbeec78802fe4aa2a77a2dd83a9577ab'}, '10027645': {'sample_name': '10027645', 'L1': '10027645_L1_ds.5049f1d2af754fd28541d3e205a9ce09', 'L2': '10027645_L2_ds.e7ff8b297bf94e29b7174eccd2a3298c'}}
    
    for sample in samples:
        logger.info(f"Processing {sample}")
        logger.debug(samples[sample]["L1"])
        logger.debug(os.path.join(FILES_PATH, samples[sample]["L1"]))
        L1_files = os.listdir(os.path.join(FILES_PATH, samples[sample]["L1"]))
        L2_files = os.listdir(os.path.join(FILES_PATH, samples[sample]["L2"]))
        L1_files = [os.path.join(FILES_PATH, samples[sample]["L1"], f) for f in L1_files]
        L2_files = [os.path.join(FILES_PATH, samples[sample]["L2"], f) for f in L2_files]
        # merge two lists with files containing R1 and R2
        merged_files = L1_files + L2_files

        # full path to the files
        content = ["_R1_", "_R2_"]
        for c in content:
            merge = [f for f in merged_files if f.find(c) > 0]
            with open(os.path.join(FILES_PATH, f"{sample}_L001{c}001.fastq.gz"), 'wb') as outfile:
                for file in merge:
                    with open(file, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile)
                        
    # for sample_name in samples:
    #     # Limpiar por si hay espacios en blanco
    #     sample_name = sample_name.strip()
    #     logger.info("Processing %s", sample_name)
        
    #     # Crear los paths de entrada y salida
    #     input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq.gz")
    #     input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq.gz")
    #     if not os.path.exists(input_r1_path) or not os.path.exists(input_r2_path):
    #         logger.error(f"The fastq.gz file for {sample_name} don't exist")
    #         logger.error(f"One of this files does not exist:\n {input_r1_path}\n {input_r2_path}")

    #     # Renombrar los ficheros de salida de 1P y 2P a R1_001 y R2_001
    #     for suffix, new_suffix in [("1P", "R1"), ("2P", "R2")]:
    #         old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #         new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{new_suffix}.fastq.gz")
    #         shutil.move(old_file_path, new_file_path)
    #         logger.info(f"Renaming file {new_file_path}")
    #         os.system(f"gunzip -f {new_file_path}")
    #         logger.info(f"Unzip file {new_file_path}")

        
    #     for suffix in ["1U", "2U"]:
    #         old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #         new_file_path = os.path.join(UNPAIRED_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #         shutil.move(old_file_path, new_file_path)
    #         logger.info(f"Storing unpaired file {new_file_path}")


    #     # Example of output_files = /home/micro/Analysis/Trimmomatic/lineage/sample/{line}.trimmed.1P.fastq.gz
    #     output_files = [os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{file}.fastq.gz") for file in ["1P", "1U", "2P", "2U"]]
    #     if all([os.path.exists(file) for file in output_files]):
    #         logger.info(f"Files for {sample_name} already exist, skipping")
    #     else:
    #         # Ejecutar Trimmomatic
    #         command = ["java", "-jar", TRIMMOMATIC_JAR_PATH, "PE",
    #                 input_r1_path, input_r2_path] + output_files + TRIMMOMATIC_OPTIONS
            
    #         result = execute_command(command)

    #         if result:
    #             # Renombrar los ficheros de salida de 1P y 2P a R1_001 y R2_001
    #             for suffix, new_suffix in [("1P", "R1"), ("2P", "R2")]:
    #                 old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #                 new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{new_suffix}.fastq.gz")
    #                 shutil.move(old_file_path, new_file_path)
    #                 logger.info(f"Renaming file {new_file_path}")
    #                 os.system(f"gunzip -f {new_file_path}")
    #                 logger.info(f"Unzip file {new_file_path}")

                
    #             for suffix in ["1U", "2U"]:
    #                 old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #                 new_file_path = os.path.join(UNPAIRED_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
    #                 shutil.move(old_file_path, new_file_path)
    #                 logger.info(f"Storing unpaired file {new_file_path}")

    #         else:
    #             logger.error("There is an error in Trimmomatic execution, check the log files")


if __name__ == "__main__":
    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config)   
    
    config["force"] = args.force
    
    configure_logs(PROJECT_NAME, "novasec", config)
    logger = logging.getLogger(__name__)
    logger.info(config)
    novasec_run(PROJECT_NAME, config)

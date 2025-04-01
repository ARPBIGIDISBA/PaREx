"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/

    This script is used to modify novaseac files to the format output to be compatible with the pipeline
"""

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


def novasec_run(project_name, config=config, extra_config={"force": False, "keep_output": False}):
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
    folders = [f for f in files if os.path.isdir(os.path.join(FILES_PATH, f))]

    if len(folders) == 0:
        logger.warning(f"No files found in {FILES_PATH}")
        logger.warning(f"The novasec program needs the files in the following format: sample_L1_ds and sample_L2_ds")
        logger.warning(f"Please check the files in {FILES_PATH}")
        return
    
    pattern = re.compile(r'(?P<sample_name>[^_]+)_L(?P<lane>\d)_ds\..+')
    samples = {}
    for item in folders:
        match = pattern.match(item)
        if match:
            sample_name = match.group('sample_name')
            lane = f"L{match.group('lane')}"
            if sample_name not in samples:
                samples[sample_name] = {}
            samples[sample_name]["sample_name"] = sample_name
            samples[sample_name][lane] = item    
    
    for sample in samples:
        if not (os.path.join(FILES_PATH, f"{sample}_L001_R1_001.fastq.gz")):    
            logger.info(f"Processing {sample}")
            logger.debug(os.path.join(FILES_PATH, samples[sample]["L1"]))
            try:
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
                
            except Exception:
                logger.warning("Error in sample")
                logger.warning(samples[sample])
        else:
            logger.warning(f"Sample already proccessed {sample}")


if __name__ == "__main__":
    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Force the execution of the program')
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config)   
    
    config["force"] = args.force
    configure_logs(PROJECT_NAME, "novasec", config)

    logger = logging.getLogger(__name__)
    novasec_run(PROJECT_NAME, config, extra_config={"force": args.force, "keep_output": args.keep_output})

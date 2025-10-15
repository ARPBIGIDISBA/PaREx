"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the Trimmomatic program to clean the reads of the samples
"""
import os
import shutil
import argparse
import logging
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs


logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "trimmomatic.json", required_keys=["TRIMMOMATIC_JAR_PATH", "TRIMMOMATIC_OPTIONS"])


def trimmomatic_run(project_name, config=config, extra_config={"force": False, "keep_output": True}):
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

    # Leer las muestras y ficheros de configuracion
    samples = read_args(project_name, config)

    # Parametros de configuración de este script
    PROJECTS_PATH = config['PROJECTS_PATH']
    TRIMMOMATIC_JAR_PATH = config['TRIMMOMATIC_JAR_PATH']

    # Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
    TRIMMOMATIC_OPTIONS = config['TRIMMOMATIC_OPTIONS']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "FASTQ_Trimmomatic")
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        sample_name = sample_name.strip()
        logger.info("Processing %s", sample_name)
        
        # Crear los paths de entrada y salida
        input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq.gz")
        input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq.gz")
        zipped = True
        if not os.path.exists(input_r1_path) or not os.path.exists(input_r2_path):
            input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq")
            input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq")
            zipped = False
            if not os.path.exists(input_r1_path) or not os.path.exists(input_r2_path):
                logger.error(f"The fastq.gz or fastq file for {sample_name} don't exist")
                logger.error(f"One of this files does not exist:\n {input_r1_path}\n {input_r2_path} \n ")
                continue
        
        if zipped:
            output_files = [os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{file}.fastq.gz") for file in ["1P", "1U", "2P", "2U"]]
        else:
            output_files = [os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{file}.fastq") for file in ["1P", "1U", "2P", "2U"]]

        
        if all([os.path.exists(file) for file in output_files]) and not extra_config["force"]:
            logger.info(f"Files for {sample_name} already exist, skipping")
        else:
            # Ejecutar Trimmomatic
            command = ["java", "-jar", TRIMMOMATIC_JAR_PATH, "PE",
                    input_r1_path, input_r2_path] + output_files + TRIMMOMATIC_OPTIONS
            
            result = execute_command(command)

            if result:
                # Renombrar los ficheros de salida de 1P y 2P a R1_001 y R2_001
                for suffix, new_suffix in [("1P", "R1"), ("2P", "R2")]:

                    if zipped:
                        old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{suffix}.fastq.gz")
                        new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{new_suffix}.fastq.gz")
                    else:
                        old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{suffix}.fastq")
                        new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{new_suffix}.fastq")
                    
                    logger.info(f"Renaming file {new_file_path}")
                    shutil.move(old_file_path, new_file_path)
                    if zipped:
                        os.system(f"gunzip -f {new_file_path}")
                        logger.info(f"Unzip file {new_file_path}")
                    

                # Remove the unpaired files or move them to a different folder
                UNPAIRED_PATH = os.path.join(OUTPUT_PATH, "UNPAIRED")
                for sufix in ["1U", "2U"]:
                    # Borrar los ficheros
                    if zipped:
                        file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{sufix}.fastq.gz")
                    else:
                        file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{sufix}.fastq")
                    if os.path.exists(file_path):
                        os.remove(file_path)      
                        logger.info(f"Deleting unpaired file {file_path}")

                # if extra_config["keep_output"]:
                #     UNPAIRED_PATH = os.path.join(OUTPUT_PATH, "UNPAIRED")
                #     os.makedirs(UNPAIRED_PATH, exist_ok=True)
                #     logger.info(f"Unpaired files will be stored in {UNPAIRED_PATH}")
                #     for suffix in ["1U", "2U"]:
                #         if zipped:
                #             old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{suffix}.fastq.gz")
                #             new_file_path = os.path.join(UNPAIRED_PATH, f"{sample_name}_trim_{suffix}.fastq.gz")
                #         else:
                #             old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{suffix}.fastq")
                #             new_file_path = os.path.join(UNPAIRED_PATH, f"{sample_name}_trim_{suffix}.fastq")

                #         shutil.move(old_file_path, new_file_path)
                #         logger.info(f"Storing unpaired file {new_file_path}")
                #         if os.path.exists(old_file_path) and not extra_config["keep_output"]:
                #             os.remove(old_file_path)      

            else:
                logger.error("There is an error in Trimmomatic execution, check the log files")


if __name__ == "__main__":
    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default="INFO")
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', default=True, help='Keep the output files')

    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config, required_keys=["TRIMMOMATIC_JAR_PATH", "TRIMMOMATIC_OPTIONS"])   
    
    
    configure_logs(PROJECT_NAME, "trimmomatic", config, extra_config={"force": args.force, "keep_output": args.keep_output})
    logger = logging.getLogger(__name__)
    logger.info(config)
    trimmomatic_run(PROJECT_NAME, config, extra_config={"force": args.force, "keep_output": args.keep_output})

'''
    Generic functions for the scripts
'''
import subprocess
import logging
import json
import os
import sys
import shutil

logger = logging.getLogger(__name__)


def validate_config(config, required_keys, config_json):
    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        logger.error(f"Missing required config keys: {', '.join(missing_keys)}")
        logger.error(f"Check the config file: {config_json}")
        sys.exit(1)
    return True

def read_config(config_json, required_keys = []):
    # Check if the .json config file exists
    if not os.path.exists(config_json):
        # If not, look for the .json.sample file
        sample_path = f"{config_json}.sample"
        
        if os.path.exists(sample_path):
            logger.warning(f"Config file '{config_json}' not found. Creating it from '{sample_path}'")
            try:
                shutil.copyfile(sample_path, config_json)
            except IOError as e:
                logger.error(f"Error creating config file from sample: {e}")
                sys.exit(1)
        else:
            logger.error(f"Neither config file '{config_json}' nor sample '{sample_path}' found.")
            sys.exit(1)
            
    with open(config_json, 'r') as file:
        config = json.load(file)
        validate_config(config, required_keys, config_json)
        return config

def check_project(project_path):
    if not os.path.exists(project_path):
        print("Project not found")
        print(f"Creating project root {project_path}")
        # Command line question if you are sure to continue
        print("Are you sure to continue? (y/n)")
        answer = input()
        if answer =="y":
            os.makedirs(project_path, exist_ok=True)
        else:
            print("Operation cancelled")
            sys.exit(1)
    if not os.path.exists(os.path.join(project_path, f"Logs")):
        os.makedirs(os.path.join(project_path, f"Logs"), exist_ok=True)

def read_args(project_name, config):
    ''' Función para leer los argumentos de la línea de comandos, el fichero de samples y la configuración del json
        arguments:
            project_name (str): Nombre del proyecto, carpeta con el mismo nombre en el PROJECTS_PATH definido en geneal.json
            fichero con el listado de samples y carpeta ANALYSIS_{project_name} en el PROJECTS_PATH con los fasta.gz de las muestras
            config (dict): Diccionario con los valores de configuración del json mezcla de general.json y bowtie.json

    '''

    # Create project directory in case it is not created
    project_path = os.path.join(config["PROJECTS_PATH"], project_name)
    os.makedirs(project_path, exist_ok=True)

    sample_file = os.path.join(project_path, f"SAMPLES_LIST_{project_name}")
    if os.path.exists(sample_file):
        with open(sample_file, 'r') as file:
            samples = file.readlines()

        return samples
    else:
        logger.error("Problem reading the sample list file ")
        logger.error(f"Not found {sample_file}")
        logger.error("Check the correct spelling of the project name")
        exit(1)


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
            logger.error(output.strip())  # Log errors separately

        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return False


def init_configs(script_directory, config_json=None):
    '''
        This function is used to initialize the config files
        parameters:
            script_directory (str): path to the directory where the script is
            config_json (str): path to the config file
        results:
            config (dict): dictionary with the config values of the specific json plus 
            the general values of the general.json file
    '''
    general_config = os.path.join(script_directory, os.path.join("configs","general.json"))
    config_general = read_config(general_config)
    if config_json:
        default_config_json = os.path.join(script_directory, os.path.join("configs", config_json))
        config = read_config(default_config_json)
    else:
        config = {}
    config["PROJECTS_PATH"] = config_general["PROJECTS_PATH"]
    config["REFERENCE_PATH"] = config_general["REFERENCE_PATH"]
    return config


def configure_logs(project_name, script_name, config, log_mode="w", log_level="INFO"):
    log_levels = {
        "INFO": logging.INFO,
        "DEBUG": logging.DEBUG,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL
    }
    level = log_levels.get(log_level.upper(), logging.INFO)
    log_path = os.path.join(config["PROJECTS_PATH"], project_name, "Logs")
    os.makedirs(log_path, exist_ok=True)
    log_file = os.path.join(log_path, f'{project_name}_{script_name}.log')

    logging.basicConfig(
        level=level,
        format='%(asctime)s-%(levelname)s- %(message)s - %(filename)s:%(lineno)d',
        datefmt='%H:%M:%S',
        handlers=[
            logging.FileHandler(log_file, mode=log_mode),
            logging.StreamHandler()
        ]
    )
    logger.info(f"Log file created: {log_file}")
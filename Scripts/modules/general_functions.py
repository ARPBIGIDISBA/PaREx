'''
    Generic functions for the scripts
'''
import subprocess
import logging
import json
import os


logger = logging.getLogger(__name__)


def read_config(config_json):
    # Leer el archivo catde configuración
    with open(config_json, 'r') as file:
        config = json.load(file)
        return config


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
        logger.error("Problem readingh the sample list file ")
        logger.error(f"Not found {sample_file}")
        logger.error("Check the correct spelling of the project name")
        exit(1)


# Execute a command and log the output
def execute_command(command):
    # Enseñar el comando que se va a ejecutar
    logger.info(f"Executing command line: {' '.join(command)}")

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    # Capturar las salidas del comando y añadirlas al fichero de salida
    # stdout, stderr = process.communicate()

    while True:
        output = process.stdout.readline()
        if output and output != "":
            logger.info(output.strip())
        # Check for termination
        return_code = process.poll()
        if return_code is not None:
            for output in process.stdout.readlines():
                logger.info(output.strip())
            for output in process.stderr.readlines():
                logger.info(output.strip())
            return (return_code == 0)


def init_configs(script_directory, config_json):
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
    default_config_json = os.path.join(script_directory, os.path.join("configs", config_json))
    config = read_config(default_config_json)
    config["PROJECTS_PATH"] = config_general["PROJECTS_PATH"]
    config["LOGS_PATH"] = config_general["LOGS_PATH"]
    config["REFERENCE_PATH"] = config_general["REFERENCE_PATH"]
    return config


def configure_logs(project_name, script_name, config, log_mode="w"):

    LOG_MODE = log_mode  # "a" to append or "w" to overwrite
    LOG_PATH = os.path.join(config["LOGS_PATH"], project_name)
    os.makedirs(LOG_PATH, exist_ok=True)
    LOG_NAME = os.path.join(LOG_PATH, f'{project_name}_{script_name}.log')                       
    LOG_FORMAT = '%(asctime)s-%(levelname)s- %(message)s - %(filename)s:%(lineno)d'
    logging.basicConfig(level=logging.DEBUG,
                        format=LOG_FORMAT,
                        datefmt='%H:%M:%S',
                        handlers=[
                            logging.FileHandler(LOG_NAME, mode=LOG_MODE),
                            logging.StreamHandler()
                        ])
    logging.info("Log file: %s", LOG_NAME)
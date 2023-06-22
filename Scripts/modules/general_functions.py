'''
    Generic functions for the scripts
'''
import subprocess
import logging
import json
import os
import traceback

logger = logging.getLogger(__name__)


def read_config(config_json):
    # Leer el archivo catde configuración
    with open(config_json, 'r') as file:
        config = json.load(file)
        return config


def read_args(project_name, config):
    ''' Función para leer los argumentos de la línea de comandos, el fichero de samples y la configuración del json'''

    # Create project directory in case it is not created
    project_path = os.path.join(config["PROJECTS_PATH"], project_name)
    os.makedirs(project_path, exist_ok=True)

    sample_file = os.path.join(project_path, f"SAMPLES_LIST_{project_name}")
    try:
        with open(sample_file, 'r') as file:
            samples = file.readlines()

        return samples
    except Exception:
        traceback.print_exc()
        logger.error("Problem readingh the sample list file ")
        logger.error(f"Not found {sample_file}")
        exit(1)


# Execute a command and log the output
def execute_command(command):
    # Enseñar el comando que se va a ejecutar
    logger.info(f"Ejecutando el comando: {' '.join(command)}")

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Capturar las salidas del comando y añadirlas al fichero de salida
    stdout, stderr = process.communicate()

    # Log the output
    if stdout:
        logger.info(stdout.decode())
    if stderr:
        logger.info(stderr.decode())

    # Return True if the process was successful
    return (process.returncode == 0)


def configure_logs(project_name, script_name, config, log_mode="w"):
    LOG_MODE = log_mode  # "a" to append or "w" to overwrite
    print(config["LOGS_PATH"], project_name)
    LOG_PATH = os.path.join(config["LOGS_PATH"], project_name)
    os.makedirs(LOG_PATH, exist_ok=True)
    LOG_NAME = os.path.join(LOG_PATH, f'{project_name}_{script_name}.log')                       
    LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s - %(name)s - %(pathname)s:%(lineno)d'
    logging.basicConfig(level=logging.DEBUG,
                        format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[
                            logging.FileHandler(LOG_NAME, mode=LOG_MODE),
                            logging.StreamHandler()
                        ])
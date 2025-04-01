"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/

This is the api.py file that will be used to create the API to execute the pipeline from a web interface.

"""
from flask import Flask, request, Response
from flask_restx import Api, Resource, Namespace, reqparse
from flask_cors import CORS
import select
import json

import os
import sys
import subprocess
import time
import logging
from werkzeug.datastructures import FileStorage 
import uuid
import shutil  # ja que ho necessitarem per eliminar directoris

# Add path Scripts to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'programs_scripts')))

logger = logging.getLogger("api_pipeline")
logger.setLevel(logging.INFO)

logger = logging.getLogger("api_pipeline")
logger.setLevel(logging.DEBUG)  # Capture all logs


# Read projects path from environment variable
PROJECTS_PATH = os.getenv("PROJECTS_PATH") or "/home/mbonet/microbiologia/Projects/"  
TEMP_FOLDER = os.getenv("TEMP_FOLDER") or "/tmp"    

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s - %(filename)s:%(lineno)d',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(
                            os.path.join(TEMP_FOLDER, f'execute_pipeline.log'),
                            mode="w"),
                        logging.StreamHandler()
                    ])


# Configuración de Flask y API
app = Flask(__name__)
cors_config = {
        'ORIGINS': [
            os.getenv("CORS_ORIGIN") or "http://localhost:5173",
        ]}

CORS(app, resources={r'/api*': {'origins': cors_config['ORIGINS']}}, supports_credentials=True)

api = Api(app, version="1.0", title="Pipeline API", prefix="/api",
          description="API para ejecutar el pipeline de análisis de datos")

ns_pipeline = Namespace('pipeline', description='Operaciones del pipeline')
api.add_namespace(ns_pipeline)


UPLOAD_FOLDER = "/tmp"  # Directorio para almacenar archivos temporales
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Definir los parámetros del formulario (incluye archivo)
pdc_parser = reqparse.RequestParser()
pdc_parser.add_argument('direct_sequence', type=str, required=False, help='Secuencia fasta en texto')
pdc_parser.add_argument('file', type=FileStorage, location='files', required=False, help='Archivo FASTA')

# ✅ Store the latest process
processes = {}

@ns_pipeline.route('/pdc')
class PDCRun(Resource):
    @ns_pipeline.expect(pdc_parser)
    @ns_pipeline.doc(
        description="Upload a FASTA file or enter a direct sequence for analysis.",
        consumes="multipart/form-data",
    )
    def post(self):
        global processes
        args = pdc_parser.parse_args()
        direct_sequence = args['direct_sequence']
        file_path = None
        process_id = str(uuid.uuid4())

        if 'file' in request.files:
            uploaded_file = request.files['file']
            if uploaded_file.filename == '':
                return {"error": "No se ha seleccionado un archivo"}, 400
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], process_id + ".fasta")
            uploaded_file.save(file_path)
            logger.info(f"📂 Archivo recibido: {uploaded_file.filename}")

        elif direct_sequence:
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], "temp_sequence.fasta")
            with open(file_path, "w") as f:
                f.write(direct_sequence)
            logger.info(f"📝 Secuencia recibida: {file_path}")

        if not file_path:
            return {"error": "Se requiere un archivo FASTA o una secuencia en texto"}, 400

        project_name = "temp"
        command = ["python3", "execute_pipeline.py", project_name, "PDC", "--file", file_path, "--force"]
        logger.info(f"🚀 Ejecutando: {' '.join(command)}")

        # ✅ Start subprocess and store process
        current_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        processes[process_id] = {
                                    "file_path": file_path,
                                    "process":  current_process
                                }
        
        logger.info(f"🔧 Procés iniciat amb ID: {process_id}")
        return {"message": "PDC execution started", "process_id": process_id}

@app.route('/api/pipeline/pdc/logs')
def stream_logs():
    global processes
    process_id = request.args.get("process_id")
    process = processes.get(process_id)
    current_process = process["process"] if process else None
    file_path = process["file_path"] if process else None
    
    def generate_logs():
        
        # try 5 seconds to get the process
        for _ in range(5):
            yield "data: Waiting for process to start\n\n"
            if current_process:
                break
            time.sleep(1)
        if not current_process:
            logger.info("No process running")
            yield "error: No process running\n\n"
        else:
            while True:
                reads = [current_process.stdout.fileno(), current_process.stderr.fileno()]
                ready, _, _ = select.select(reads, [], [])
                for fd in ready:
                    if fd == current_process.stdout.fileno():
                        output = current_process.stdout.readline().strip()
                        if output:
                            yield f"data: {output}\n\n"

                    if fd == current_process.stderr.fileno():
                        error = current_process.stderr.readline().strip()
                        if error:
                            yield f"data: {error}\n\n"

                if current_process.poll() is not None:
                    break


            output_file = os.path.join(PROJECTS_PATH, "temp/ANALYSIS_temp/PDC_results", "temp_PDC_results.csv")
            if not os.path.exists(output_file):
                yield "error: No results found\n\n"
                return
            with open(output_file, "r") as f:
                lines = f.readlines()
                headers = lines[0].strip().split(";")
                for line in lines[1:]:
                    data = line.strip().split(";")
                    data_dict = dict(zip(headers, data))
                    yield f"data:Result:{json.dumps(data_dict)}\n\n"

            # cleanup: eliminar els resultats temporals
            project_analysis_path = os.path.join(PROJECTS_PATH, "temp/ANALYSIS_temp")
            try:
                shutil.rmtree(project_analysis_path)
                logger.info(f"🧹 S'ha eliminat {project_analysis_path}")
            except Exception as e:
                logger.warning(f"⚠️ No s'ha pogut netejar {project_analysis_path}: {str(e)}")            
            
            try:
                os.remove(file_path)
                logger.info(f"🗑️ S'ha eliminat {file_path}")
            except Exception as e:
                logger.warning(f"⚠️ No s'ha pogut eliminar {file_path}: {str(e)}")

            del processes[process_id]
            logger.info(f"🔧 Procés finalitzat amb ID: {process_id}"
                        )
    return Response(generate_logs(), mimetype='text/event-stream')


@ns_pipeline.route('/pdc/stop/<string:process_id>')
class StopProcess(Resource):
    def post(self, process_id):
        global processes
        process_info = processes.get(process_id)

        if not process_info:
            return {"error": "Process not found"}, 404

        process = process_info if isinstance(process_info, subprocess.Popen) else process_info["process"]
        del processes[process_id]
        if process.poll() is None:  # Encara actiu
            process.terminate()
            logger.info(f"🛑 Procés {process_id} aturat per l'usuari.")
            return {"message": f"Process {process_id} terminated"}, 200
        else:
            return {"message": "Process already completed"}, 200
        
if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0", port=5000)

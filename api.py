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
CORS(app)
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
current_process = None

@ns_pipeline.route('/pdc')
class PDCRun(Resource):
    @ns_pipeline.expect(pdc_parser)
    @ns_pipeline.doc(
        description="Upload a FASTA file or enter a direct sequence for analysis.",
        consumes="multipart/form-data",
    )
    def post(self):
        global current_process
        args = pdc_parser.parse_args()
        direct_sequence = args['direct_sequence']
        file_path = None

        if 'file' in request.files:
            uploaded_file = request.files['file']
            if uploaded_file.filename == '':
                return {"error": "No se ha seleccionado un archivo"}, 400
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], uploaded_file.filename)
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

        return {"message": "PDC execution started"}

@app.route('/api/pipeline/pdc/logs')
def stream_logs():
    global current_process
    if not current_process:
        return {"error": "No process running"}, 400

    def generate_logs():
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
        # read final output /home/mbonet/microbiologia/Projects/temp/ANALYSIS_temp/PDC_results/temp_PDC_results.csv
        # sample_name;PDC;PDC_REFERENCE;bit_score;gaps;identity
        # PA01-DK_S26_L001.SPAdes.denovoassembly;;PDC-1;807.364;0;100.0

        output_file = os.path.join(PROJECTS_PATH, "temp/ANALYSIS_temp/PDC_results", "temp_PDC_results.csv")
        with open(output_file, "r") as f:
            lines = f.readlines()
            headers = lines[0].strip().split(";")
            for line in lines[1:]:
                data = line.strip().split(";")
                data_dict = dict(zip(headers, data))
                yield f"data:Result:{json.dumps(data_dict)}\n\n"
        

    return Response(generate_logs(), mimetype='text/event-stream')
if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0", port=5000)

from flask import Flask, request, Response
from flask_restx import Api, Resource, Namespace, reqparse
import os
import sys
import subprocess
import logging
from werkzeug.datastructures import FileStorage 

# Add path Scripts to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'programs_scripts')))
# Importar funciones del pipeline
from programs_scripts.PDC_run import PDC_run



logger = logging.getLogger("api_pipeline")
logger.setLevel(logging.INFO)

logger = logging.getLogger("api_pipeline")
logger.setLevel(logging.DEBUG)  # Capture all logs

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s - %(filename)s:%(lineno)d',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(
                            os.path.join("/tmp", f'execute_pipeline.log'),
                            mode="w"),
                        logging.StreamHandler()
                    ])


# Configuración de Flask y API
app = Flask(__name__)
api = Api(app, version="1.0", title="Pipeline API",
          description="API para ejecutar el pipeline de análisis de datos")

ns_pipeline = Namespace('pipeline', description='Operaciones del pipeline')
api.add_namespace(ns_pipeline)


UPLOAD_FOLDER = "/tmp"  # Directorio para almacenar archivos temporales
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER



# Definir los parámetros del formulario (incluye archivo)
pdc_parser = reqparse.RequestParser()
pdc_parser.add_argument('direct_sequence', type=str, required=False, help='Secuencia fasta en texto')
pdc_parser.add_argument('file', type=FileStorage, location='files', required=False, help='Archivo FASTA')

@ns_pipeline.route('/pdc')
class PDCRun(Resource):
    @ns_pipeline.expect(pdc_parser)
    @ns_pipeline.doc(
        description="Upload a FASTA file with a project name. This will show a file upload button in Swagger UI.",
        consumes="multipart/form-data",
    )
    def post(self):
        """
        Ejecuta el análisis PDC con un archivo subido o una secuencia en texto.
        """
        args = pdc_parser.parse_args()
        direct_sequence = args['direct_sequence']
        
        file_path = None

        # Si se sube un archivo, guardarlo en /tmp
        if 'file' in request.files:
            uploaded_file = request.files['file']
            
            if uploaded_file.filename == '':
                return {"error": "No se ha seleccionado un archivo"}, 400

            file_path = os.path.join(app.config['UPLOAD_FOLDER'], uploaded_file.filename)
            uploaded_file.save(file_path)
            logger.info(f"Archivo recibido: {file_path}")

        # Si se recibe una secuencia en texto, guardarla en un archivo temporal
        elif direct_sequence:
            # Write direct sequence to a file in /tmp
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], "temp_sequence.fasta")
            with open(file_path, "w") as f:
                f.write(direct_sequence)
            logger.info(f"Secuencia recibida: {file_path}")

        # Si no hay archivo ni secuencia, devolver error
        if not file_path:
            return {"error": "Se requiere un archivo FASTA o una secuencia en texto"}, 400

        # Ejecutar PDC_run con el archivo
        extra_config = {"force": True, "keep_output": False}
        project_name = "temp"
        # Execute PDC_run in a subprocess
        def stream_logs():
            print("hola")
            logger.info("Executing PDC analysis")
            command = ["python3", "execute_pipeline.py", "temp", "PDC", "--file", file_path, "--force"]
            logger.info(" ".join(command))
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Stream logs in real-time
            for line in process.stdout:
                print(line.strip())
                yield f"data: {line.strip()}\n\n"

            process.stdout.close()
            process.wait()

            # Final result message
            yield f"data: PDC analysis completed successfully.\n\n"

        return Response(stream_logs(), mimetype='text/event-stream')


if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0", port=5000)

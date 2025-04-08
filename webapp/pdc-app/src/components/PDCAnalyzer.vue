<script setup>
import { ref, computed, onMounted, onUnmounted, nextTick } from "vue";
import axios from "axios";

const API_URL = import.meta.env.VITE_API_URL;

const fileInput = ref(null);
const directSequence = ref("");
const logMessages = ref([]);
const eventSource = ref(null);
const isStreaming = ref(false);
const result = ref(null); 
const logContainer = ref(null);
const processId = ref(null);
const runing_process = ref(null);
const totalSamples = ref(0)
const isUploading = ref(false)

const logText = computed(() => logMessages.value.join("\n"));

const parseLog = (logLine) => {
  const pattern = /\((\d+)\/(\d+)\)\) gaps: (\d+) identity:([\d.]+)\s+pdc\s+(PDC-\d+)/;
  const match = logLine.match(pattern);
  if (match) {
    const [, x, y, gaps, identity, pdcName] = match;
    return {
      pdc: pdcName,
      gaps: parseInt(gaps, 10),
      identity: parseFloat(identity),
      x: parseInt(x, 10),
      y: parseInt(y, 10),
      progress: `${x}/${y}`
    };
  } else {
    return { error: "No match found" };
  }
};

const lastLogMessage = computed(() => {
  if (logMessages.value.length === 0) return "";
  return parseLog(logMessages.value[logMessages.value.length - 1])
});

const scrollToBottom = () => {
  nextTick(() => {
    if (logContainer.value) {
      logContainer.value.scrollTop = logContainer.value.scrollHeight;
    }
  });
};

const startLogStream = (pid) => {
  if (isStreaming.value) return;
  logMessages.value = [];
  isStreaming.value = true;
  processId.value = pid;

  eventSource.value = new EventSource(`${API_URL}/pipeline/pdc/logs?process_id=${pid}`);
  eventSource.value.onmessage = (event) => {
    if (event.data.includes("Result:")) {
      result.value = event.data.split("Result:")[1];
      stopLogStream();
    } else {
      logMessages.value.push(event.data);
      scrollToBottom();
    }
  };

  eventSource.value.onerror = (error) => {
    logMessages.value.push("⚠️ Error streaming logs.");
    stopLogStream();
  };
};

const stopLogStream = () => {
  if (eventSource.value) {
    fetchStats()
    eventSource.value.close();
    isStreaming.value = false;
    runing_process.value = null;
  }
};

const stopProcess = async () => {
  try {
    await axios.post(`${API_URL}/pipeline/pdc/stop/${processId.value}`);
    stopLogStream();
    result.value = null;
    logMessages.value = [];
  } catch (error) {
    console.error("❌ Error stopping the process:", error);
  }
};
const uploadFile = async () => {
  const file = fileInput.value;
  if (!file) {
    alert("📂 Select a FASTA file first");
    return;
  }

  const formData = new FormData();
  formData.append("file", file);

  try {
    isUploading.value = true
    result.value = null;
    const res = await axios.post(`${API_URL}/pipeline/pdc`, formData, {
      headers: { "Content-Type": "multipart/form-data" }
    });
    runing_process.value = res.data.process_id;
    startLogStream(res.data.process_id);
  } catch (error) {
    console.error("❌ Error uploading the file:", error);
  } finally {
    isUploading.value = false
  }
};

const sendDirectSequence = async () => {
  if (!directSequence.value.trim()) {
    alert("Enter a valid FASTA sequence");
    return;
  }

  try {
    result.value = null;
    const res = await axios.post(`${API_URL}/pipeline/pdc`, {
      direct_sequence: directSequence.value
    });
    startLogStream(res.data.process_id);
  } catch (error) {
    console.error("❌ Error sending the sequence:", error);
  }
};
const fetchStats = async () => {
  try {
    const res = await axios.get(`${API_URL}/pipeline/pdc/stats`)
    totalSamples.value = Object.keys(res.data).length
  } catch (err) {
    console.error("❌ Error loading stats:", err)
  }
}

onMounted(() => {
  fetchStats()
  const interval = setInterval(fetchStats, 5000) // cada 10 segons

  onUnmounted(() => clearInterval(interval)) // neteja l'interval si el component es desmunta
})

onUnmounted(() => stopLogStream());
</script>

<template>
  <v-container>
    <div  style="position: relative;">
      <v-toolbar color="light-grey"><v-spacer></v-spacer><h3 style="margin-right:10px"> 🔬 Total different sample analyzed: <strong>{{ totalSamples }} </strong>  </h3></v-toolbar>
      <v-overlay persistent :model-value="runing_process!=null" contained absolute opacity="0.3">
      </v-overlay>
      <v-row>
        <v-col>
          <v-card>
            <v-card-text v-if="isUploading">
              <p>Uploading...</p>
              <v-progress-linear indeterminate color="primary" />
            </v-card-text>
            <v-card-text>
              <p>Upload a FASTA file to analyze</p>
              <v-file-input
                label="Upload FASTA File"
                @change="fileInput = $event.target.files[0]"
                accept=".fasta"
                outlined
              />
              <v-btn v-if="fileInput" class="mt-2" @click="uploadFile" color="primary" :disabled="runing_process!=null">Submit File</v-btn>
            </v-card-text>
            
          </v-card>
        </v-col>
        <v-col>
          <v-card>
            <v-card-text>
              <p>Enter FASTA sequence</p>
              <v-textarea
                v-model="directSequence"
                label="Enter FASTA sequence"
                placeholder=""
                rows="4"
                outlined
              />
              <v-btn class="mt-2" @click="sendDirectSequence" color="primary" v-if="directSequence" :disabled="runing_process!=null">Analyze Sequence</v-btn>
            </v-card-text>
          </v-card>
        </v-col>
      </v-row>
    </div>
    <v-row  v-if="logText.length > 0 || result">
      <v-col>
        <v-card color="surface" elevation="10" class="mt-4" style="height: 300px; overflow-y: auto;" ref="logContainer">
          <v-card-text v-if="result" >
            <h2>Detected: {{ JSON.parse(result)["PDC_REFERENCE"] }}</h2>
            <v-table>  
              <tbody>
                <tr>
                  <th><b>Difference with PDC-1</b></th>
                  <td><b>Observations</b></td>
                </tr>
                <tr v-for="(value, key) in JSON.parse(result)['PDC'].split(',')" :key="key">
                  <td>{{ value }}</td>
                  <td></td>
                </tr>
              </tbody>
            </v-table>
            <h2> Sample name: {{ JSON.parse(result)["sample_name"] }}</h2>
          </v-card-text>
          <v-card-text v-if="!result && logText.length > 0" width="300px">
            <v-btn @click="stopProcess" color="error" class="float-right">
              <v-icon>close</v-icon>
            </v-btn>
            <v-table width="300px" style="background-color: white !important;">
              <tbody>
                <tr>
                  <th>Comparing to</th>
                  <td>{{ lastLogMessage["pdc"] }}</td>
                </tr>
                <tr>
                  <th>Progress</th>
                  <td>{{ lastLogMessage["progress"] }}</td>
                </tr>
                <tr>
                  <th>Identity</th>
                  <td>{{ lastLogMessage["identity"] }}%</td>
                </tr>
              </tbody>
            </v-table>
          </v-card-text>
        </v-card>
      </v-col>
    </v-row>
    
  </v-container>
</template>

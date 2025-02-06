<script setup>
import { ref, computed, onUnmounted, nextTick } from "vue";
import axios from "axios";

const API_URL = import.meta.env.VITE_API_URL;

const fileInput = ref(null);
const directSequence = ref("");
const logMessages = ref([]);
const eventSource = ref(null);
const isStreaming = ref(false);
const result = ref(null); 
const logContainer = ref(null);

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

const startLogStream = () => {
  if (isStreaming.value) return;

  logMessages.value = [];
  isStreaming.value = true;

  console.debug("📡 Connecting to real-time logs...");

  eventSource.value = new EventSource(`${API_URL}/pipeline/pdc/logs`);
  eventSource.value.onmessage = (event) => {
    console.debug("📩 Log received:", event.data);
    // check if result is ready then finish event.restult
    if (event.data.includes("Result:")) {
      console.debug("🎉 Result ready!");
      result.value = event.data.split("Result:")[1];
      stopLogStream();
    }else {
      logMessages.value.push(event.data);
      nextTick( () => {
        scrollToBottom();
      });
    }
  };

  eventSource.value.onerror = (error) => {
    console.error("❌ Error streaming logs:", error);
    logMessages.value.push("⚠️ Error streaming logs.");
    stopLogStream();
  };
};

const stopLogStream = () => {
  if (eventSource.value) {
    eventSource.value.close();
    isStreaming.value = false;
    console.debug("🔌 Stream closed.");
  }
};

const uploadFile = async () => {
  const file = fileInput.value.files[0];
  if (!file) {
    alert("📂 Select a FASTA file first");
    return;
  }

  const formData = new FormData();
  formData.append("file", file);
 
  try {
    result.value = null;
    await axios.post(`${API_URL}/pipeline/pdc`, formData, {
      headers: { "Content-Type": "multipart/form-data" }
    });
  } catch (error) {
    console.error("❌ Error uploading the file:", error);
  }

  startLogStream(); // ✅ Start log stream BEFORE sending request

};

const sendDirectSequence = async () => {
  if (!directSequence.value.trim()) {
    alert("Enter a valid FASTA sequence");
    return;
  }

  startLogStream(); // ✅ Start log stream BEFORE sending request

  try {
    await axios.post(`${API_URL}/pipeline/pdc`, {
      direct_sequence: directSequence.value
    });
  } catch (error) {
    console.error("❌ Error sending the sequence:", error);
  }
};

onUnmounted(() => stopLogStream());
</script>

<template>
  <div class="container">
    <h1>🔬 PDC Analyzer</h1>

    <div class="upload-section">
      <label for="file">Upload FASTA File:</label>
      <input type="file" ref="fileInput" accept=".fasta" />
      <button @click="uploadFile">📤 Submit File</button>
    </div>

    <div class="text-input-section">
      <label for="sequence">Enter FASTA Sequence:</label>
      <textarea v-model="directSequence" placeholder=">Example\nATGCATGC..." rows="4"></textarea>
      <button @click="sendDirectSequence">🚀 Analyze Sequence</button>
    </div>
    <div class="result-section" v-if="result">
      <h2>Result:</h2>
      <table>
        <tbody>
          <tr v-for="(value, key) in JSON.parse(result)">
            <th>{{ key }}</th>
            <td>{{ value }}</td>
          </tr>
        </tbody>
      </table>
    </div>
    <div class="progress-section" v-if="!result && logText.length > 0">
      <table>
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
      </table>
    </div>

    <div class="logs-section">
      <h2>Real-Time Logs:</h2>  
      <div ref="logContainer" class="logs">
        <pre>{{ logText }}</pre>
      </div>
    </div>
  </div>
</template>

<style scoped>
.progress-section table{
  margin:0 auto;
}
.result-section table{
  margin:0 auto;
}
.result-section table th{
  text-align:right;
  background: #00334c;
  color: white;
  font-weight: bold;
}
.result-section table td{
  background: #ffffff;
  text-align:left;
}
.result-section table th, .result-section table td{
  padding: 10px;
}
.result-section tbody tr:nth-child(even) {
  background: #f1f1f1;
}

.result-section table th{
  font-weight: bold;
}
.result-section table td{
  border-bottom: 1px solid #ddd;
}


.container {
  width: 900px;
  margin: auto;
  text-align: center;
}

.upload-section, .text-input-section {
  margin-bottom: 20px;
}

textarea {
  width: 100%;
  height: 100px;
}

.logs {
  width: 100%;
  height: 200px;
  background: black;
  color: green;
  font-family: monospace;
  overflow-y: auto;
  resize: none;
}
</style>

<script setup>
import { ref, computed, onUnmounted, nextTick } from "vue";
import axios from "axios";

const fileInput = ref(null);
const directSequence = ref("");
const logMessages = ref([]);
const eventSource = ref(null);
const isStreaming = ref(false);
const result = ref(null); 
const logContainer = ref(null);

const logText = computed(() => logMessages.value.join("\n"));

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

  eventSource.value = new EventSource("http://localhost:5000/pipeline/pdc/logs");

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
    await axios.post("http://localhost:5000/pipeline/pdc", formData, {
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
    await axios.post("http://localhost:5000/pipeline/pdc", {
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
      <!--{"sample_name": "PA01-DK_S26_L001.SPAdes.denovoassembly", "PDC": "", "PDC_REFERENCE": "PDC-1", "bit_score": "807.364", "gaps": "0", "identity": "100.0"}-->
      <table>
        <tr v-for="(value, key) in JSON.parse(result)">
          <th>{{ key }}</th>
          <td>{{ value }}</td>
        </tr>
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

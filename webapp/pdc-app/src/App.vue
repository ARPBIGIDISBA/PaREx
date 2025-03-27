<script setup>
import { ref } from 'vue'
import PDCAnalyzer from './components/PDCAnalyzer.vue'
import logo from './assets/pdcanalyzer.png'

const showIntro = ref(true)

const onLogoLoaded = () => {
  setTimeout(() => {
    showIntro.value = false
  }, 1000)
}
</script>

<template>
  <v-app>
    <!-- capa del logo centrat -->
    <transition name="fade" appear>
      <div v-if="showIntro" class="intro-logo">
        <img :src="logo" @load="onLogoLoaded" class="logo-big" />
      </div>
    </transition>

    <!-- layout principal sempre visible -->
    <v-layout class="rounded rounded-md border" style="height: 100vh;width: 100vw;">
      <v-app-bar color="primary" name="app-bar" title="PDC Analyzer" dark />
      <v-main>
        <v-container>
          <v-row>
            <v-col><PDCAnalyzer /></v-col>
          </v-row>
          <v-row>
            <v-col cols="3">
              <img :src="logo" alt="PDC Analyzer" height="150px" />
            </v-col>
            <v-col cols="9">
              <v-card elevation="10">
                <v-card-text>
                  Text for the citation of the article. Team ArpBigidisBA. (2024). PDC Analyzer. https://arpbigidisba.com
                  <br />
                  <a href="https://arpbigidisba.com/wp-content/uploads/2024/04/pseudomonas-aeruginosa-derived-cephalosporin-database-arpbig_190424.xlsx">
                    Download the database in Excel format
                  </a>
                </v-card-text>
              </v-card>
            </v-col>
          </v-row>
        </v-container>
      </v-main>
    </v-layout>
  </v-app>
</template>

<style scoped>
.intro-logo {
  position: fixed;
  top: 0;
  left: 0;
  z-index: 10000;
  width: 100vw;
  height: 100vh;
  background-color: white;
  display: flex;
  align-items: center;
  justify-content: center;
}

.logo-big {
  width: 1200px;
  height: auto;
  animation: zoomFade 0.5s ease-in-out;
}

/* Animació del logo gran */
@keyframes zoomFade {
  from {
    opacity: 0;
    transform: scale(1.5);
  }
  to {
    opacity: 1;
    transform: scale(1);
  }
}

/* Animació fade-out del div */
.fade-leave-active {
  transition: opacity 0.5s ease;
}
.fade-leave-to {
  opacity: 0;
}
</style>

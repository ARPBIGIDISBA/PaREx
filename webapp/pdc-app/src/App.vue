<script setup>
import PDCAnalyzer from './components/PDCAnalyzer.vue';
import logo from './assets/pdcanalyzer.png'
</script>

<template>
  <v-app>
    <v-layout class="rounded rounded-md border" style="height: 100vh;width: 100vw;">
      <v-app-bar color="primary" name="app-bar" title = "PDC Analyzer" dark>      
      </v-app-bar>
      <v-main style="margin-top:0px">
        <v-container>
          <v-row>
            <v-col><PDCAnalyzer /></v-col>
          </v-row>           
          <v-row>
            <v-col cols="3" ><img :src="logo" alt="PDC Analyzer" height="250px" /></v-col>
            <v-col cols="9">
              <v-card elevation-10>
                  <v-card-text>
                    Text for the citation of the article. Team ArpBigidisBA. (2024). PDC Analyzer. https://arpbigidisba.com
                    <br/>
                    <a href="https://arpbigidisba.com/wp-content/uploads/2024/04/pseudomonas-aeruginosa-derived-cephalosporin-database-arpbig_190424.xlsx" >
                      Download the database in excell format
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

<style>
</style>

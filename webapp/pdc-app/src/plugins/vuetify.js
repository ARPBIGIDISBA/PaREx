import 'vuetify/styles' 
import '@mdi/font/css/materialdesignicons.css'
import 'material-design-icons-iconfont/dist/material-design-icons.css'

import { createVuetify } from 'vuetify'
import * as components from "vuetify/components";
import * as labsComponents from 'vuetify/labs/components'
import * as directives from "vuetify/directives";
import { aliases, md } from 'vuetify/iconsets/md'
import { mdi } from 'vuetify/iconsets/mdi'

const themes_list = {
    pdc_theme: {
        dark: false,
        colors: {
          primary: "#00334c",     // Blau fosc del text
          secondary: "#2dc6d6",   // Turquesa de l’hèlix
          accent: "#2dc6d6",      // Mismo turquesa per accents
          background: "#ffffff",  // Fons blanc
          surface: "#f9f9f9",     // Ligerament gris per contrast
          success: "#28a745",
          info: "#2dc6d6",
          warning: "#ffc107",
          error: "#e53935",
          dark: "#212121",
          light: "#f4f4f4",
          anchor: "#2dc6d6"
        }
      }
}

const vuetify = createVuetify({
    theme: {
        defaultTheme: 'pdc_theme',
        themes: themes_list
    },
    icons: {
        defaultSet: 'md', // Establece el conjunto de iconos predeterminado a 'md'
        aliases,
        sets: {
            md,
            mdi
        }
    },
    components: {
        ...components,
        ...labsComponents
    },
    directives
})

export default vuetify


import react from '@vitejs/plugin-react';
import { defineConfig } from 'vite';

/* Num project site o app pode não viver na raiz; fica em variável de ambiente
   para o build local continuar servindo de "/" e outro host poder usar o seu. */
const base = process.env.BASE_PATH ?? '/';

export default defineConfig({
  base,
  plugins: [react()],
  /* Duas páginas: o aplicativo e o tutorial. O tutorial é entrada de verdade, e
     não arquivo em public/, porque de public/ o servidor de desenvolvimento cai
     no fallback de SPA e devolve o APP em /tutorial/ — o endereço só
     funcionaria depois de publicado, que é a pior hora de descobrir. */
  build: {
    rollupOptions: { input: { app: 'index.html', tutorial: 'tutorial/index.html' } },
  },
  worker: { format: 'es' },
  server: { port: 5173, open: false },
  test: {
    globals: true,
    include: ['src/**/*.test.ts', 'src/**/*.test.tsx'],
    environment: 'node',
  },
} as Parameters<typeof defineConfig>[0]);

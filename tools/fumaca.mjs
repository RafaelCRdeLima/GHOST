#!/usr/bin/env node
/* fumaca.mjs — o aplicativo abre, propaga e DESENHA?
 *
 * Os testes de fixture garantem que os números e o mapeamento estão certos.
 * Este garante que a figura existe: mede se o canvas tem pixels não-pretos.
 * Sem essa checagem, um mapa preto com séries corretas passa despercebido.
 *
 * ALCANCE, para ninguém confiar demais: ele verifica que a figura FINAL existe
 * e que a norma se conservou. Não verifica a animação quadro a quadro — se só
 * o desenho por quadro quebrasse, o mapa ficaria parado até o fim da
 * propagação e este teste ainda passaria. Quem garante o CONTEÚDO da figura
 * são as fixtures de src/nucleo/indices.test.ts, não este arquivo.
 */
import { spawn } from 'node:child_process';

const ALVO = process.argv[2] ?? 'http://localhost:4173/';

const EXPR = `(() => {
  const f = document.querySelector('footer')?.innerText ?? '';
  if (!/pronto|done|terminé|fatto/.test(f)) return false;
  /* O rodape mostra o desvio da norma como (Δ x): e o diagnostico permanente
     de CONVENTIONS.md 10.3, e o teste le exatamente o que o usuario le. */
  const m = f.match(/Δ ([0-9.eE+-]+)/);
  if (!m) return 'FALHA: rodape sem diagnostico de norma';
  const desvio = Number(m[1]);
  if (!(desvio < 1e-9)) return 'FALHA: |1-norma| = ' + m[1];
  const c = document.querySelector('canvas');
  if (!c) return 'FALHA: sem canvas';
  const gl = c.getContext('webgl2');
  if (!gl) return 'FALHA: sem contexto webgl2';
  /* CANVAS INTEIRO, e o criterio e CONTRASTE.
     cor(0) soma 16+26+36 = 78: um limiar de "pixels vivos" abaixo disso conta
     o proprio mapa desenhado todo no zero da rampa, mede 100%, e o guarda de
     "mapa quase preto" nao pode disparar no caso que ele nomeia.
     Mas exigir uma FRACAO alta acima de cor(0) tambem esta errado, e por
     motivo fisico: no primeiro quadro a excitacao esta num sitio so, e o mapa
     DEVE estar quase todo em cor(0); esse criterio reprovaria uma corrida
     legitima.
     O que distingue "desenhou" de "nao desenhou" nos dois casos e o contraste.
     E so o maximo serve como criterio: a FRACAO acima de cor(0) depende de qual
     quadro da animacao esta na tela quando o teste dispara (0,29% e 4,4% em
     duas corridas iguais). Um limiar sobre ela seria um teste que reprova por
     acaso, que e pior que um que nao reprova nunca. Ela continua no relatorio
     como informacao, sem ser porta.
     Medido na corrida padrao: min = 78 (exatamente cor(0)) e max = 677. Um mapa
     inteiro no zero da rampa daria max = 78 e reprova. */
  const L = gl.drawingBufferWidth, A = gl.drawingBufferHeight;
  const px = new Uint8Array(4 * L * A);
  gl.readPixels(0, 0, L, A, gl.RGBA, gl.UNSIGNED_BYTE, px);
  let vivos = 0, maxS = 0;
  for (let i = 0; i < px.length; i += 4) {
    const v = px[i] + px[i+1] + px[i+2];
    if (v > maxS) maxS = v;
    if (v > 80) ++vivos;
  }
  const frac = vivos / (L * A);
  if (maxS < 200) return 'FALHA: mapa sem contraste, maximo ' + maxS + ' contra cor(0)=78';
  return 'ok  |1-norma|=' + m[1] + '  acima de cor(0)=' + (100*frac).toFixed(2) + '%  max=' + maxS;
})()`;

/* O navegador considera "pronto" qualquer valor verdadeiro, e "FALHA: ..." é
   uma string verdadeira: o processo sairia com codigo 0 anunciando a falha. O
   veredito tem de ser lido, nao herdado. */
const p = spawn(process.execPath, [
  new URL('navegador.mjs', import.meta.url).pathname, ALVO, EXPR,
], { stdio: ['ignore', 'pipe', 'inherit'] });

let saida = '';
p.stdout.on('data', (b) => { saida += b; process.stdout.write(b); });
p.on('exit', (codigo) => {
  const ok = codigo === 0 && /PRONTO: "ok/.test(saida);
  if (!ok) console.error('  fumaca: FALHOU');
  process.exit(ok ? 0 : 1);
});

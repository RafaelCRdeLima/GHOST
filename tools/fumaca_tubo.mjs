#!/usr/bin/env node
/* fumaca_tubo.mjs — a vista do cilindro desenha, e GIRAR muda a figura.
 *
 * A segunda metade e que e o teste. "O tubo aparece" passa tambem com uma
 * projecao que ignore profundidade: uma faixa de pontos parece um tubo visto de
 * frente, e ninguem descobre olhando. Se a projecao fosse plana, inclinar o
 * eixo nao mudaria nada — entao a afirmacao "isto e tridimensional" e medida
 * como "arrastar reorganiza a figura", e nao como "ha um botao escrito 3D".
 *
 * Uso: node tools/fumaca_tubo.mjs <url>
 */
import { spawn } from 'node:child_process';
import { writeFileSync } from 'node:fs';

const URL_ALVO = process.argv[2];
const PNG = process.argv[3];
if (!URL_ALVO) { console.error('uso: fumaca_tubo.mjs <url> [saida.png]'); process.exit(2); }

const PORTA = 9335;
const chrome = spawn('google-chrome', [
  '--headless=new', '--disable-gpu', '--enable-unsafe-swiftshader',
  '--hide-scrollbars', '--window-size=1500,950', '--no-first-run',
  `--remote-debugging-port=${PORTA}`, '--user-data-dir=/tmp/daedalus-chrome-tubo',
  'about:blank',
], { stdio: 'ignore' });
const dorme = (ms) => new Promise((r) => setTimeout(r, ms));
const fim = (cod, msg) => { console.log(msg); chrome.kill(); process.exit(cod); };

let wsUrl;
for (let i = 0; i < 60 && !wsUrl; ++i) {
  try {
    const alvos = await (await fetch(`http://127.0.0.1:${PORTA}/json/list`)).json();
    const p = alvos.find((a) => a.type === 'page');
    if (p) wsUrl = p.webSocketDebuggerUrl;
  } catch { /* subindo */ }
  if (!wsUrl) await dorme(200);
}
if (!wsUrl) fim(1, 'FALHA: Chrome nao abriu a porta de depuracao');

const sock = new WebSocket(wsUrl);
await new Promise((r) => sock.addEventListener('open', r));
let id = 0; const pend = new Map();
sock.addEventListener('message', (e) => {
  const m = JSON.parse(e.data);
  if (m.id && pend.has(m.id)) { pend.get(m.id)(m); pend.delete(m.id); }
});
const cmd = (method, params = {}) => new Promise((r) => {
  const i = ++id; pend.set(i, r); sock.send(JSON.stringify({ id: i, method, params }));
});
const aval = async (expr) => {
  const r = await cmd('Runtime.evaluate', { expression: expr, returnByValue: true, awaitPromise: true });
  if (r.result?.exceptionDetails) throw new Error(r.result.exceptionDetails.text);
  return r.result?.result?.value;
};

await cmd('Page.enable'); await cmd('Runtime.enable');
await cmd('Page.navigate', { url: URL_ALVO });

/* 1. O aplicativo fica pronto na vista padrao. */
let pronto = false;
for (let i = 0; i < 120 && !pronto; ++i) {
  pronto = await aval(`/pronto|done|terminé|fatto/.test(document.querySelector('footer')?.innerText ?? '')`);
  if (!pronto) await dorme(500);
}
if (!pronto) fim(1, 'FALHA: o aplicativo nao ficou pronto');

/* 2. Trocar para a vista do tubo. */
const trocou = await aval(`(() => {
  const b = [...document.querySelectorAll('.vistas button')]
    .find((x) => /tub[oe]/i.test(x.textContent || ''));
  if (!b) return 'sem botao';
  b.click(); return 'ok';
})()`);
if (trocou !== 'ok') fim(1, `FALHA: ${trocou} — a vista do tubo nao esta oferecida`);
await dorme(1200);

/* Assinatura da figura: luminancia ACIMA DO FUNDO, media em blocos 8x8 do
   canvas inteiro.
   O fundo (--dd-canvas) soma 11+18+25 = 54, entao um limiar de 40 daria "100%
   de pixels vivos" medindo o proprio fundo. O limiar tem de ficar acima de 54,
   e o sinal do bloco tem de ser o EXCESSO sobre o fundo, senao a media de cada bloco e dominada pelo fundo e
   a assinatura para de distinguir figuras. */
const AMOSTRA = `(() => {
  const c = document.querySelector('.mapa canvas');
  if (!c) return { erro: 'sem canvas' };
  const gl = c.getContext('webgl2');
  if (!gl) return { erro: 'sem webgl2' };
  const L = gl.drawingBufferWidth, A = gl.drawingBufferHeight;
  const px = new Uint8Array(4 * L * A);
  gl.readPixels(0, 0, L, A, gl.RGBA, gl.UNSIGNED_BYTE, px);
  const FUNDO = 54;
  let vivos = 0;
  const bloco = new Float64Array(64), conta = new Float64Array(64);
  for (let y = 0; y < A; ++y) for (let x = 0; x < L; ++x) {
    const i = 4 * (y * L + x);
    const acima = Math.max(0, px[i] + px[i + 1] + px[i + 2] - FUNDO);
    if (acima > 36) ++vivos;
    const b = Math.min(7, (y * 8 / A) | 0) * 8 + Math.min(7, (x * 8 / L) | 0);
    bloco[b] += acima; conta[b] += 1;
  }
  return {
    vivos: vivos / (L * A),
    bloco: Array.from(bloco, (v, i) => v / Math.max(1, conta[i])),
  };
})()`;

const captura = async (sufixo) => {
  if (!PNG) return;
  const r = await cmd('Page.captureScreenshot', { format: 'png' });
  const arq = PNG.replace(/\.png$/, `${sufixo}.png`);
  writeFileSync(arq, Buffer.from(r.result.data, 'base64'));
  console.error(`  captura -> ${arq}`);
};

const antes = await aval(AMOSTRA);
await captura('');
if (antes?.erro) fim(1, `FALHA: ${antes.erro}`);
if (!(antes.vivos > 0.004)) fim(1, `FALHA: o tubo nao desenhou (${(100 * antes.vivos).toFixed(2)}% de pixels acima do fundo)`);

/* 3. Arrastar na horizontal: aponta o EIXO do tubo para dentro da tela. Numa
      projecao plana isso nao teria efeito nenhum. */
await aval(`(() => {
  const c = document.querySelector('.mapa canvas');
  const r = c.getBoundingClientRect();
  const x0 = r.left + r.width / 2, y0 = r.top + r.height / 2;
  const ev = (tipo, x, y) => c.dispatchEvent(new PointerEvent(tipo, {
    bubbles: true, cancelable: true, pointerId: 1, isPrimary: true,
    clientX: x, clientY: y, buttons: 1,
  }));
  ev('pointerdown', x0, y0);
  for (let k = 1; k <= 24; ++k) ev('pointermove', x0 + k * 5, y0);
  ev('pointerup', x0 + 120, y0);
  return true;
})()`);
await dorme(600);

await captura('-girado');

const depois = await aval(AMOSTRA);
if (depois?.erro) fim(1, `FALHA: ${depois.erro}`);

let mudou = 0, pior = 0;
for (let i = 0; i < 64; ++i) {
  const d = Math.abs(antes.bloco[i] - depois.bloco[i]);
  if (d > pior) pior = d;
  if (d > 1) ++mudou;
}

const linha = `acima do fundo ${(100 * antes.vivos).toFixed(2)}% -> ${(100 * depois.vivos).toFixed(2)}%` +
              `  blocos alterados pelo giro ${mudou}/64  (maior mudanca ${pior.toFixed(1)})`;
if (mudou < 16) fim(1, `FALHA: girar quase nao mudou a figura (${linha}) — projecao sem profundidade?`);
fim(0, `PRONTO: ok  ${linha}`);

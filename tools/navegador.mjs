#!/usr/bin/env node
/* navegador.mjs — teste de fumaça da interface, sem dependências.
 *
 * Dirige o Chrome pelo protocolo DevTools usando o WebSocket embutido do Node.
 * Existe porque `--dump-dom` e `--virtual-time-budget` não são confiáveis para
 * páginas que fazem trabalho assíncrono em Web Worker: dão silêncio quando o
 * worker ainda não respondeu, e silêncio é indistinguível de falha.
 *
 * Uso: node tools/navegador.mjs <url> <expressão-de-prontidão> [saida.png]
 *   A expressão é avaliada na página até devolver algo verdadeiro.
 */
import { spawn } from 'node:child_process';
import { writeFileSync } from 'node:fs';

const [, , url, prontidao, saida] = process.argv;
if (!url) { console.error('uso: navegador.mjs <url> <expr> [png]'); process.exit(2); }

const PORTA = 9333;
const chrome = spawn('google-chrome', [
  '--headless=new', '--disable-gpu', '--enable-unsafe-swiftshader',
  '--hide-scrollbars', '--window-size=1500,950', '--no-first-run',
  `--remote-debugging-port=${PORTA}`, '--user-data-dir=/tmp/daedalus-chrome',
  'about:blank',
], { stdio: 'ignore' });

const dorme = (ms) => new Promise((r) => setTimeout(r, ms));

async function alvo() {
  for (let i = 0; i < 60; ++i) {
    try {
      const r = await fetch(`http://127.0.0.1:${PORTA}/json/list`);
      const alvos = await r.json();
      const p = alvos.find((a) => a.type === 'page');
      if (p) return p.webSocketDebuggerUrl;
    } catch { /* ainda subindo */ }
    await dorme(200);
  }
  throw new Error('Chrome nao abriu a porta de depuracao');
}

const ws = new WebSocket(await alvo());
await new Promise((r) => { ws.onopen = r; });

let id = 0;
const pendentes = new Map();
const console_ = [];
ws.onmessage = (ev) => {
  const m = JSON.parse(ev.data);
  if (m.id && pendentes.has(m.id)) { pendentes.get(m.id)(m); pendentes.delete(m.id); }
  if (m.method === 'Runtime.consoleAPICalled') {
    console_.push(m.params.args.map((a) => a.value ?? a.description ?? '').join(' '));
  }
  if (m.method === 'Target.attachedToTarget') {
    ws.send(JSON.stringify({ id: ++id, sessionId: m.params.sessionId,
                             method: 'Runtime.enable', params: {} }));
  }
  if (m.method === 'Runtime.exceptionThrown') {
    console_.push('EXCECAO: ' + (m.params.exceptionDetails.exception?.description
                  ?? m.params.exceptionDetails.text));
  }
};
const cmd = (method, params = {}) => new Promise((res) => {
  const i = ++id;
  pendentes.set(i, res);
  ws.send(JSON.stringify({ id: i, method, params }));
});

await cmd('Runtime.enable');
await cmd('Page.enable');
await cmd('Log.enable');
/* Workers tem contexto de execucao proprio: sem auto-attach, o console deles
   nao aparece, e um worker travado fica indistinguivel de um worker ausente. */
await cmd('Target.setAutoAttach', { autoAttach: true, waitForDebuggerOnStart: false,
                                    flatten: true });
await cmd('Page.navigate', { url });

let pronto = false, valor = null;
for (let i = 0; i < 120; ++i) {           /* até 30 s de tempo REAL */
  await dorme(250);
  const r = await cmd('Runtime.evaluate', {
    expression: prontidao ?? 'true', returnByValue: true, awaitPromise: false,
  });
  valor = r.result?.result?.value;
  if (valor) { pronto = true; break; }
}

if (saida) {
  const s = await cmd('Page.captureScreenshot', { format: 'png' });
  if (s.result?.data) writeFileSync(saida, Buffer.from(s.result.data, 'base64'));
}

const fim = await cmd('Runtime.evaluate', {
  expression: "document.querySelector('footer')?.innerText ?? '(sem rodape)'",
  returnByValue: true,
});
for (const l of console_) console.log('  [pagina]', l);
console.log('  [rodape]', (fim.result?.result?.value ?? '').replace(/\n/g, ' | '));
console.log(pronto ? `  PRONTO: ${JSON.stringify(valor)}` : '  NAO FICOU PRONTO em 30 s');
ws.close();
chrome.kill();
process.exit(pronto ? 0 : 1);

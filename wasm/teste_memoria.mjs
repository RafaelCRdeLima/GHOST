/* teste_memoria.mjs — a armadilha da fronteira, demonstrada.
 *
 * Faz o heap do WASM crescer no meio do uso e verifica as DUAS metades:
 *
 *   1. que uma view guardada de fato morre — se ela sobrevivesse, este teste
 *      não provaria nada e o contrato seria superstição;
 *   2. que a view refeita a cada leitura continua certa.
 *
 * É o teste que pega o defeito "funciona até o usuário aumentar N".
 */
import { Daedalus } from './build/daedalus.mjs';

let falhas = 0;
const ok = (cond, texto) => {
  if (cond) console.log(`   . ${texto}`);
  else { console.log(`   x ${texto}`); ++falhas; }
};

const d = await Daedalus.criar();
console.log(`== fronteira de memoria  (daedalus ${d.versao}, nucleo ${d.hashNucleo})`);

const spec = (nPar, alvo) => ({
  format_version: 2, seed: 1,
  graph: { generator: 'microtubule',
           params: { n_par: nPar, n_perp: 13, seam_shift: 3, n_modules: 2 } },
  hamiltonian: { kind: 'adjacency', gamma: 1, normalization: 'none', lanczos_steps: 0 },
  initial: { site: 0 }, time: { t1: 5, nt: 20 },
  observables: { target: alvo, module_concurrence: false },
});

/* rede pequena, alguns passos, e uma leitura */
d.carregar(spec(20, 259));
d.avancar(5);

const bufferAntes = d.n_.wasmMemory.buffer;
const viewGuardada = d.populacao();          /* <- exatamente o que NAO se faz */
const somaAntes = viewGuardada.reduce((a, b) => a + b, 0);
/* O buffer de populacao e f32 DE PROJETO — ele alimenta a textura WebGL, e so
   a propagacao e f64 (CONVENTIONS.md, parte 9). Somar N valores f32 acumula
   ~N*6e-8; cobrar tolerancia de f64 aqui e erro de categoria no teste, nao
   defeito no nucleo. Quem vigia a norma em f64 e a serie temporal, abaixo. */
ok(Math.abs(somaAntes - 1) < 1e-5, `populacao (f32) soma 1 antes de crescer (${somaAntes})`);

/* agora o usuario aumenta N: o C aloca muito e o heap cresce */
d.carregar(spec(4000, -1));
const bufferDepois = d.n_.wasmMemory.buffer;

ok(bufferDepois !== bufferAntes,
   'o heap realmente cresceu (ArrayBuffer trocado) — sem isso o teste nao prova nada');
ok(viewGuardada.byteLength === 0 || viewGuardada.buffer.detached === true,
   `a view guardada morreu, como tinha de morrer (byteLength=${viewGuardada.byteLength})`);

/* a leitura refeita continua correta */
d.avancar(20);
const viewNova = d.populacao();
const somaDepois = viewNova.reduce((a, b) => a + b, 0);
ok(viewNova.length === 4000 * 13, `view refeita tem o tamanho novo (${viewNova.length})`);
ok(Math.abs(somaDepois - 1) < 1e-4,
   `populacao (f32) refeita soma 1 depois de crescer (${somaDepois})`);

/* e as demais leituras tambem se refazem */
const s = d.series();
const norma = s[(20 - 1) * 4 + 0];
ok(Math.abs(norma - 1) < 1e-12, `serie temporal legivel apos crescimento (norma=${norma})`);

d.destruir();
console.log(falhas === 0 ? `   PASSOU  6 verificacoes` : `   FALHOU  ${falhas} de 6`);
process.exit(falhas === 0 ? 0 : 1);

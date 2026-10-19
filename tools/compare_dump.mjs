#!/usr/bin/env node
/* compare_dump.mjs — confronta duas saidas de native/dump.c.
 *
 * Uso: node tools/compare_dump.mjs a.txt b.txt [tolerancia]
 *
 * Token nao-numerico tem de ser IDENTICO — inclusive a impressao digital
 * inteira do grafo, que e o primeiro lugar onde um gerador divergente
 * apareceria, antes de qualquer discussao sobre ponto flutuante.
 */
import { readFileSync } from 'node:fs';

const [, , fa, fb, tolArg, minExatoArg] = process.argv;
const TOL = tolArg ? Number(tolArg) : 1e-14;
/* Piso de identidade bit a bit. Os dois alvos dao 100%: o -msimd128 nao
   reassocia soma de ponto flutuante, entao o SpMV esparso soma na mesma ordem.
   Exigir o piso e vigiar o MECANISMO: uma queda para "dentro da tolerancia mas
   reassociado" e sinal de flag de compilacao nova, e aparece aqui antes de
   virar discrepancia fisica. */
const MIN_EXATO = minExatoArg ? Number(minExatoArg) : 0;

const toks = (p) => readFileSync(p, 'utf8').split('\n')
  .filter((l) => !l.startsWith('#')).join(' ').split(/\s+/).filter(Boolean);

const A = toks(fa), B = toks(fb);
if (A.length !== B.length) {
  console.error(`ESTRUTURA DIFERENTE: ${A.length} contra ${B.length} tokens`);
  process.exit(1);
}

let pior = 0, piorOnde = '', numericos = 0, exatos = 0;
for (let i = 0; i < A.length; ++i) {
  const na = Number(A[i]), nb = Number(B[i]);
  const ehNum = A[i] !== '' && !Number.isNaN(na) && !Number.isNaN(nb);
  if (!ehNum) {
    if (A[i] !== B[i]) {
      console.error(`TOKEN DIFERENTE em ${i}: "${A[i]}" contra "${B[i]}"`);
      process.exit(1);
    }
    continue;
  }
  ++numericos;
  if (A[i] === B[i]) ++exatos;
  const esc = Math.max(Math.abs(na), Math.abs(nb), 1e-300);
  const rel = Math.abs(na - nb) / esc;
  if (rel > pior) { pior = rel; piorOnde = `${i}: ${A[i]} contra ${B[i]}`; }
}

const pctExato = ((100 * exatos) / numericos).toFixed(2);
console.log(`  ${numericos} numeros comparados, ${pctExato}% bit a bit identicos`);
console.log(`  pior desvio relativo: ${pior.toExponential(3)}  (tolerancia ${TOL.toExponential(0)})`);
if (pior > TOL) {
  console.error(`  FALHOU em ${piorOnde}`);
  process.exit(1);
}
if (Number(pctExato) < MIN_EXATO) {
  console.error(`  FALHOU: identidade bit a bit caiu para ${pctExato}%, piso ${MIN_EXATO}%`);
  console.error('  (dentro da tolerancia, mas alguma flag passou a reassociar ponto flutuante)');
  process.exit(1);
}
console.log('  ok');

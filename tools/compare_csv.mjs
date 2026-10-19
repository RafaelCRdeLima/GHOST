#!/usr/bin/env node
/* compare_csv.mjs — confronta dois CSV do Daedalus.
 *
 * ORDEM DA COMPARAÇÃO IMPORTA. A impressão digital do grafo vem primeiro,
 * antes de qualquer observável: "o emissor serializou errado" e "o numérico
 * divergiu" têm causas e correções completamente diferentes, e sem a digital as
 * duas chegam como "os números não batem".
 *
 * ANTI-VACUIDADE. Dois arquivos idênticos de uma simulação degenerada passam em
 * qualquer comparação — norma se conserva trivialmente quando nada evolui, e um
 * delta parado bate consigo mesmo em todas as colunas. Por isso o comparador
 * exige que a corrida tenha ACONTECIDO: norma em 1, estado espalhado.
 *
 * Uso: compare_csv.mjs a.csv b.csv [tolerancia]
 */
import { readFileSync } from 'node:fs';

const [, , fa, fb, tolArg] = process.argv;
const TOL = tolArg ? Number(tolArg) : 1e-14;
let falhou = false;
const erro = (m) => { console.error(`  FALHOU: ${m}`); falhou = true; };

function ler(caminho) {
  const bruto = readFileSync(caminho, 'utf8').split('\n');
  const meta = new Map();
  const secoes = [];
  let atual = null;
  for (const linha of bruto) {
    if (linha.startsWith('#!')) {
      const t = linha.slice(2).trim();
      const i = t.indexOf(' ');
      meta.set(i < 0 ? t : t.slice(0, i), i < 0 ? '' : t.slice(i + 1));
      continue;
    }
    if (linha.startsWith('#') || linha.trim() === '') { atual = null; continue; }
    const campos = linha.split(',');
    if (!atual) { atual = { cab: campos, linhas: [] }; secoes.push(atual); continue; }
    atual.linhas.push(campos);
  }
  return { meta, secoes };
}

const A = ler(fa), B = ler(fb);

/* --- 1. impressão digital e identidade do núcleo, antes de tudo --- */
for (const chave of ['graph_fingerprint', 'core_hash', 'n', 'nnz', 'nmod']) {
  const a = A.meta.get(chave), b = B.meta.get(chave);
  if (a === undefined || b === undefined) continue;
  if (a !== b) {
    erro(`${chave} diverge: "${a}" contra "${b}"`);
    if (chave === 'graph_fingerprint') {
      console.error('  >> os dois lados construiram GRAFOS DIFERENTES.');
      console.error('  >> isso e o emissor ou o gerador, nao o propagador:');
      console.error('  >> compare os spec canonicos antes de olhar qualquer numero.');
    }
  }
}
if (falhou) process.exit(1);

/* --- 2. metadados numéricos compartilhados --- */
let piorMeta = 0, ondeMeta = '';
let comuns = 0;
for (const [k, va] of A.meta) {
  if (!B.meta.has(k) || k === 'spec') continue;
  const a = Number(va), b = Number(B.meta.get(k));
  if (Number.isNaN(a) || Number.isNaN(b)) {
    if (va !== B.meta.get(k)) erro(`meta ${k}: "${va}" contra "${B.meta.get(k)}"`);
    continue;
  }
  ++comuns;
  const rel = Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-300);
  if (rel > piorMeta) { piorMeta = rel; ondeMeta = k; }
}
if (A.meta.has('spec') && B.meta.has('spec') && A.meta.get('spec') !== B.meta.get('spec')) {
  erro('o spec.json canonico difere entre os dois lados');
}

/* --- 3. tabelas: só as colunas presentes nos dois --- */
let pior = 0, onde = '';
let comparados = 0;
for (let s = 0; s < Math.min(A.secoes.length, B.secoes.length); ++s) {
  const sa = A.secoes[s], sb = B.secoes[s];
  if (sa.linhas.length !== sb.linhas.length) {
    erro(`secao ${s}: ${sa.linhas.length} linhas contra ${sb.linhas.length}`);
    continue;
  }
  const pares = [];
  for (let i = 0; i < sa.cab.length; ++i) {
    const j = sb.cab.indexOf(sa.cab[i]);
    if (j >= 0) pares.push([sa.cab[i], i, j]);
  }
  if (pares.length === 0) erro(`secao ${s}: nenhuma coluna em comum`);
  for (let r = 0; r < sa.linhas.length; ++r) {
    for (const [nome, i, j] of pares) {
      const ta = sa.linhas[r][i], tb = sb.linhas[r][j];
      if (ta === undefined || tb === undefined) continue;
      const a = Number(ta), b = Number(tb);
      if (Number.isNaN(a) && Number.isNaN(b)) continue;   /* nan == nan aqui */
      if (Number.isNaN(a) !== Number.isNaN(b)) {
        erro(`${nome} linha ${r}: "${ta}" contra "${tb}" (um e nan)`);
        continue;
      }
      ++comparados;
      const rel = Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-300);
      if (rel > pior) { pior = rel; onde = `${nome} linha ${r}: ${ta} contra ${tb}`; }
    }
  }
}

/* --- 4. anti-vacuidade: a corrida aconteceu? --- */
function diagnostico(D, nome) {
  const tab = D.secoes[0];
  if (!tab || tab.linhas.length === 0) { erro(`${nome}: sem tabela`); return; }
  const iN = tab.cab.indexOf('norm'), iI = tab.cab.indexOf('ipr');
  const ult = tab.linhas[tab.linhas.length - 1];
  const norma = Number(ult[iN]), ipr = Number(ult[iI]);
  if (!(Math.abs(norma - 1) < 1e-9)) erro(`${nome}: norma final = ${norma}`);
  if (!(ipr < 0.5)) {
    erro(`${nome}: IPR final = ${ipr} — o pacote nao se espalhou, e uma`
       + ' comparacao entre duas simulacoes paradas passa sem testar nada');
  }
}
diagnostico(A, 'a');
diagnostico(B, 'b');

console.log(`  meta: ${comuns} campos comuns, pior desvio ${piorMeta.toExponential(2)}`
          + (ondeMeta ? ` (${ondeMeta})` : ''));
console.log(`  dados: ${comparados} numeros, pior desvio relativo ${pior.toExponential(3)}`
          + `  (tolerancia ${TOL.toExponential(0)})`);
if (pior > TOL) erro(onde);
if (piorMeta > TOL) erro(`metadado ${ondeMeta} fora da tolerancia`);
if (!falhou) console.log('  ok');
process.exit(falhou ? 1 : 0);

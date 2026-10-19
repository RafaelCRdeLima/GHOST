/* indices.ts — o mapeamento entre sítio e texel, num lugar só.
 *
 * A partir da interface, erro deixa de aparecer como número errado e passa a
 * aparecer como IMAGEM PLAUSÍVEL. Um heatmap com os eixos m/n trocados, ou com
 * um off-by-one na coluna, é indistinguível de física interessante para quem
 * está olhando — e vira figura de artigo.
 *
 * Grafo SEM rede desenrolada não entra aqui: ele vai para o layout espectral
 * (dae_layout_espectral). Não há aqui uma "tira por ordem de índice":
 * desenhar os sítios na ordem em que foram rotulados desenha a ROTULAGEM, não
 * a rede, e para SBM ou lista importada a ordem não significa nada.
 *
 * Por isso a transposição mora aqui, sozinha, exportada, e é o alvo direto de
 * src/nucleo/indices.test.ts: as fixtures de simetria e o instantâneo
 * congelado atacam ESTA função, não o pixel na tela.
 */

/** Índice do sítio (m, q) no vetor que o núcleo devolve. Espelha
 *  `j = m * N_perp + q` de dae_graph.c — ver CONVENTIONS.md, parte 5. */
export function indiceSitio(m: number, q: number, nPerp: number): number {
  return m * nPerp + q;
}

/** Empacota a população na ordem linha-maior que a textura WebGL espera:
 *  linha `q` (protofilamento), coluna `m` (posição ao longo do eixo).
 *
 *  É a única transposição do projeto. Se ela inverter, o heatmap continua
 *  bonito e passa a mentir. */
export function empacotarLattice(
  pop: Float32Array, nPar: number, nPerp: number, destino?: Float32Array,
): Float32Array {
  const out = destino ?? new Float32Array(nPar * nPerp);
  for (let q = 0; q < nPerp; ++q) {
    for (let m = 0; m < nPar; ++m) out[q * nPar + m] = pop[m * nPerp + q];
  }
  return out;
}

/** Do texel de volta ao sítio — usada pelo clique. Inversa de
 *  `empacotarLattice`, e testada como inversa, não por inspeção. */
export function sitioNoTexel(
  m: number, q: number, nPar: number, nPerp: number,
): number {
  if (m < 0 || q < 0 || m >= nPar || q >= nPerp) return -1;
  return indiceSitio(m, q, nPerp);
}

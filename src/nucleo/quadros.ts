/* quadros.ts — quais passos viram quadro de mapa.
 *
 * A série temporal tem um ponto por passo; o MAPA não precisa disso, e guardar
 * um quadro de N floats por passo estoura a memória depressa. Existe um passo
 * de amostragem (`pop_stride`), e a interface precisa CONHECÊ-LO: o sintoma de
 * não conhecer é a animação parar antes do fim, no ponto em que os quadros
 * acabam, com o deslizador ainda mostrando o total de passos.
 *
 * Duas regras que parecem detalhe e não são:
 *
 *   1. O ÚLTIMO passo entra sempre, seja ele múltiplo do passo de amostragem ou
 *      não. Sem isso, `nt = 400` com passo 30 pararia em 390 — a animação
 *      terminaria antes do estado final, que costuma ser o que se quer ver.
 *   2. O passo 0 (estado inicial) é o quadro zero. Sem ele não há de onde a
 *      animação partir.
 */

/** Passos em que um quadro é emitido, incluindo 0 e nt. */
export function passosDosQuadros(nt: number, passo: number): number[] {
  const p = Math.max(1, Math.floor(passo));
  const out = [0];
  for (let k = p; k < nt; k += p) out.push(k);
  if (nt > 0) out.push(nt);
  return out;
}

/** Bytes que os quadros ocupam na interface (f32 por sítio). */
export function bytesDosQuadros(nt: number, passo: number, n: number): number {
  return passosDosQuadros(nt, passo).length * n * 4;
}

/** Passo de amostragem mínimo para caber num teto de bytes. */
export function passoParaCaber(nt: number, n: number, teto: number): number {
  if (n <= 0 || nt <= 0) return 1;
  const maxQuadros = Math.max(2, Math.floor(teto / (n * 4)));
  return Math.max(1, Math.ceil(nt / (maxQuadros - 1)));
}

/** Índice do quadro cujo passo é o mais próximo de `passo`. */
export function quadroDoPasso(passos: readonly number[], passo: number): number {
  let melhor = 0, dist = Infinity;
  for (let i = 0; i < passos.length; ++i) {
    const d = Math.abs(passos[i] - passo);
    if (d < dist) { dist = d; melhor = i; }
  }
  return melhor;
}

/* paleta.ts — as escalas de cor do Daedalus, e as únicas.
 *
 * São os mesmos valores que o matplotlib e o Wolfram usam nos mapas
 * exportados, para que a figura da tela e a figura do artigo tenham a mesma
 * leitura.
 *
 * O shader monta a LUT a partir DESTA função, então o teste de monotonicidade
 * em paleta.test.ts vale para o que aparece na tela. Um colormap invertido é
 * um dos dois jeitos de o heatmap ficar bonito e mentir.
 */

/* Probabilidade: rampa sequencial com luminância monotônica, para sobreviver à
   impressão em escala de cinza. */
const PROB_ESCURO = [
  [0x10, 0x1a, 0x24], [0x1c, 0x3a, 0x56], [0x17, 0x55, 0x7c], [0x3f, 0x7c, 0x74],
  [0x9c, 0x9a, 0x3e], [0xe5, 0xa8, 0x3f], [0xfc, 0xe9, 0xc0],
] as const;

/* Fase: rampa CÍCLICA — 0 e 2π fecham na mesma cor, sem a descontinuidade
   artificial que um mapa sequencial criaria num observável periódico. */
const FASE = [
  [0x17, 0x55, 0x7c], [0x2f, 0x7d, 0x6a], [0x7e, 0x9a, 0x3c], [0xd2, 0xa0, 0x3a],
  [0xc0, 0x63, 0x2b], [0xa8, 0x45, 0x2c], [0x7c, 0x3c, 0x68], [0x45, 0x40, 0x8f],
  [0x17, 0x55, 0x7c],
] as const;

function amostrar(paradas: ReadonlyArray<readonly number[]>, t: number): [number, number, number] {
  const u = t <= 0 ? 0 : t >= 1 ? 1 : t;
  const x = u * (paradas.length - 1);
  const i = Math.min(paradas.length - 2, Math.floor(x));
  const f = x - i;
  const a = paradas[i], b = paradas[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

/** Probabilidade, t em [0,1] → [r,g,b] em 0..255. */
export function cor(t: number): [number, number, number] { return amostrar(PROB_ESCURO, t); }

/** Fase, t em [0,1] representando 0..2π. */
export function corFase(t: number): [number, number, number] { return amostrar(FASE, t); }

export function luminancia([r, g, b]: [number, number, number]): number {
  return (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255;
}

/** LUT 256x1 RGB para o shader. Mesma função, mesma escala. */
export function lut(n = 256): Uint8Array {
  const out = new Uint8Array(n * 3);
  for (let i = 0; i < n; ++i) {
    const [r, g, b] = cor(i / (n - 1));
    out[3 * i] = Math.round(r);
    out[3 * i + 1] = Math.round(g);
    out[3 * i + 2] = Math.round(b);
  }
  return out;
}

/** Gradiente CSS da barra de escala, a partir das mesmas paradas. */
export function gradienteCss(): string {
  const p = PROB_ESCURO.map(([r, g, b], i) =>
    `rgb(${r},${g},${b}) ${((100 * i) / (PROB_ESCURO.length - 1)).toFixed(1)}%`);
  return `linear-gradient(90deg, ${p.join(', ')})`;
}

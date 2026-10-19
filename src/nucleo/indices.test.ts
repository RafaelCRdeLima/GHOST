/* indices.test.ts — as fixtures visuais com resposta conhecida.
 *
 * A partir da interface, erro vira imagem plausivel. Estes casos tem resposta
 * que da para conferir a olho UMA vez e depois congelar: um delta tem de cair
 * num texel especifico, e a rede SEM costura tem de produzir padrao simetrico
 * em torno do sitio de origem. Assimetria ali e indexacao, nao fisica.
 */
import { beforeAll, describe, expect, it } from 'vitest';
import { assinatura } from './assinatura';
import { empacotarLattice, indiceSitio, sitioNoTexel } from './indices';

import { Daedalus } from '../../wasm/build/daedalus.mjs';

const N_PAR = 41, N_PERP = 13, M0 = 20, Q0 = 6;

let d: any;
beforeAll(async () => { d = await Daedalus.criar(); });

/* `passos = 0` devolve o estado inicial sem propagar — a grade ainda precisa
   de pelo menos um ponto, entao nt = 1 e nao se avanca. */
function propagar(seamShift: number, t1: number, passos: number): Float32Array {
  /* O mesmo spec.json que o .cpp exportado embute: a fixture exercita o caminho
     de verdade, nao um atalho de teste. */
  d.carregar({
    format_version: 2, seed: 1,
    graph: { generator: 'microtubule',
             params: { n_par: N_PAR, n_perp: N_PERP, seam_shift: seamShift, n_modules: 1 } },
    hamiltonian: { kind: 'adjacency', gamma: 1, normalization: 'none', lanczos_steps: 0 },
    initial: { site: indiceSitio(M0, Q0, N_PERP) },
    time: { t1, nt: Math.max(1, passos) },
    observables: { target: -1, module_concurrence: false },
  });
  if (passos > 0) d.avancar(passos);
  return empacotarLattice(d.populacao(), N_PAR, N_PERP);
}

describe('mapeamento sitio -> texel', () => {
  it('um delta cai no texel certo — pega transposicao na hora', () => {
    const tex = propagar(0, 1, 0);
    let iMax = -1, vMax = -1;
    for (let i = 0; i < tex.length; ++i) if (tex[i] > vMax) { vMax = tex[i]; iMax = i; }
    expect(vMax).toBeCloseTo(1, 12);
    expect(iMax).toBe(Q0 * N_PAR + M0);
    /* com N_par != N_perp, uma transposicao mudaria este indice */
    expect(N_PAR).not.toBe(N_PERP);
  });

  it('sitioNoTexel e a inversa de empacotarLattice', () => {
    for (const [m, q] of [[0, 0], [1, 0], [0, 1], [40, 12], [20, 6]]) {
      const j = sitioNoTexel(m, q, N_PAR, N_PERP);
      const pop = new Float32Array(N_PAR * N_PERP);
      pop[j] = 1;
      const tex = empacotarLattice(pop, N_PAR, N_PERP);
      expect(tex[q * N_PAR + m]).toBe(1);
    }
    expect(sitioNoTexel(-1, 0, N_PAR, N_PERP)).toBe(-1);
    expect(sitioNoTexel(0, N_PERP, N_PAR, N_PERP)).toBe(-1);
  });
});

describe('fixture de simetria: rede SEM costura', () => {
  /* seam_shift = 0 faz o grafo ser o produto cartesiano P_41 x C_13, entao a
     dinamica fatoriza e o padrao e simetrico nos dois eixos. O pacote nao
     alcanca as bordas em t = 3 (contribuicao ~ J_20(6)^2 ~ 1e-20). */
  let tex: Float32Array;
  beforeAll(() => { tex = propagar(0, 3, 30); });

  /* ANTI-VACUIDADE: com o estado parado — por exemplo, se o wrapper JS zerasse
     j_par e j_perp — TODA verificacao de simetria passa trivialmente, porque
     um delta e simetrico. Um heatmap
     parado tambem nao parece bug na tela: parece localizacao. */
  it('o pacote de fato se espalhou', () => {
    let soma = 0, ipr = 0, maxPop = 0;
    for (let i = 0; i < tex.length; ++i) {
      soma += tex[i]; ipr += tex[i] * tex[i];
      maxPop = Math.max(maxPop, tex[i]);
    }
    expect(soma).toBeCloseTo(1, 6);
    expect(maxPop).toBeLessThan(0.5);   /* nao ficou no sitio de origem */
    expect(ipr).toBeLessThan(0.1);      /* delocalizado                  */
    expect(tex[Q0 * N_PAR + M0]).toBeLessThan(0.2);
  });

  it('simetrico no eixo longitudinal em torno de m0', () => {
    let pior = 0;
    for (let q = 0; q < N_PERP; ++q)
      for (let dm = 1; dm <= 15; ++dm)
        pior = Math.max(pior, Math.abs(tex[q * N_PAR + M0 + dm] - tex[q * N_PAR + M0 - dm]));
    expect(pior).toBeLessThan(1e-12);
  });

  it('simetrico no eixo transversal em torno de q0', () => {
    let pior = 0;
    for (let m = 0; m < N_PAR; ++m)
      for (let dq = 1; dq <= 6; ++dq) {
        const a = tex[(((Q0 + dq) % N_PERP) * N_PAR) + m];
        const b = tex[((((Q0 - dq + N_PERP) % N_PERP)) * N_PAR) + m];
        pior = Math.max(pior, Math.abs(a - b));
      }
    expect(pior).toBeLessThan(1e-12);
  });
});

describe('fixture de assimetria: a costura quebra a simetria transversal', () => {
  /* Sem este caso o teste acima seria vacuo: um empacotamento que zerasse tudo
     tambem passaria por "simetrico". */
  it('seam_shift = 3 quebra o eixo transversal e preserva o longitudinal', () => {
    const tex = propagar(3, 3, 30);
    let piorTransversal = 0;
    for (let m = 0; m < N_PAR; ++m)
      for (let dq = 1; dq <= 6; ++dq) {
        const a = tex[(((Q0 + dq) % N_PERP) * N_PAR) + m];
        const b = tex[((((Q0 - dq + N_PERP) % N_PERP)) * N_PAR) + m];
        piorTransversal = Math.max(piorTransversal, Math.abs(a - b));
      }
    expect(piorTransversal).toBeGreaterThan(1e-6);

    let soma = 0;
    for (let i = 0; i < tex.length; ++i) soma += tex[i];
    expect(soma).toBeCloseTo(1, 6);
  });
});

describe('assinatura da figura', () => {
  /* Quase tudo aqui e PREVISAO ANALITICA, nao instantaneo medido — o unico
     numero congelado esta marcado. Ver o cabecalho de assinatura.ts para
     por que a estatistica marginal em q nao serve (ela e cega a costura). */
  const T = 3, GAMA = 1;

  it('sem costura: a rede e o produto P_41 x C_13 e o eixo longitudinal e livre', () => {
    const a = assinatura(propagar(0, T, 30), N_PAR, N_PERP);
    expect(a.soma).toBeCloseTo(1, 6);
    /* centroide na origem, por simetria */
    expect(a.centroM).toBeCloseTo(M0, 9);
    /* variancia da caminhada quantica livre na linha: 2 (gama t)^2 */
    expect(a.varM).toBeCloseTo(2 * (GAMA * T) ** 2, 5);
    /* cada anel ve a MESMA dinamica longitudinal: o produto cartesiano fatoriza */
    for (let q = 0; q < N_PERP; ++q) expect(a.centroPorLinha[q]).toBeCloseTo(M0, 8);
    /* massa igual em todos os aneis, por simetria do C_13 em torno de q0 */
    for (let q = 0; q < N_PERP; ++q)
      expect(a.massaPorLinha[q]).toBeCloseTo(a.massaPorLinha[(N_PERP + 2 * Q0 - q) % N_PERP], 9);
  });

  it('com costura: o perfil longitudinal por anel fica antissimetrico', () => {
    const a = assinatura(propagar(3, T, 30), N_PAR, N_PERP);
    expect(a.soma).toBeCloseTo(1, 6);
    expect(a.centroM).toBeCloseTo(M0, 9);          /* globalmente ainda simetrico */

    /* A costura liga (m, 12) a (m+3, 0): ela transporta massa em m ao ser
       cruzada, nos dois sentidos, e o desvio de cada anel e antissimetrico em
       torno do anel de origem. Previsao estrutural, nao numero medido. */
    for (let dq = 1; dq <= 6; ++dq) {
      const acima = a.centroPorLinha[(Q0 + dq) % N_PERP] - M0;
      const abaixo = a.centroPorLinha[(N_PERP + Q0 - dq) % N_PERP] - M0;
      expect(acima).toBeCloseTo(-abaixo, 8);
    }
    /* o anel mais distante e o que mais desloca — unico numero congelado */
    expect(a.centroPorLinha[0] - M0).toBeCloseTo(0.65233, 4);
    /* e o eixo longitudinal ganha espalhamento extra: a costura e um atalho */
    expect(a.varM).toBeGreaterThan(2 * (GAMA * T) ** 2 + 0.3);
  });
});

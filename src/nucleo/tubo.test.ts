/* tubo.test.ts — o embutimento helicoidal.
 *
 * A afirmação que interessa é geométrica e tem consequência física: com
 * z = m + seam·q/nq toda ligação lateral tem o mesmo comprimento, INCLUSIVE a
 * da costura. É isso que faz a costura ser descontinuidade de rede e não de
 * forma, e é isso que o desenho tem de mostrar.
 *
 * Cada afirmação vem com a companheira que prova que ela poderia quebrar: sem
 * o empilhamento ingênuo do lado, "as ligações têm o mesmo comprimento" seria
 * verdade por construção do teste, não por propriedade do embutimento.
 */
import { describe, expect, it } from 'vitest';
import { compressaoParaCaber, coordenadasTubo, passoDaHelice } from './tubo.ts';

const NP = 20, NQ = 13, SEAM = 3;

const dist = (a: Float32Array, i: number, j: number) =>
  Math.hypot(a[3 * i] - a[3 * j], a[3 * i + 1] - a[3 * j + 1], a[3 * i + 2] - a[3 * j + 2]);

/** O empilhamento ingênuo: anéis planos, z = m, a costura ignorada. */
function coordenadasPlanas(np: number, nq: number): Float32Array {
  const xyz = new Float32Array(3 * np * nq);
  const R = nq / (2 * Math.PI);
  for (let m = 0; m < np; ++m)
    for (let q = 0; q < nq; ++q) {
      const j = m * nq + q, th = (2 * Math.PI * q) / nq;
      xyz[3 * j] = m; xyz[3 * j + 1] = R * Math.cos(th); xyz[3 * j + 2] = R * Math.sin(th);
    }
  return xyz;
}

/** Comprimentos de todas as ligações laterais, a da costura por último. */
function lateraisEDaCostura(xyz: Float32Array, np: number, nq: number, seam: number) {
  const comuns: number[] = [];
  const costuras: number[] = [];
  for (let m = 0; m < np; ++m)
    for (let q = 0; q < nq; ++q) {
      const j = m * nq + q;
      if (q + 1 < nq) comuns.push(dist(xyz, j, m * nq + q + 1));
      else if (m + seam < np && m + seam >= 0) costuras.push(dist(xyz, j, (m + seam) * nq));
    }
  return { comuns, costuras };
}

describe('cilindro do microtubulo', () => {
  it('todo vertice existe, tres coordenadas cada', () => {
    expect(coordenadasTubo(NP, NQ, SEAM).length).toBe(3 * NP * NQ);
  });

  it('a costura tem o MESMO comprimento das outras ligacoes laterais', () => {
    const xyz = coordenadasTubo(NP, NQ, SEAM);
    const { comuns, costuras } = lateraisEDaCostura(xyz, NP, NQ, SEAM);
    expect(comuns.length).toBeGreaterThan(100);
    expect(costuras.length).toBeGreaterThan(10);   /* anti-vacuidade: ha costura */
    const ref = comuns[0];
    for (const d of comuns) expect(d).toBeCloseTo(ref, 5);
    for (const d of costuras) expect(d).toBeCloseTo(ref, 5);
  });

  /* COMPANHEIRA: sem ela o teste acima passaria com qualquer embutimento que
     tratasse a costura como as demais, inclusive um que estivesse errado. */
  it('e o empilhamento ingenuo NAO tem — a costura sai esticada', () => {
    const xyz = coordenadasPlanas(NP, NQ);
    const { comuns, costuras } = lateraisEDaCostura(xyz, NP, NQ, SEAM);
    const ref = comuns[0];
    for (const d of costuras) expect(d).toBeGreaterThan(ref * 2);
    /* e o excesso e exatamente o salto de `seam` unidades no eixo */
    expect(Math.hypot(SEAM, ref)).toBeCloseTo(costuras[0], 5);
  });

  it('ligacoes longitudinais valem 1, a unidade do eixo', () => {
    const xyz = coordenadasTubo(NP, NQ, SEAM);
    for (let m = 0; m + 1 < NP; ++m)
      for (let q = 0; q < NQ; ++q)
        expect(dist(xyz, m * NQ + q, (m + 1) * NQ + q)).toBeCloseTo(1, 5);
  });

  it('o raio poe o passo lateral na mesma escala do longitudinal', () => {
    const xyz = coordenadasTubo(NP, NQ, 0);
    /* corda entre vizinhos laterais ~ 1: sem isso o tubo sai espremido ou
       inflado, e a razao entre j_par e j_perp apareceria distorcida */
    expect(dist(xyz, 0, 1)).toBeGreaterThan(0.95);
    expect(dist(xyz, 0, 1)).toBeLessThan(1.0);
  });

  it('sem costura os aneis sao planos, com costura nao', () => {
    const plano = coordenadasTubo(NP, NQ, 0);
    for (let q = 0; q < NQ; ++q) expect(plano[3 * q]).toBeCloseTo(0, 6);
    const helice = coordenadasTubo(NP, NQ, SEAM);
    const eixos = Array.from({ length: NQ }, (_, q) => helice[3 * q]);
    expect(Math.max(...eixos) - Math.min(...eixos)).toBeCloseTo(SEAM * (NQ - 1) / NQ, 5);
  });

  it('o passo da helice e seam por volta completa', () => {
    expect(passoDaHelice(NQ, SEAM) * NQ).toBeCloseTo(SEAM, 10);
    expect(passoDaHelice(NQ, 0)).toBe(0);
  });

  /* A compressão do eixo é distorção ASSUMIDA — ela quebra a comparação entre
     distância longitudinal e lateral, e é por isso que o fator vai para a tela.
     O que ela não pode quebrar é a igualdade das laterais entre si: se
     quebrasse, a costura voltaria a aparecer como cicatriz e a vista passaria a
     afirmar uma coisa falsa sobre a rede. */
  it('com o eixo comprimido a costura CONTINUA igual as outras laterais', () => {
    const xyz = coordenadasTubo(NP, NQ, SEAM, 0.1);
    const { comuns, costuras } = lateraisEDaCostura(xyz, NP, NQ, SEAM);
    const ref = comuns[0];
    for (const d of comuns) expect(d).toBeCloseTo(ref, 5);
    for (const d of costuras) expect(d).toBeCloseTo(ref, 5);
    /* e ela de fato encurtou o eixo — senao o teste acima passaria por nada */
    const metrico = coordenadasTubo(NP, NQ, SEAM, 1);
    expect(xyz[3 * (NP - 1) * NQ]).toBeCloseTo(metrico[3 * (NP - 1) * NQ] * 0.1, 5);
  });

  it('a compressao mira a proporcao pedida e nunca estica', () => {
    /* 160 x 13 e o padrao do aplicativo: proporcao 39:1, um fio */
    const c = compressaoParaCaber(160, 13, 3.5);
    expect(c).toBeLessThan(0.2);
    const xyz = coordenadasTubo(160, 13, 3, c);
    let x0 = Infinity, x1 = -Infinity;
    for (let i = 0; i < 160 * 13; ++i) { const x = xyz[3 * i]; if (x < x0) x0 = x; if (x > x1) x1 = x; }
    expect((x1 - x0) / (13 / Math.PI)).toBeGreaterThan(3);
    expect((x1 - x0) / (13 / Math.PI)).toBeLessThan(5);
    /* tubo curto: escala metrica, sem esticar */
    expect(compressaoParaCaber(6, 13, 3.5)).toBe(1);
  });
});

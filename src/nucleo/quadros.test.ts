/* quadros.test.ts — a animação chega ao fim?
 *
 * "A animação chega ao fim" passa TRIVIALMENTE quando existe um quadro só: o
 * primeiro é o último. Por isso toda asserção de completude aqui vem com a
 * companheira que exige mais de um quadro, e há um caso degenerado explícito
 * onde a completude é verdadeira e vazia ao mesmo tempo.
 */
import { describe, expect, it } from 'vitest';
import { bytesDosQuadros, passoParaCaber, passosDosQuadros, quadroDoPasso } from './quadros.ts';

describe('plano de quadros', () => {
  it('o ultimo passo entra mesmo sem ser multiplo do passo de amostragem', () => {
    /* nt=400, passo=30: sem esta regra pararia em 390, e a animacao terminaria
       antes do estado final — que costuma ser justamente o que se quer ver. */
    const p = passosDosQuadros(400, 30);
    expect(p[p.length - 1]).toBe(400);
    expect(p.length).toBeGreaterThan(2);          /* companheira: nao e trivial */
    expect(p[0]).toBe(0);
  });

  it('conta certo em varios passos de amostragem', () => {
    for (const [nt, passo, esperado] of [[400, 1, 401], [400, 10, 41], [400, 30, 15],
                                         [100, 7, 16], [50, 50, 2]] as const) {
      const p = passosDosQuadros(nt, passo);
      expect(p.length, `nt=${nt} passo=${passo}`).toBe(esperado);
      expect(p[p.length - 1], `fim nt=${nt} passo=${passo}`).toBe(nt);
    }
  });

  it('CASO DEGENERADO: passo maior que nt da um quadro alem do inicial', () => {
    /* Aqui "chega ao fim" e verdadeiro E vazio: o unico quadro ja e o ultimo.
       O caso existe para deixar claro que a assercao de completude, sozinha,
       nao prova animacao nenhuma. */
    const p = passosDosQuadros(400, 100000);
    expect(p).toEqual([0, 400]);
    expect(p[p.length - 1]).toBe(400);            /* passa... */
    expect(p.length).toBe(2);                     /* ...e nao ha o que animar */
  });

  it('o passo nunca some: valores invalidos caem em 1', () => {
    expect(passosDosQuadros(10, 0).length).toBe(11);
    expect(passosDosQuadros(10, -5).length).toBe(11);
  });

  it('memoria cresce como o numero de quadros, e o passo a controla', () => {
    const n = 2080;
    const cheio = bytesDosQuadros(400, 1, n);
    const ralo = bytesDosQuadros(400, 10, n);
    expect(cheio).toBe(401 * n * 4);
    expect(ralo).toBeLessThan(cheio / 9);
  });

  it('passoParaCaber respeita o teto e nunca devolve menos que 1', () => {
    const n = 2080, teto = 6_000_000;
    const p = passoParaCaber(400, n, teto);
    expect(bytesDosQuadros(400, p, n)).toBeLessThanOrEqual(teto * 1.01);
    expect(passoParaCaber(10, 10, 1e12)).toBe(1);
  });

  it('quadroDoPasso mapeia passo -> quadro mais proximo', () => {
    const p = passosDosQuadros(400, 10);          /* 0, 10, 20, ... 400 */
    expect(p[quadroDoPasso(p, 0)]).toBe(0);
    expect(p[quadroDoPasso(p, 137)]).toBe(140);
    expect(p[quadroDoPasso(p, 400)]).toBe(400);
    /* companheira: um mapeamento que devolvesse sempre 0 passaria no primeiro
       caso e falha aqui */
    expect(quadroDoPasso(p, 400)).toBe(p.length - 1);
  });
});

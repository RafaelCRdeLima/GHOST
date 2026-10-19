import { describe, expect, it } from 'vitest';
import { cor, corFase, gradienteCss, luminancia, lut } from './paleta';

describe('paleta', () => {
  it('e monotonica em luminosidade — "mais claro = mais populacao" tem de ser verdade', () => {
    let anterior = -1;
    let quedas = 0;
    for (let i = 0; i <= 512; ++i) {
      const L = luminancia(cor(i / 512));
      if (L < anterior - 1e-9) ++quedas;
      anterior = L;
    }
    expect(quedas).toBe(0);
  });

  it('nao esta invertida: 0 e o extremo escuro, 1 o claro', () => {
    expect(luminancia(cor(0))).toBeLessThan(0.15);
    expect(luminancia(cor(1))).toBeGreaterThan(0.8);
    expect(luminancia(cor(1))).toBeGreaterThan(luminancia(cor(0)) + 0.6);
  });

  it('a rampa de fase e CICLICA: 0 e 2pi fecham na mesma cor', () => {
    /* Um mapa sequencial num observavel periodico inventa uma descontinuidade
       onde a fisica nao tem nenhuma. */
    const a = corFase(0), b = corFase(1);
    expect(a[0]).toBeCloseTo(b[0], 6);
    expect(a[1]).toBeCloseTo(b[1], 6);
    expect(a[2]).toBeCloseTo(b[2], 6);
  });

  it('a barra de escala usa as mesmas paradas do shader', () => {
    const g = gradienteCss();
    const [r, gg, bb] = cor(0);
    expect(g).toContain(`rgb(${Math.round(r)},${Math.round(gg)},${Math.round(bb)})`);
  });

  it('satura fora de [0,1] em vez de dar a volta', () => {
    expect(cor(-5)).toEqual(cor(0));
    expect(cor(5)).toEqual(cor(1));
  });

  it('a LUT do shader e a mesma funcao', () => {
    const t = lut(256);
    for (const i of [0, 1, 77, 128, 255]) {
      const [r, g, b] = cor(i / 255);
      expect(t[3 * i]).toBe(Math.round(r));
      expect(t[3 * i + 1]).toBe(Math.round(g));
      expect(t[3 * i + 2]).toBe(Math.round(b));
    }
  });
});

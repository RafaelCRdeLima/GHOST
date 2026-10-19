/* trajetorias.test.ts — o custo mostrado é o custo pago.
 *
 * Fixture com resposta conhecida: os números de referência para a varredura
 * abaixo, calculados à mão.
 * Se o modelo derivar, o usuário escolhe o modo com base num número errado — e
 * a escolha só se mostra errada meses depois, quando a RAM acabar ou o disco
 * encher, que é tarde.
 */
import { describe, expect, it } from 'vitest';
import {
  amostrasDeTempo, bytesAcumularRho, bytesArquivarPsi, formatarBytes,
} from './trajetorias.ts';

/* A escala da varredura: N = 40 x 13 = 520, 40 amostras, 200 trajetorias. */
const N = 520, AMOSTRAS = 40, NTRAJ = 200;

describe('custo dos modos de saida', () => {
  it('reproduz os 173 MB por thread de acumular_rho', () => {
    const b = bytesAcumularRho(N, AMOSTRAS);
    expect(b / 1e6).toBeGreaterThan(170);
    expect(b / 1e6).toBeLessThan(175);
    expect(formatarBytes(b)).toBe('173 MB');
  });

  it('reproduz os 66 MB por celula de arquivar_psi', () => {
    const b = bytesArquivarPsi(N, AMOSTRAS, NTRAJ);
    expect(b / 1e6).toBeGreaterThan(66);
    expect(b / 1e6).toBeLessThan(67);
    expect(formatarBytes(b)).toBe('67 MB');
    /* e as 100 celulas do plano dao 6,6 GB */
    expect(formatarBytes(b * 100)).toBe('6.7 GB');
  });

  /* ANTI-VACUIDADE: um modelo que devolvesse constante passaria nos dois testes
     acima. Estes exigem que ele RESPONDA aos parametros, e com a lei certa. */
  it('acumular_rho cresce com N ao QUADRADO', () => {
    const a = bytesAcumularRho(100, 10), b = bytesAcumularRho(200, 10);
    expect(b / a).toBeCloseTo(4, 10);
  });

  it('arquivar_psi cresce com N LINEARMENTE, e com as trajetorias tambem', () => {
    expect(bytesArquivarPsi(200, 10, 5) / bytesArquivarPsi(100, 10, 5)).toBeCloseTo(2, 10);
    expect(bytesArquivarPsi(100, 10, 10) / bytesArquivarPsi(100, 10, 5)).toBeCloseTo(2, 10);
  });

  /* E a diferenca de LEI e o que decide a escolha: em rede pequena arquivar
     custa mais, em rede grande acumular custa mais, e o ponto de virada depende
     do numero de trajetorias. Se os dois crescessem igual, nao haveria escolha
     a fazer e o radio seria decoracao. */
  it('a escolha vira: qual e mais caro depende de N', () => {
    expect(bytesArquivarPsi(50, 40, 200)).toBeGreaterThan(bytesAcumularRho(50, 40));
    expect(bytesArquivarPsi(4000, 40, 200)).toBeLessThan(bytesAcumularRho(4000, 40));
  });

  it('amostras de tempo espelham dae_traj_amostras', () => {
    expect(amostrasDeTempo(400, 10)).toBe(41);
    expect(amostrasDeTempo(20, 5)).toBe(5);
    expect(amostrasDeTempo(0, 5)).toBe(0);
    expect(amostrasDeTempo(400, 0)).toBe(0);
  });
});

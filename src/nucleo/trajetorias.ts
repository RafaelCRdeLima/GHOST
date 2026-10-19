/* trajetorias.ts — o custo dos dois modos de saída, calculado dos parâmetros.
 *
 * Os dois modos resolvem problemas diferentes e a escolha é irreversível na
 * prática: quem escolheu `accumulate_rho` e depois mudou a definição do
 * observável não tem os ψ para recalcular, e a varredura inteira volta para a
 * fila. Por isso não há padrão implícito, e por isso o custo aparece ao lado
 * de cada opção em vez de numa nota de rodapé: um número que muda quando o
 * usuário mexe em N ou no número de trajetórias é lido; uma nota, não.
 */

export type ModoSaida = 'accumulate_rho' | 'archive_psi';

/** Amostras de tempo do acumulador: t = 0 e depois uma a cada `stride` passos.
 *  Espelha `dae_traj_amostras` em C — e é O MESMO cálculo, não uma segunda
 *  versão dele, porque o custo mostrado tem de ser o custo pago. */
export function amostrasDeTempo(nt: number, stride: number): number {
  if (!(nt > 0) || !(stride > 0)) return 0;
  return 1 + Math.floor(nt / stride);
}

const BYTES_COMPLEXO = 16;   /* dois f64 */

/** RAM por THREAD: uma matriz N×N complexa por amostra de tempo. */
export function bytesAcumularRho(n: number, amostras: number): number {
  return amostras * n * n * BYTES_COMPLEXO;
}

/** Disco por CÉLULA do plano: N amplitudes por amostra, por trajetória. */
export function bytesArquivarPsi(n: number, amostras: number, nTraj: number): number {
  return amostras * n * nTraj * BYTES_COMPLEXO;
}

export function formatarBytes(b: number): string {
  if (!Number.isFinite(b) || b <= 0) return '0 B';
  if (b < 1e3) return `${b.toFixed(0)} B`;
  if (b < 1e6) return `${(b / 1e3).toFixed(0)} kB`;
  if (b < 1e9) return `${(b / 1e6).toFixed(0)} MB`;
  return `${(b / 1e9).toFixed(1)} GB`;
}

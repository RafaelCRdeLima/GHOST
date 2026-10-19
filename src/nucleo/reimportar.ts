/* reimportar.ts — o CSV que volta do cluster.
 *
 * ESTA É A ÚNICA ENTRADA DO SISTEMA QUE NÃO NASCE AQUI DENTRO, e por isso é a
 * que mais precisa desconfiar. Um CSV de uma versão anterior do núcleo, ou com
 * as colunas em outra ordem, plotaria — e plotaria errado, sem nada quebrar.
 *
 * Três defesas:
 *
 *   1. PROCEDÊNCIA OBRIGATÓRIA. Sem `#! spec` e `#! core_hash` o arquivo é
 *      recusado. Não é possível saber que simulação ele descreve, e um gráfico
 *      sem essa resposta é pior que nenhum gráfico.
 *   2. O SPEC PASSA PELO PARSER ESTRITO em C, o mesmo do resto do sistema, via
 *      WASM. Se ele não passa, o CSV não entra.
 *   3. COLUNAS POR NOME, nunca por posição. Ordem diferente tem de funcionar;
 *      coluna ausente tem de aparecer como ausente, não como outra coluna.
 *
 * `core_hash` diferente é AVISO, não bloqueio: reler resultados antigos é
 * legítimo. Mas o aviso vai para a tela, não para o console — e a interface
 * marca o painel inteiro em bronze, porque o número deixou de ser calculado
 * aqui.
 */
import type { Params } from '../i18n/index.ts';

export interface Aviso { chave: string; params?: Params; grave: boolean; }

export interface Reimportado {
  nome: string;
  meta: Map<string, string>;
  colunas: string[];
  /** Por NOME de coluna. */
  dados: Map<string, Float64Array>;
  estado?: { re: Float64Array; im: Float64Array };
  avisos: Aviso[];
  utilizavel: boolean;
}

/* QUEM calculou, que é outra pergunta e o `core_hash` não responde: o mesmo
   core_hash pode ter produzido o número pelo núcleo em C ou pelo pacote
   Wolfram, que usa outro método (decomposição espectral, não Chebyshev).
   Um resultado do Wolfram não deve se passar por um do núcleo.

   Campo ausente é caso real — CSV de versão anterior — e vira `desconhecida`,
   nunca `c`: supor a origem seria atribuir ao arquivo uma procedência que ele
   não afirma. Valor não previsto sobrevive como literal, porque inventar um
   rótulo para ele seria a mesma suposição por outro caminho. */
export type Origem =
  | { tipo: 'c' } | { tipo: 'wolfram' } | { tipo: 'desconhecida' }
  | { tipo: 'literal'; literal: string };

export function origemImplementacao(v: string | undefined): Origem {
  if (v === 'c') return { tipo: 'c' };
  if (v === 'wolfram') return { tipo: 'wolfram' };
  if (v !== undefined && v.trim().length > 0) return { tipo: 'literal', literal: v.trim() };
  return { tipo: 'desconhecida' };
}

/** Devolve null se o spec for válido, ou a mensagem do parser em C. */
export type ValidarSpec = (json: string) => string | null;

export function reimportar(
  nome: string, texto: string, hashAtual: string, validarSpec: ValidarSpec,
): Reimportado {
  const avisos: Aviso[] = [];
  const meta = new Map<string, string>();
  const secoes: Array<{ cab: string[]; linhas: string[][] }> = [];
  let atual: { cab: string[]; linhas: string[][] } | null = null;

  for (const linha of texto.split('\n')) {
    if (linha.startsWith('#!')) {
      const t = linha.slice(2).trim();
      const i = t.indexOf(' ');
      meta.set(i < 0 ? t : t.slice(0, i), i < 0 ? '' : t.slice(i + 1));
      continue;
    }
    if (linha.startsWith('#') || linha.trim() === '') { atual = null; continue; }
    const campos = linha.split(',').map((c) => c.trim());
    if (!atual) { atual = { cab: campos, linhas: [] }; secoes.push(atual); continue; }
    atual.linhas.push(campos);
  }

  const spec = meta.get('spec');
  const hashArquivo = meta.get('core_hash');
  const falhar = (chave: string, params?: Params): Reimportado => {
    avisos.push({ chave, params, grave: true });
    return { nome, meta, colunas: [], dados: new Map(), avisos, utilizavel: false };
  };

  if (!spec || !hashArquivo) return falhar('av_sem_procedencia');

  const erro = validarSpec(spec);
  if (erro) return falhar('av_spec_invalido', { erro });

  if (hashArquivo !== hashAtual) {
    /* Aviso, não bloqueio: reler resultado antigo é legítimo. Mas visível. */
    avisos.push({ chave: 'av_hash_diferente',
                  params: { arquivo: hashArquivo, atual: hashAtual }, grave: false });
  }

  const tabela = secoes[0];
  if (!tabela || tabela.linhas.length === 0) return falhar('av_sem_procedencia');

  /* Por NOME. A ordem das colunas no arquivo é irrelevante por construção. */
  const dados = new Map<string, Float64Array>();
  tabela.cab.forEach((nomeCol, i) => {
    const v = new Float64Array(tabela.linhas.length);
    for (let r = 0; r < tabela.linhas.length; ++r) {
      const bruto = tabela.linhas[r][i];
      v[r] = bruto === undefined || bruto === 'nan' ? Number.NaN : Number(bruto);
    }
    dados.set(nomeCol, v);
  });

  let estado: Reimportado['estado'];
  const sec2 = secoes[1];
  if (sec2 && sec2.cab.includes('re') && sec2.cab.includes('im')) {
    const ire = sec2.cab.indexOf('re'), iim = sec2.cab.indexOf('im');
    const re = new Float64Array(sec2.linhas.length);
    const im = new Float64Array(sec2.linhas.length);
    for (let r = 0; r < sec2.linhas.length; ++r) {
      re[r] = Number(sec2.linhas[r][ire]);
      im[r] = Number(sec2.linhas[r][iim]);
    }
    estado = { re, im };
  }

  avisos.push({ chave: 'av_reimportado',
                params: { arquivo: nome, linhas: tabela.linhas.length }, grave: false });
  return { nome, meta, colunas: tabela.cab, dados, estado, avisos, utilizavel: true };
}

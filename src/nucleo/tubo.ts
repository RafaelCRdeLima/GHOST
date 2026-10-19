/* tubo.ts — o microtúbulo enrolado, em três dimensões.
 *
 * O núcleo guarda a rede DESENROLADA: `xy[j] = (m, q)`, com m ao longo do eixo
 * e q em volta da circunferência. O cilindro é derivado aqui, e de propósito —
 * o grafo não tem geometria tridimensional, tem topologia; o tubo é uma
 * LEITURA dela, e leitura é trabalho da interface.
 *
 * ---------------------------------------------------------------------------
 * O EMBUTIMENTO É HELICOIDAL, E ISSO NÃO É ENFEITE
 *
 * A tentação é empilhar anéis planos: z = m, e cada anel de nq vértices num
 * plano. Fica errado, e erra exatamente onde a física está.
 *
 * A costura liga (m, nq−1) a (m + seam, 0). Com anéis planos essa ligação salta
 * `seam` unidades no eixo enquanto todas as outras ligações laterais não saltam
 * nada — a costura aparece como uma cicatriz geométrica, uma aresta esticada
 * atravessando o tubo. Mas ela não é isso. No microtúbulo real a costura é uma
 * descontinuidade de REDE (contato α–α onde o resto tem α–β), não de forma: os
 * monômeros continuam encaixados como todos os outros.
 *
 * O embutimento correto é a hélice de `seam` partidas:
 *
 *     z = m + seam · q / nq
 *
 * Com ele cada passo lateral sobe seam/nq, inclusive o da costura, e TODA
 * ligação lateral fica com o mesmo comprimento. A costura some da forma e
 * permanece na topologia, que é onde ela mora. `tubo.test.ts` mede isso, e mede
 * junto que o empilhamento ingênuo falharia — senão a afirmação passa sozinha.
 *
 * ---------------------------------------------------------------------------
 * O RAIO
 *
 * R = nq / 2π faz o arco entre vizinhos laterais valer 1, a mesma distância dos
 * vizinhos longitudinais. Sem isso o tubo sai espremido ou inflado, e a
 * proporção entre as duas direções de acoplamento — que é justamente o que
 * j_par e j_perp controlam — apareceria distorcida no desenho.
 *
 * ---------------------------------------------------------------------------
 * E MESMO ASSIM O EIXO PRECISA ENCURTAR, O QUE É UMA DISTORÇÃO ASSUMIDA
 *
 * Um microtúbulo é MUITO mais comprido que largo. No padrão do aplicativo são
 * 160 anéis de 13 protofilamentos: em escala métrica o tubo tem proporção 39:1
 * e vira um fio de cabelo atravessando o painel. Não há ângulo de câmera que
 * resolva — encurtar por perspectiva exigiria olhar quase pelo cano, e aí os
 * 160 anéis se sobrepõem e não se vê mais onde a excitação está.
 *
 * Então `compressao` encolhe o eixo. Isso QUEBRA de propósito a igualdade de
 * escala que o parágrafo acima defende: com o eixo comprimido, distância ao
 * longo do tubo e distância em volta dele deixam de ser comparáveis, e a razão
 * j_par/j_perp não se lê mais na figura.
 *
 * Por isso o fator aparece na tela em vez de ficar embutido. A vista existe
 * para mostrar como a rede se FECHA — que o retângulo desenrolado é um cilindro
 * e onde está a costura — e essa é uma pergunta topológica; para a pergunta
 * métrica existe o valor 1, que devolve a escala verdadeira.
 *
 * O que a compressão NÃO estraga: ela é um fator uniforme no eixo, então todas
 * as ligações laterais continuam com o mesmo comprimento entre si, costura
 * inclusive. `tubo.test.ts` verifica isso com compressão ligada.
 */

/** Coordenadas do cilindro, 3 floats por vértice: (eixo, y, z).
 *  `compressao` = 1 é a escala métrica verdadeira; abaixo disso o eixo encurta.
 */
export function coordenadasTubo(nPar: number, nPerp: number, seamShift: number,
                                compressao = 1): Float32Array {
  const np = Math.max(1, nPar | 0), nq = Math.max(1, nPerp | 0);
  const c = compressao > 0 ? compressao : 1;
  const xyz = new Float32Array(3 * np * nq);
  const R = nq / (2 * Math.PI);
  for (let m = 0; m < np; ++m) {
    for (let q = 0; q < nq; ++q) {
      const j = m * nq + q;
      const th = (2 * Math.PI * q) / nq;
      xyz[3 * j] = (m + (seamShift * q) / nq) * c;
      xyz[3 * j + 1] = R * Math.cos(th);
      xyz[3 * j + 2] = R * Math.sin(th);
    }
  }
  return xyz;
}

/** Passo lateral em unidades do eixo — quanto a hélice sobe por protofilamento. */
export function passoDaHelice(nPerp: number, seamShift: number): number {
  return seamShift / Math.max(1, nPerp | 0);
}

/** Compressão que põe o tubo na proporção `alvo` (comprimento : diâmetro).
 *  Nunca ESTICA: um tubo já curto fica em escala métrica, com fator 1. */
export function compressaoParaCaber(nPar: number, nPerp: number, alvo = 3.5): number {
  const diametro = Math.max(1, nPerp | 0) / Math.PI;
  const comprimento = Math.max(1, (nPar | 0) - 1);
  return Math.min(1, (alvo * diametro) / comprimento);
}

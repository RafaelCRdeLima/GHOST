import type { Catalog } from './types.ts';

/** Português — a língua em que o programa foi escrito, e a régua das outras. */
export const pt: Catalog = {
  modo_local: () => 'Laboratório · local',
  modo_reimportado: () => 'Reimportado',
  impl_c: () => 'núcleo C',
  impl_wolfram: () => 'Wolfram',
  impl_desconhecida: () => 'origem não declarada',

  sec_gerador: () => 'Gerador',
  sec_hamiltoniano: () => 'Hamiltoniano',
  sec_tempo: () => 'Tempo',
  sec_sitios: () => 'Sítios',
  sec_exportar: () => 'Exportar',
  sec_varredura: () => 'Varredura',
  sec_reimportar: () => 'Reimportar',

  ger_microtubule: () => 'microtúbulo',
  ger_sbm: () => 'blocos estocásticos',
  ger_path: () => 'linha',
  ger_cycle: () => 'ciclo',
  ger_grid2d: () => 'grade 2D',
  ger_hypercube: () => 'hipercubo',
  ger_complete: () => 'completo',

  c_npar: () => 'N∥',
  c_nperp: () => 'N⊥',
  c_costura: () => 'costura',
  c_fechar: () => 'fechar as pontas',
  c_acoplamento: () => 'j∥ / j⊥',
  c_modulos: () => 'módulos',
  c_n: () => 'N',
  c_pinpout: () => 'p_in / p_out',
  c_religacao: () => 'religação p',
  c_religar_fixo: () => 'religar mantendo |E| fixo',
  c_semente: () => 'semente',
  c_gamma: () => 'γ',
  c_normalizacao: () => 'normalização de ‖H‖',
  c_tfinal: () => 't final',
  c_pontos: () => 'pontos no tempo',
  c_realizacoes: () => 'realizações',
  c_passos: () => 'passos em p',

  c_passo_quadro: () => 'passo do mapa',
  q_quadros: (p) => `${p['k']} quadros · ${p['mb']} MB`,
  q_incompleto: (p) => `animação parou em ${p['tem']} de ${p['esperado']} quadros`,
  ham_adjacency: () => 'Adjacência',
  ham_laplacian: () => 'Laplaciano',
  norm_spectral: () => 'raio espectral',
  norm_mean_degree: () => 'grau médio',
  norm_none: () => 'nenhuma',

  a_gerar: () => 'gerar rede',
  a_propagar: () => 'propagar',
  a_propagando: () => 'propagando…',
  a_parar: () => 'parar',
  a_escala_fixa: () => 'escala fixa',
  a_varrer: () => 'varrer',
  a_csv: () => 'CSV dos resultados',
  a_abrir_csv: () => 'abrir CSV do cluster',
  a_limpar_alvo: () => 'limpar alvo',

  s_inicial: (p) => `inicial ${p['j']}`,
  s_alvo: (p) => `alvo ${p['j']}`,
  s_dica: () => 'Clique no mapa para marcar o sítio selecionado acima.',

  av_acrescentar: () =>
    'Acrescentar arestas aumenta |E| e melhora o transporte trivialmente: '
    + 'o resultado passa a medir o número de arestas, não a topologia.',
  av_sem_normalizar: () =>
    'Sem normalizar, "mais coerência" pode ser apenas "hopping maior" ao '
    + 'comparar topologias diferentes.',
  av_nao_lattice: () =>
    'Este grafo não tem rede desenrolada: o mapa mostra os sítios em ordem de '
    + 'índice, dobrados em linhas. Não é o layout do grafo.',
  av_exportar: () =>
    'O C++ é autocontido e regenera a rede a partir do spec.json. Wolfram e '
    + 'Python recebem a lista de arestas e servem como oráculo do propagador.',
  av_sem_procedencia: () =>
    'Este CSV não tem o cabeçalho de procedência (#! spec e #! core_hash). '
    + 'Não dá para saber que simulação ele descreve, então não será plotado.',
  av_hash_diferente: (p) =>
    `Este resultado veio do núcleo ${p['arquivo']}; aqui está rodando o `
    + `${p['atual']}. Os números continuam válidos, mas não são reproduzíveis `
    + 'por este binário — para reproduzir, use aquela versão do núcleo.',
  av_spec_invalido: (p) => `O spec.json embutido no CSV não passa no parser: ${p['erro']}`,
  av_reimportado: (p) => `${p['arquivo']} · ${p['linhas']} passos`,
  av_scrub_grande: () => 'rede grande: só ao vivo, sem voltar no tempo',

  r_norma: () => 'Σ|ψ|²',
  r_vertices: (p) => `${p['n']} vértices`,
  r_arestas: (p) => `${p['m']} arestas`,
  r_grau: () => '⟨d⟩',
  r_lambda2_alto: () => 'não convergiu · limite superior',
  r_componentes: (p) => `${p['k']} componentes`,
  r_duplicadas: (p) => `${p['k']} duplicadas descartadas`,
  r_religacoes: (p) => `${p['k']} religações sem destino`,

  f_iniciando: () => 'iniciando',
  f_gerando: () => 'gerando rede',
  f_pronta: () => 'rede pronta',
  f_propagando: () => 'propagando',
  f_pronto: () => 'pronto',
  f_cancelado: () => 'cancelado',
  f_erro: () => 'erro',
  f_varrendo: (p) => `varrendo ${p['i']}/${p['n']}`,

  sec_vista: () => 'Vista',
  v_desenrolada: () => 'rede desenrolada',
  v_espectral: () => 'layout espectral',
  v_propria: () => 'geometria do gerador',
  v_espectral_dica: () => 'Autovetores 2 e 3 da laplaciana como coordenadas. Para rede modular, os módulos se separam no plano.',
  v_tubo: () => 'Tubo 3D',
  v_tubo_dica: () => 'Arraste para girar. A costura não aparece como emenda porque o embutimento é helicoidal: toda ligação lateral tem o mesmo comprimento, inclusive a dela. A costura é descontinuidade de rede — contato α–α onde o resto tem α–β — não de forma.',
  v_tubo_fechado: () => 'Com o tubo fechado longitudinalmente, as ligações que unem as duas pontas atravessam o cilindro no desenho: o grafo é periódico e o cilindro não pode ser.',
  v_tubo_metrico: () => 'escala métrica',
  v_tubo_escala: (p) => `Eixo comprimido ×${p['fator']}: ao longo do tubo as distâncias NÃO estão na mesma escala das laterais. Um microtúbulo é longo demais para caber em escala — em ×1 a escala é verdadeira e o tubo vira um fio.`,
  sec_trajetorias: () => 'trajetórias',
  tr_desligado: () => 'desligado',
  tr_n: () => 'trajetórias',
  tr_gamma: () => 'γ defasagem',
  tr_stride: () => 'amostragem de ρ',
  tr_amostras: () => 'amostras',
  tr_rho: () => 'acumular ρ',
  tr_psi: () => 'arquivar ψ',
  tr_por_thread: () => 'por thread',
  tr_por_celula: () => 'por célula',
  tr_plano: () => 'no plano de 100',
  tr_sem_modo_curto: () => 'sem modo',
  tr_rho_quando: () => 'ρ é somado durante a execução e cada ψ é descartado. Escolha para varredura ampla, com os observáveis já definidos.',
  tr_psi_quando: () => 'os ψ são gravados e ρ é montado na análise. Escolha quando a definição do observável ainda puder mudar.',
  tr_sem_modo: () => 'Escolha um dos dois modos. Não há padrão: eles resolvem problemas diferentes, e quem escolheu acumular ρ e depois mudou a definição do observável não tem os ψ para recalcular.',
  v_sem_espectral: () => 'rede grande demais para o layout espectral',
  v_sem_arestas: () => 'arestas ocultas acima de 2000 vértices',
  se_ipr: () => 'IPR',
  se_coerencia: () => 'coerência ℓ₁',
  se_palvo: () => 'p no alvo',
  se_passo: () => 'passo',
  v_eixo: () => 'p (fração religada)',
  v_media: () => 'p̄ no alvo (média ± desvio)',
  v_dica: () => 'Varredura pequena roda aqui. Para varreduras longas, exporte o C++ com realizações.',
};

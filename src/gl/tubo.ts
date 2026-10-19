/* tubo.ts — o microtúbulo enrolado, desenhado.
 *
 * A geometria vem de src/nucleo/tubo.ts; aqui só se projeta e se desenha.
 *
 * ---------------------------------------------------------------------------
 * A PROJEÇÃO É FEITA EM JS, NÃO NO SHADER
 *
 * Parece desperdício mandar posições já projetadas para a GPU em vez de mandar
 * a matriz. Não é: `sitioNoEvento` precisa saber onde cada vértice foi parar na
 * tela, e se o shader projetasse, a mesma projeção existiria duas vezes — uma
 * em GLSL e outra em TypeScript. Duas implementações de uma conta é o arranjo
 * em que elas divergem, e a divergência aqui apareceria como "clico num vértice
 * e seleciona o vizinho", que ninguém lê como bug de matriz.
 *
 * A conta roda uma vez por mudança de orientação, não por quadro. Com N na casa
 * dos milhares isso é nada.
 *
 * ---------------------------------------------------------------------------
 * O QUE CODIFICA PROFUNDIDADE, E O QUE NÃO PODE CODIFICAR
 *
 * Cor e raio do vértice já significam |ψ_j|², e essa redundância é deliberada —
 * é ela que deixa a figura legível sem cor. Então nem cor
 * nem raio podem carregar profundidade: um vértice do outro lado do tubo
 * ficaria menor, e menor já quer dizer menos população.
 *
 * Sobram três cues, e nenhum deles mente sobre o valor:
 *
 *   perspectiva na POSIÇÃO — o lado de cá se abre, o de lá se fecha;
 *   oclusão — teste de profundidade, o vértice da frente cobre o de trás;
 *   as ARESTAS desbotam com a distância. Elas já eram "contexto, não conteúdo"
 *     na nuvem espectral, então podem carregar isso sem disputar significado
 *     com nada.
 *
 * E o cue mais forte é o giro, que é do usuário: arrastar resolve qualquer
 * ambiguidade que uma projeção estática deixe.
 */
import { lut } from '../nucleo/paleta';

const VS = `#version 300 es
in vec2 aPos;
in float aProf;
in float aVal;
uniform vec2 uEscala;
uniform vec2 uCentro;
uniform float uMax;
uniform float uRaio;
out float vT;
void main() {
  vec2 p = (aPos - uCentro) * uEscala;
  /* aProf em [0,1], 0 = perto. Vira z do NDC para o teste de profundidade
     ocluir o fundo; NÃO entra em gl_PointSize. */
  gl_Position = vec4(p.x, -p.y, aProf * 2.0 - 1.0, 1.0);
  vT = uMax > 0.0 ? clamp(aVal / uMax, 0.0, 1.0) : 0.0;
  gl_PointSize = uRaio * (0.55 + 1.35 * sqrt(vT));
}`;

const FS = `#version 300 es
precision highp float;
uniform sampler2D uLut;
in float vT;
out vec4 cor;
void main() {
  vec2 d = gl_PointCoord - vec2(0.5);
  float r = length(d);
  if (r > 0.5) discard;
  float borda = smoothstep(0.5, 0.42, r);
  cor = vec4(texture(uLut, vec2(vT, 0.5)).rgb, borda);
}`;

const VS_ARESTA = `#version 300 es
in vec2 aPos;
in float aProf;
uniform vec2 uEscala;
uniform vec2 uCentro;
out float vProf;
void main() {
  vec2 p = (aPos - uCentro) * uEscala;
  gl_Position = vec4(p.x, -p.y, aProf * 2.0 - 1.0, 1.0);
  vProf = aProf;
}`;

const FS_ARESTA = `#version 300 es
precision highp float;
in float vProf;
out vec4 cor;
/* Cinza-cera da identidade. O alfa cai com a profundidade: é este gradiente que
   faz o cilindro parecer cilindro numa figura parada. A faixa é estreita de
   propósito — some demais e o outro lado do tubo desaparece, e o tubo vira
   meia-cana.

   Os valores são baixos porque são MUITAS arestas: o microtúbulo padrão tem
   4144, e com alfa alto elas se somam até virar uma superfície opaca que cobre
   os vértices. Aí a figura fica bonita e perde a função — quem tem de ser lido
   é |ψ_j|², e a malha é só o andaime que diz onde ele está. */
void main() { cor = vec4(0.486, 0.533, 0.573, mix(0.155, 0.03, vProf)); }`;

function programa(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
  const compilar = (tipo: number, fonte: string) => {
    const s = gl.createShader(tipo)!;
    gl.shaderSource(s, fonte);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(`shader: ${gl.getShaderInfoLog(s)}`);
    return s;
  };
  const p = gl.createProgram()!;
  gl.attachShader(p, compilar(gl.VERTEX_SHADER, vs));
  gl.attachShader(p, compilar(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(`programa: ${gl.getProgramInfoLog(p)}`);
  return p;
}

/** Teto próprio, e mais alto que o da nuvem espectral de propósito: lá as
 *  arestas são contexto e acima de 2000 vértices viram névoa; aqui elas são o
 *  que revela a forma, e a estrutura regular do tubo não vira névoa. */
export const MAX_ARESTAS_DESENHADAS = 8000;

/** Inclinação inicial: o tubo levemente de perfil, mostrando que é redondo.
 *  De frente ele viraria uma faixa e o 3D não se anunciaria. */
export const ELEVACAO_PADRAO = 0.42;

export class Tubo {
  private gl: WebGL2RenderingContext;
  private progPto: WebGLProgram;
  private progAre: WebGLProgram;
  private bufPos: WebGLBuffer;
  private bufProf: WebGLBuffer;
  private bufVal: WebGLBuffer;
  private bufAre: WebGLBuffer;
  private texLut: WebGLTexture;

  private n = 0;
  private nAre = 0;
  private modelo: Float32Array<ArrayBufferLike> = new Float32Array(0);
  private ligacoes: Array<[number, number]> = [];
  private centroModelo: [number, number, number] = [0, 0, 0];
  private foco = 1;

  private azim = 0;
  private elev = ELEVACAO_PADRAO;

  private pos = new Float32Array(0);
  private prof = new Float32Array(0);
  private centro: [number, number] = [0, 0];
  private meia: [number, number] = [1, 1];

  constructor(private canvas: HTMLCanvasElement) {
    const gl = canvas.getContext('webgl2',
      { antialias: true, alpha: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL2 indisponivel neste navegador');
    this.gl = gl;
    this.progPto = programa(gl, VS, FS);
    this.progAre = programa(gl, VS_ARESTA, FS_ARESTA);
    this.bufPos = gl.createBuffer()!;
    this.bufProf = gl.createBuffer()!;
    this.bufVal = gl.createBuffer()!;
    this.bufAre = gl.createBuffer()!;
    this.texLut = gl.createTexture()!;
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB8, 256, 1, 0, gl.RGB, gl.UNSIGNED_BYTE, lut(256));
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  }

  /** `xyz` são 3 floats por vértice, de coordenadasTubo. */
  rede(xyz: Float32Array, n: number, arestas?: Array<[number, number, number]> | null) {
    this.n = n;
    this.modelo = xyz;
    this.pos = new Float32Array(2 * n);
    this.prof = new Float32Array(n);

    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; ++i) { cx += xyz[3 * i]; cy += xyz[3 * i + 1]; cz += xyz[3 * i + 2]; }
    this.centroModelo = [cx / n, cy / n, cz / n];

    let raio = 1e-9;
    for (let i = 0; i < n; ++i) {
      const d = Math.hypot(xyz[3 * i] - this.centroModelo[0],
                           xyz[3 * i + 1] - this.centroModelo[1],
                           xyz[3 * i + 2] - this.centroModelo[2]);
      if (d > raio) raio = d;
    }
    /* Distância focal proporcional ao tamanho do objeto: a perspectiva fica
       igual num tubo de 20 anéis e num de 200. Fixa em unidades absolutas, o
       tubo longo receberia uma deformação absurda e o curto, nenhuma. */
    this.foco = 2.6 * raio;

    this.ligacoes = [];
    if (arestas && n <= MAX_ARESTAS_DESENHADAS) {
      for (const [i, j] of arestas) if (i !== j) this.ligacoes.push([i, j]);
    }
    this.projetar();
  }

  /** Orientação absoluta, em radianos. */
  orientar(azim: number, elev: number) {
    this.azim = azim;
    /* Trava a elevação antes do polo: passar de 90° inverte o sentido do
       arrasto e o usuário perde a referência de onde está o eixo do tubo. */
    this.elev = Math.max(-1.45, Math.min(1.45, elev));
    this.projetar();
  }

  orientacao(): [number, number] { return [this.azim, this.elev]; }

  /* Uma só implementação da projeção, e é esta. Ver o cabeçalho. */
  private projetar() {
    const n = this.n;
    if (n === 0) return;
    const gl = this.gl;
    const [ox, oy, oz] = this.centroModelo;
    const ca = Math.cos(this.azim), sa = Math.sin(this.azim);
    const ce = Math.cos(this.elev), se = Math.sin(this.elev);

    const px = new Float32Array(n), py = new Float32Array(n), pz = new Float32Array(n);
    let zlo = Infinity, zhi = -Infinity;
    for (let i = 0; i < n; ++i) {
      const x = this.modelo[3 * i] - ox;
      const y = this.modelo[3 * i + 1] - oy;
      const z = this.modelo[3 * i + 2] - oz;
      /* gira em torno do EIXO DO TUBO (x do modelo): escolhe que protofilamento
         fica de frente */
      const y1 = y * ca - z * sa;
      const z1 = y * sa + z * ca;
      /* inclina em torno da vertical da tela: é isto que aponta o eixo do tubo
         para dentro ou para fora, e é o que faz a figura ser tridimensional em
         vez de uma faixa */
      const x2 = x * ce + z1 * se;
      const z2 = -x * se + z1 * ce;
      px[i] = x2; py[i] = y1; pz[i] = z2;
      if (z2 < zlo) zlo = z2;
      if (z2 > zhi) zhi = z2;
    }

    const alcance = Math.max(zhi - zlo, 1e-9);
    let xlo = Infinity, xhi = -Infinity, ylo = Infinity, yhi = -Infinity;
    for (let i = 0; i < n; ++i) {
      const w = this.foco / (this.foco + (pz[i] - zlo));
      const X = px[i] * w, Y = py[i] * w;
      this.pos[2 * i] = X; this.pos[2 * i + 1] = Y;
      this.prof[i] = (pz[i] - zlo) / alcance;
      if (X < xlo) xlo = X; if (X > xhi) xhi = X;
      if (Y < ylo) ylo = Y; if (Y > yhi) yhi = Y;
    }
    this.centro = [(xlo + xhi) / 2, (ylo + yhi) / 2];
    this.meia = [Math.max((xhi - xlo) / 2, 1e-9), Math.max((yhi - ylo) / 2, 1e-9)];

    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
    gl.bufferData(gl.ARRAY_BUFFER, this.pos, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufProf);
    gl.bufferData(gl.ARRAY_BUFFER, this.prof, gl.DYNAMIC_DRAW);

    this.nAre = 0;
    if (this.ligacoes.length > 0) {
      const dados = new Float32Array(this.ligacoes.length * 6);
      let k = 0;
      for (const [i, j] of this.ligacoes) {
        dados[k++] = this.pos[2 * i]; dados[k++] = this.pos[2 * i + 1]; dados[k++] = this.prof[i];
        dados[k++] = this.pos[2 * j]; dados[k++] = this.pos[2 * j + 1]; dados[k++] = this.prof[j];
      }
      this.nAre = this.ligacoes.length * 2;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufAre);
      gl.bufferData(gl.ARRAY_BUFFER, dados, gl.DYNAMIC_DRAW);
    }
  }

  private ajuste(): { escala: [number, number]; k: number; L: number; A: number } {
    const cw = this.canvas.clientWidth || 1, ch = this.canvas.clientHeight || 1;
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    if (this.canvas.width !== Math.round(cw * dpr) || this.canvas.height !== Math.round(ch * dpr)) {
      this.canvas.width = Math.round(cw * dpr);
      this.canvas.height = Math.round(ch * dpr);
    }
    const L = this.canvas.width, A = this.canvas.height;
    const margem = 22 * dpr;
    /* Escala UNIFORME, como na nuvem: a projeção já distorce distância pela
       perspectiva, e esticar por eixo em cima disso tiraria qualquer leitura
       de forma que ainda restasse. */
    const k = Math.min((L / 2 - margem) / this.meia[0], (A / 2 - margem) / this.meia[1]);
    return { escala: [(2 * k) / L, (2 * k) / A], k, L, A };
  }

  desenhar(valores: Float32Array, max: number) {
    const gl = this.gl;
    const { escala, L, A } = this.ajuste();
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    gl.viewport(0, 0, L, A);
    gl.clearColor(0.043, 0.071, 0.098, 1);
    gl.clearDepth(1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    if (this.n === 0) return;

    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);

    if (this.nAre > 0) {
      /* Arestas NÃO escrevem profundidade: se escrevessem, cada linha recortaria
         os vértices que passam atrás dela e a figura ficaria picotada. Elas
         leem o buffer, então somem atrás dos vértices da frente, que é o que se
         quer. */
      gl.depthMask(false);
      gl.useProgram(this.progAre);
      gl.uniform2f(gl.getUniformLocation(this.progAre, 'uEscala'), escala[0], escala[1]);
      gl.uniform2f(gl.getUniformLocation(this.progAre, 'uCentro'), this.centro[0], this.centro[1]);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufAre);
      const aP = gl.getAttribLocation(this.progAre, 'aPos');
      gl.enableVertexAttribArray(aP);
      gl.vertexAttribPointer(aP, 2, gl.FLOAT, false, 12, 0);
      const aD = gl.getAttribLocation(this.progAre, 'aProf');
      gl.enableVertexAttribArray(aD);
      gl.vertexAttribPointer(aD, 1, gl.FLOAT, false, 12, 8);
      gl.drawArrays(gl.LINES, 0, this.nAre);
      gl.depthMask(true);
    }

    gl.useProgram(this.progPto);
    gl.uniform2f(gl.getUniformLocation(this.progPto, 'uEscala'), escala[0], escala[1]);
    gl.uniform2f(gl.getUniformLocation(this.progPto, 'uCentro'), this.centro[0], this.centro[1]);
    gl.uniform1f(gl.getUniformLocation(this.progPto, 'uMax'), max);
    gl.uniform1f(gl.getUniformLocation(this.progPto, 'uRaio'),
                 Math.max(3, Math.min(16, 340 / Math.sqrt(this.n))) * dpr);
    gl.uniform1i(gl.getUniformLocation(this.progPto, 'uLut'), 0);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);

    const aPos = gl.getAttribLocation(this.progPto, 'aPos');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
    gl.enableVertexAttribArray(aPos);
    gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

    const aProf = gl.getAttribLocation(this.progPto, 'aProf');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufProf);
    gl.enableVertexAttribArray(aProf);
    gl.vertexAttribPointer(aProf, 1, gl.FLOAT, false, 0, 0);

    const aVal = gl.getAttribLocation(this.progPto, 'aVal');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufVal);
    gl.bufferData(gl.ARRAY_BUFFER, valores, gl.DYNAMIC_DRAW);
    gl.enableVertexAttribArray(aVal);
    gl.vertexAttribPointer(aVal, 1, gl.FLOAT, false, 0, 0);

    gl.drawArrays(gl.POINTS, 0, this.n);
    gl.disable(gl.DEPTH_TEST);
  }

  /** Vértice mais próximo do clique, ou -1. Desempata pelo MAIS PERTO do
   *  observador: no tubo dois vértices caem sobre o mesmo pixel o tempo todo, e
   *  escolher o de trás seria escolher o que está escondido. */
  sitioNoEvento(ev: { clientX: number; clientY: number }): number {
    if (this.n === 0) return -1;
    const r = this.canvas.getBoundingClientRect();
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    const { k, L, A } = this.ajuste();
    const px = ((ev.clientX - r.left) / r.width) * L;
    const py = ((ev.clientY - r.top) / r.height) * A;
    const lim = (26 * dpr) ** 2;
    let melhor = -1, melhorProf = Infinity;
    for (let i = 0; i < this.n; ++i) {
      const x = L / 2 + (this.pos[2 * i] - this.centro[0]) * k;
      const y = A / 2 + (this.pos[2 * i + 1] - this.centro[1]) * k;
      if ((x - px) ** 2 + (y - py) ** 2 > lim) continue;
      if (this.prof[i] < melhorProf) { melhorProf = this.prof[i]; melhor = i; }
    }
    return melhor;
  }
}

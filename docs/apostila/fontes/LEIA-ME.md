# Fontes da identidade visual

As três famílias usadas na apostila — **Newsreader** (serifada),
**Manrope** (sem serifa) e **DM Mono** (monoespaçada) — ficam aqui em vez de
serem exigidas do sistema, para que o build da apostila não dependa do que cada
máquina tem instalado.

Todas são distribuídas sob a **SIL Open Font License 1.1**, que permite
redistribuição desde que a licença acompanhe os arquivos. Ela está em
`OFL-newsreader.txt`, `OFL-manrope.txt` e `OFL-dmmono.txt`.

## Por que instâncias estáticas

Newsreader e Manrope são publicadas como **fontes variáveis**, e o XeLaTeX usa
apenas a instância padrão de uma fonte variável — o que daria um documento
inteiro num único peso. Os arquivos aqui foram gerados com
`fontTools.varLib.instancer` a partir das originais:

| arquivo | eixos fixados |
|---|---|
| `Newsreader-Regular.ttf` | `opsz=16, wght=400` |
| `Newsreader-Medium.ttf` | `opsz=16, wght=500` |
| `Newsreader-SemiBold.ttf` | `opsz=16, wght=600` |
| `Newsreader-Display.ttf` | `opsz=36, wght=500` — para a capa |
| `Newsreader-Italic.ttf` | `opsz=16, wght=400` |
| `Newsreader-SemiBoldItalic.ttf` | `opsz=16, wght=600` |
| `Manrope-{Light,Regular,Medium,SemiBold,Bold}.ttf` | `wght=300…700` |

DM Mono já é publicada em instâncias estáticas e foi usada como veio.

## Letras gregas

**Newsreader e DM Mono não têm gregas.** Por isso o preâmbulo carrega STIX Two
Math via `unicode-math`, e toda letra grega no documento aparece em modo
matemático. Nunca escreva `ν` no texto corrido: sai um retângulo vazio.

## Regenerar

```python
from fontTools.varLib import instancer
from fontTools.ttLib import TTFont
f = TTFont("Newsreader[opsz,wght].ttf")
instancer.instantiateVariableFont(f, {"opsz": 16, "wght": 600}, inplace=True,
                                  updateFontNames=True)
f.save("Newsreader-SemiBold.ttf")
```

As variáveis originais vêm de `github.com/google/fonts/ofl/{newsreader,manrope,dmmono}`.

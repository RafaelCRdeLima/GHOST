#!/bin/sh
# wasm/build.sh — o mesmo nucleo, agora para o navegador.
#
#   ./wasm/build.sh modulo   -> daedalus.mjs + .wasm (a ponte, para a interface)
#   ./wasm/build.sh cli      -> dump.js (teste 6) e a suite de testes sob Node
#   ./wasm/build.sh          -> os dois
#
# Sem pthreads e sem SharedArrayBuffer de proposito: exigiriam cabecalhos
# COOP/COEP e complicariam a hospedagem estatica simples.
#
# -ffp-contract=off e -fno-fast-math sao explicitos. -ffast-math quebraria o NaN
# de p_alvo e a reprodutibilidade bit a bit, e quebraria em silencio; o portao
# de verdade e `make -C native wasm-test`, que roda a suite inteira sob WASM.
set -eu
RAIZ=$(cd "$(dirname "$0")/.." && pwd)
SAIDA="$RAIZ/wasm/build"
ALVO=${1:-tudo}

command -v emcc >/dev/null 2>&1 || {
  echo "emcc nao encontrado. Ative o emsdk:  source ~/emsdk/emsdk_env.sh" >&2
  exit 1
}

COMUM="-O3 -msimd128 -std=c99 -Wall -Wextra -Werror -pedantic
       -ffp-contract=off -fno-fast-math -I $RAIZ/core"
NUCLEO=$(ls "$RAIZ"/core/*.c)

mkdir -p "$SAIDA"
# A raiz do projeto e "type": "module"; os artefatos de linha de comando do
# emscripten sao CommonJS.
echo '{ "type": "commonjs" }' > "$SAIDA/package.json"

node "$RAIZ/tools/hash_core.mjs" >/dev/null

if [ "$ALVO" = "modulo" ] || [ "$ALVO" = "tudo" ]; then
  # shellcheck disable=SC2086
  emcc $COMUM $NUCLEO "$RAIZ/wasm/bridge.c" -o "$SAIDA/daedalus_core.mjs" \
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=criarDaedalus \
    -sENVIRONMENT=web,worker \
    -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=33554432 \
    -sEXPORTED_RUNTIME_METHODS=wasmMemory,UTF8ToString,stringToUTF8,lengthBytesUTF8 \
    -sEXPORTED_FUNCTIONS=@"$RAIZ/wasm/exports.txt" \
    -sFILESYSTEM=0 -sASSERTIONS=0 -sSINGLE_FILE=1
  cp "$RAIZ/wasm/daedalus.mjs" "$SAIDA/"
  # SINGLE_FILE embute o .wasm em base64 dentro do .mjs. Custa ~35% de bytes e
  # elimina a questao do caminho do asset na hospedagem estatica: um import, um
  # arquivo, nenhuma configuracao de servidor.
  #
  # E e por isso que ENVIRONMENT nao inclui `node`: o unico codigo especifico de
  # Node no artefato do emscripten e o que le o .wasm do disco, e sem arquivo
  # para ler ele some. O modulo continua rodando sob Node (os testes de fixture
  # dependem disso) e o pacote do navegador para de arrastar um import de
  # node:module que o Vite tem de externalizar.
  echo "  modulo:  $(du -h "$SAIDA/daedalus_core.mjs" | cut -f1)  (wasm embutido)"
fi

if [ "$ALVO" = "cli" ] || [ "$ALVO" = "tudo" ]; then
  # shellcheck disable=SC2086
  emcc $COMUM -I "$RAIZ/native/tests" $NUCLEO \
    "$RAIZ/native/tests/scenarios.c" "$RAIZ/native/dump.c" \
    -o "$SAIDA/dump.js" -sENVIRONMENT=node -sALLOW_MEMORY_GROWTH=1 -sEXIT_RUNTIME=1
  for t in "$RAIZ"/native/tests/t*.c; do
    nome=$(basename "$t" .c)
    # shellcheck disable=SC2086
    emcc $COMUM -I "$RAIZ/native/tests" $NUCLEO \
      "$RAIZ/native/tests/harness.c" "$RAIZ/native/tests/jacobi.c" "$t" \
      -o "$SAIDA/$nome.js" -sENVIRONMENT=node -sALLOW_MEMORY_GROWTH=1 -sEXIT_RUNTIME=1
  done
  # shellcheck disable=SC2086
  emcc $COMUM $NUCLEO "$RAIZ/native/bench.c" \
    -o "$SAIDA/bench.js" -sENVIRONMENT=node -sALLOW_MEMORY_GROWTH=1 -sEXIT_RUNTIME=1
  echo "  cli:     dump + bench + $(ls "$RAIZ"/native/tests/t*.c | wc -l) testes"
fi

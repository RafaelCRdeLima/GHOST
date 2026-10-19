# native/Makefile — binário nativo e testes de aceitação.
#
# -ffast-math é PROIBIDO: quebraria o NaN de p_alvo e a reprodutibilidade.
ROOT      := ..
CC        ?= gcc
CXX       ?= g++
OPT       ?= -O2 -g
WARN      := -Wall -Wextra -Werror -pedantic
CFLAGS    := -std=c99 $(OPT) $(WARN) -I$(ROOT)/core -Itests
LDLIBS    := -lm
# Resolvido para caminho absoluto de proposito: com o emsdk ativo, ~/emsdk/node
# e um DIRETORIO na frente do PATH, e o make executa receitas simples sem shell
# — sem pular diretorios, como o shell faz. Sem isto, "node" vira
# "Permissao negada" so quando o emsdk esta carregado.
NODE      := $(shell command -v node 2>/dev/null || echo node)
# Prefixo de execução dos testes. O TSan não convive com o ASLR de 32 bits dos
# kernels 6.x ("unexpected memory mapping"), então o alvo tsan roda sob
# setarch -R. Sem isso o sanitizador aborta antes de olhar o programa.
RUN       ?=
BUILD     := build

VERSION_H := $(ROOT)/core/dae_version.generated.h
CORE_SRC  := $(sort $(wildcard $(ROOT)/core/*.c))
CORE_HDR  := $(filter-out $(VERSION_H),$(sort $(wildcard $(ROOT)/core/*.h)))
CORE_OBJ  := $(patsubst $(ROOT)/core/%.c,$(BUILD)/core_%.o,$(CORE_SRC))
TEST_SRC  := $(sort $(wildcard tests/t*.c))
TEST_BIN  := $(patsubst tests/%.c,$(BUILD)/bin/%,$(TEST_SRC))
SUP_OBJ   := $(BUILD)/harness.o $(BUILD)/jacobi.o

.PHONY: all test test-cli bench cxx-check asan ubsan tsan clean version \
        dump cli wasm wasm-test wasm-bench test6 test7
all: $(TEST_BIN) $(BUILD)/bin/bench $(BUILD)/bin/ghost $(BUILD)/bin/sonda

version:
	@$(NODE) $(ROOT)/tools/hash_core.mjs

$(VERSION_H): $(ROOT)/core/amalgam.list $(CORE_SRC) $(CORE_HDR) $(ROOT)/tools/hash_core.mjs
	@$(NODE) $(ROOT)/tools/hash_core.mjs

$(BUILD) $(BUILD)/bin:
	@mkdir -p $@

$(BUILD)/core_%.o: $(ROOT)/core/%.c $(VERSION_H) | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: tests/%.c $(VERSION_H) | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/bin/bench: bench.c $(CORE_OBJ) | $(BUILD)/bin
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/bin/dump: dump.c tests/scenarios.c $(CORE_OBJ) | $(BUILD)/bin
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/bin/sonda: sonda.c tests/jacobi.c $(CORE_OBJ) | $(BUILD)/bin
	$(CC) $(CFLAGS) -I../core -Itests $^ -o $@ -lm

$(BUILD)/bin/ghost: cli.c $(CORE_OBJ) | $(BUILD)/bin
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

dump: $(BUILD)/bin/dump
cli: $(BUILD)/bin/ghost

$(BUILD)/bin/%: $(BUILD)/%.o $(CORE_OBJ) $(SUP_OBJ) | $(BUILD)/bin
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# O CLI e casca, e os testes de aceitacao exercitam o NUCLEO. Um bug de casca
# (`traj` e `run` discordarem sobre o que e um grafo valido, por exemplo)
# nao aparece em teste de nucleo, por melhor que ele seja.
test-cli: $(BUILD)/bin/ghost
	@sh $(ROOT)/tools/teste_cli.sh

test: all test-cli
	@fail=0; for b in $(TEST_BIN); do $(RUN) ./$$b || fail=1; echo; done; \
	if [ $$fail -ne 0 ]; then echo "TESTES FALHARAM"; exit 1; fi; \
	echo "todos os testes de aceitacao passaram"

bench: $(BUILD)/bin/bench
	@./$(BUILD)/bin/bench

# O núcleo tem de compilar limpo pelos dois caminhos: é o que garante que a
# amalgamação vira um .cpp válido sem uma segunda implementação.
cxx-check: $(VERSION_H) | $(BUILD)
	@$(NODE) $(ROOT)/tools/amalgamate.mjs $(BUILD)/dae_amalgam.c
	$(CC)  -std=c99   $(WARN) -I$(ROOT)/core -fsyntax-only $(BUILD)/dae_amalgam.c
	$(CXX) -std=c++17 $(WARN) -I$(ROOT)/core -fsyntax-only -x c++ $(BUILD)/dae_amalgam.c
	@echo "nucleo amalgamado: limpo em C99 e em C++17"

asan:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" test

tsan:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O1 -g -fsanitize=thread" RUN="setarch -R" test

# ---------------------------------------------------------------- WebAssembly
# Exigem o emsdk ativo:  source ~/emsdk/emsdk_env.sh
WBUILD := $(ROOT)/wasm/build

wasm:
	@sh $(ROOT)/wasm/build.sh

# ACEITACAO 6: mesma cena, mesmo nucleo, dois compiladores. A tolerancia e
# 1e-14, mas o piso de identidade bit a bit em 99% e o que vigia o MECANISMO:
# os dois alvos batem em 100%, porque -msimd128 nao reassocia soma de
# ponto flutuante. Cair para "dentro da tolerancia mas reassociado" e sinal de
# flag nova, e tem de aparecer como falha, nao como numero um pouco pior.
test6: $(BUILD)/bin/dump
	@sh $(ROOT)/wasm/build.sh cli >/dev/null
	@./$(BUILD)/bin/dump > $(BUILD)/dump_nativo.txt
	@$(NODE) $(WBUILD)/dump.js > $(BUILD)/dump_wasm.txt
	@echo "aceitacao 6: WASM contra nativo"
	@$(NODE) $(ROOT)/tools/compare_dump.mjs $(BUILD)/dump_nativo.txt $(BUILD)/dump_wasm.txt 1e-14 99

# A suite inteira sob WASM. E o portao que pega -ffast-math herdado: t91 exige
# que p_alvo seja NaN sem alvo, e -ffast-math faria isso virar zero em silencio.
wasm-test:
	@sh $(ROOT)/wasm/build.sh >/dev/null
	@fail=0; for t in $(WBUILD)/t*.js; do $(NODE) $$t || fail=1; done; \
	$(NODE) $(ROOT)/wasm/teste_memoria.mjs || fail=1; \
	if [ $$fail -ne 0 ]; then echo "SUITE WASM FALHOU"; exit 1; fi; \
	echo "suite completa passa sob WebAssembly"

wasm-bench:
	@sh $(ROOT)/wasm/build.sh cli >/dev/null
	@$(NODE) $(WBUILD)/bench.js


clean:
	@rm -rf $(BUILD)

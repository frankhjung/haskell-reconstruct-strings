#!/usr/bin/env make

.DEFAULT_GOAL	:= default

GHC_VERSION	:= 9.6.7
TARGET		:= reconstruct-strings
CABAL		:= $(TARGET).cabal
SRCS		:= $(wildcard app/*.hs src/*.hs src/**/*.hs test/*.hs)

.PHONY: default
default: format check build test ## Run the default pipeline

.PHONY: all
all: format check build test doc exec ## Run the full pipeline

.PHONY: help
help: ## Show this help message
	@echo ""
	@echo "Default goal: ${.DEFAULT_GOAL}"
	@awk 'BEGIN {FS = ":.*##"; \
		printf "\nUsage:\n  make \033[36m<target>\033[0m\n\nTargets:\n"} \
		/^[a-zA-Z_-]+:.*?##/ { \
		printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2 }' \
		$(MAKEFILE_LIST)

.PHONY: format
format: $(SRCS) ## Format cabal file and Haskell sources
	@echo format ...
	@cabal-fmt --inplace $(CABAL)
	@stylish-haskell --inplace $(SRCS)

.PHONY: check
check: tags lint ## Run static checks

.PHONY: tags
tags: $(SRCS) ## Generate ctags for Haskell sources
	@echo tags ...
	@hasktags --ctags --extendedctag $(SRCS)

.PHONY: lint
lint: $(SRCS) ## Run hlint and cabal check
	@echo lint ...
	@hlint --cross --color --show $(SRCS)
	@cabal check

.PHONY: build
build: $(SRCS) ## Build project with cabal
	@echo build ...
	@cabal build

.PHONY: test
test: ## Run test suite
	@echo test ...
	@cabal test --test-show-details=direct

.PHONY: doc
doc: ## Build Haddock documentation
	@echo doc ...
	@cabal haddock \
		--haddock-executables \
		--haddock-quickjump \
		--haddock-hyperlink-sources \
		lib:$(TARGET) \
		exe:$(TARGET)

.PHONY: exec
exec: build ## Run sample DNA string reconstruction
	@echo "Reconstructing reads [ATGGC, GGCGT, CGTGCA] (min overlap = 2):"
	@cabal exec $(TARGET) -- -m 2 ATGGC GGCGT CGTGCA
	@echo ""
	@echo "Reconstructing reads [ABC, BCD, CDE] (min overlap = 2):"
	@cabal exec $(TARGET) -- -m 2 ABC BCD CDE

.PHONY: setup
setup: ## Init cabal config and update dependencies
ifeq (,$(wildcard ${CABAL_CONFIG}))
	-cabal user-config init
else
	@echo Using user-config from ${CABAL_CONFIG} ...
endif
	-cabal update --only-dependencies

.PHONY: ghci
ghci: ## Open GHCi via cabal repl
	@cabal repl

.PHONY: clean
clean: ## Clean build artifacts and tags
	-cabal clean
	-$(RM) tags

.PHONY: cleanall
cleanall: clean ## Purge build artifacts, cache, and tags
	-$(RM) -r dist-newstyle

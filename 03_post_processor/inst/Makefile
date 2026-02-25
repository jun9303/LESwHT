PYTHON ?= python3

.PHONY: all flib post_inst post clean

all: post_inst

flib:
	$(MAKE) -C .. flib

post_inst: flib
	$(PYTHON) post_inst.py

post: post_inst

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f ./*.pyc

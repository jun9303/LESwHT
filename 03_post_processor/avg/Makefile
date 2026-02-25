PYTHON ?= python3

.PHONY: all flib post_avg post clean

all: post_avg

flib:
	$(MAKE) -C .. flib

post_avg: flib
	$(PYTHON) post_avg.py

post: post_avg

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f ./*.pyc

NAUTY_DIR = ./nauty
CC = gcc
CFLAGS = -O3 -I$(NAUTY_DIR)
LDFLAGS = $(NAUTY_DIR)/nauty.a
REPORT = resultats.html

all: check_iso

check_iso: check_iso.c
	$(CC) $(CFLAGS) -o check_iso check_iso.c $(LDFLAGS)
 
run: check_iso
	python main.py "$(URL)"
	explorer.exe resultats.html || true

clean:
	rm -f check_iso.exe $(REPORT)
	rm -rf data/molecules/*.mol data/graphs/*.graph
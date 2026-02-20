NAUTY_VER     = 2_9_3
NAUTY_DIR     = nauty$(NAUTY_VER)
NAUTY_ARCHIVE = nauty$(NAUTY_VER).tar.gz
NAUTY_URL     = https://pallini.di.uniroma1.it/$(NAUTY_ARCHIVE)
NAUTY_LIB     = $(NAUTY_DIR)/nauty.a

CC            = gcc
PYTHON        = python3
CFLAGS        = -O3 -I$(NAUTY_DIR)
LDFLAGS       = $(NAUTY_LIB)

EXEC          = Check_iso
SRC           = Check_iso.c
REPORT        = resultats.html
ARCHIVE_NAME  = projet_code.zip
SRC_FILES     = Main.py Check_iso.c Makefile

.PHONY: run clean zip

run: $(EXEC)
	@if [ -z "$(URL)" ]; then echo "Usage: make run URL='http://...'"; exit 1; fi
	$(PYTHON) Main.py "$(URL)"

$(EXEC): $(SRC) $(NAUTY_LIB)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

$(NAUTY_LIB):
	@echo "--- Téléchargement ---"
	wget $(NAUTY_URL)
	@echo "--- Extraction de Nauty ---"
	rm -rf $(NAUTY_DIR)
	tar xvzf $(NAUTY_ARCHIVE)
	touch $(NAUTY_DIR)
	rm $(NAUTY_ARCHIVE)
	@echo "--- Configuration de Nauty ---"
	cd $(NAUTY_DIR) && ./configure
	@echo "--- Compilation de Nauty ---"
	cd $(NAUTY_DIR) && $(MAKE)

zip:
	rm -f $(EXEC) $(EXEC).exe $(REPORT)
	rm -rf data
	rm -rf $(NAUTY_DIR) $(NAUTY_ARCHIVE) $(ARCHIVE_NAME)
	rm -f $(ARCHIVE_NAME)
	zip -r $(ARCHIVE_NAME) $(SRC_FILES)

clean:
	rm -f $(EXEC) $(EXEC).exe $(REPORT)
	rm -rf data
	rm -rf $(NAUTY_DIR) $(NAUTY_ARCHIVE) $(ARCHIVE_NAME)

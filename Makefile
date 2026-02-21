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
SRC_FILES     = Main.py Check_iso.c Makefile Distance.py template.html

OPT_A         = $(if $(A),-a $(A),)

.PHONY: run clean zip

run: $(EXEC)
	@if [ -z "$(HC)" ]; then \
		echo "Erreur : Tu dois fournir le paramètre HC (nombre de clusters)."; \
		echo "Exemple : make run SDF_FILE='chemin/vers/fichier.sdf' HC=5"; \
		exit 1; \
	fi
	@if [ -n "$(URL)" ]; then \
		echo "--- Lancement avec l'URL ---"; \
		$(PYTHON) Main.py "$(URL)" $(OPT_A) -hc $(HC); \
	elif [ -n "$(SDF_FILE)" ]; then \
		echo "--- Lancement avec le fichier local ---"; \
		$(PYTHON) Main.py "$(SDF_FILE)" $(OPT_A) -hc $(HC); \
	elif [ -n "$(DATABASE)" ]; then \
		echo "--- Extraction des IDs depuis $(DATABASE) ---"; \
		IDS=$$(grep -E '^[0-9]+,' $(DATABASE) | cut -d',' -f1 | paste -sd, -); \
		echo "--- Téléchargement du fichier SDF via PubChem ---"; \
		curl -X POST -d "cid=$$IDS" "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/SDF" -o db_molecules.sdf; \
		echo "--- Lancement du script sur les molécules téléchargées ---"; \
		$(PYTHON) Main.py db_molecules.sdf $(OPT_A) -hc $(HC); \
		rm db_molecules.sdf; \
	else \
		echo "Erreur : Tu dois fournir soit URL, soit SDF_FILE, soit DATABASE."; \
		echo "Usage 1 : make run URL='http://...' HC=5 [A=0.8] [OUTPUT='dossier']"; \
		echo "Usage 2 : make run SDF_FILE='chemin/vers/fichier.sdf' HC=10 [A=0.5] [OUTPUT='dossier']"; \
		echo "Usage 3 : make run DATABASE='dataset.txt' HC=3 [A=0.9] [OUTPUT='dossier']"; \
		exit 1; \
	fi
	@if [ -n "$(OUTPUT)" ]; then \
		echo "--- Déplacement des résultats vers $(OUTPUT) ---"; \
		mkdir -p "$(OUTPUT)"; \
		rm -rf "$(OUTPUT)/data" "$(OUTPUT)/$(REPORT)"; \
		mv data "$(OUTPUT)/"; \
		mv $(REPORT) "$(OUTPUT)/"; \
	fi

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
	rm -f $(EXEC) $(EXEC).exe $(REPORT) db_molecules.sdf
	rm -rf data
	rm -rf $(NAUTY_DIR) $(NAUTY_ARCHIVE) $(ARCHIVE_NAME)
	rm -f $(ARCHIVE_NAME)
	zip -r $(ARCHIVE_NAME) $(SRC_FILES)

clean:
	rm -f $(EXEC) $(EXEC).exe $(REPORT) db_molecules.sdf
	rm -rf data
	rm -rf $(NAUTY_DIR) $(NAUTY_ARCHIVE) $(ARCHIVE_NAME)

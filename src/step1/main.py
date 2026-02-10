from requestAPI import request_chebi_api

# main qui demande en console de 1 a N motifs a chercher
if __name__ == "__main__":

    chebi_motif_list = []

    print("Entrez les motifs ChEBI à chercher (entrez 'ok' quand vous avez fini) :")
    while True:
        motif = input("Motif : ")
        if motif.lower() == 'ok':
            break
        chebi_motif_list.append(motif)

    request_chebi_api(chebi_motif_list)
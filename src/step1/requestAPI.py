# Requete vers une api chEBI
import requests
import json

def request_chebi_api(chebi_motif_list : list):
    """"
    Fonction qui récupère les molfiles des motifs ChEBI passés en argument
    (cmd) et les sauvegarde dans des fichiers .mol

    In : chebi_motif_list : list de motifs ChEBI à chercher
    Out : None
    """

    # Recuperration de l'url de la requete
    with open('src/step1/requestConfig.json') as f:
        config = json.load(f)
        url_request_motif_part1 = config["url_request_motif_part1"]
        url_request_motif_part2 = config["url_request_motif_part2"]
        url_request_id = config["url_request_id"]

    for motif in chebi_motif_list:
        request1 = url_request_motif_part1 + motif + url_request_motif_part2
        response1 = requests.get(request1)
        data1 = response1.json()
        # print(f"data for motif {motif} : {data1}")

        # On recupere l'id du premier resultat
        motif_id = data1["results"][0]["_id"]
        print(f"Motif '{motif}' a pour ChEBI ID : {motif_id}")

        # On fait une deuxieme requete pour recuperer le molfile
        request2 = url_request_id + str(motif_id)
        response2 = requests.get(request2)
        molfile_content = response2.text

        # On save le molfile dans un fichier .mol
        with open(f"src/molfiles/{motif_id}.mol", "w") as molfile:
            molfile.write(molfile_content)
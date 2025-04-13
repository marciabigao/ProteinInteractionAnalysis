import requests
import os

#download sequence in FASTA format
pdb_id = intput("Digite o ID da sequência de entrada: ")
fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"

fasta_response = requests.get(fasta_url)
if not fasta_response.ok:
    print("Error obtaining FASTA sequence.")
    exit()

fasta_lines = fasta_response.text.strip().slip('\n')
sequence = ''.join(line.strip() for line in fasta_lines if not line.startswith('>'))

#search for similar sequencies
query = {
    "query": {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
            "evalue_cutoff": 1, #lower means more precise search
            "identity_cutoff": 0.9, #minimum identity between sequencies (90%)
            "target": "pdb_protein_sequece", #compare chains directly (more precision)
            "sequence": sequence_only #sequence used as entry
        }
    },
    "retun_type": "entry", #define the results as PDB IDs
    "request_options": {
        "pager": {"start": 0, "rows": 10}, #starts from first result and limits to 10 (maximum = 100)
        "scoring_strategy": "sequence" #order the results from sequece similarity
    }
}

search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
search_response = requests.post(search_url, json=query)

if search_response.ok:
    similar_sequencies = search_response.json().get("result_set", []) #take the "result_set" list from response or an empty list
else:
    print("Search error: ", search_response.status_code)
    exit()

#download .cif files from similar_sequencies
for line in similar_sequencies:
    download_id = line["identifier"]
    cif_url = f"https://files.rcsb.org/download/{download_id}.cif"
    cif_response = requests.get(cif_url)

    folder_path = "../data"
    os.makedirs(folder_path, exist_ok=True) #make sure folder exists

    if cif_response.ok:
        file_path = os.path.join(folder_path, f"{download_id}.cif")
        with open(file_path, 'w') as f:
            f.write(cif_response.text)
    else:
        print("Download error in {download_id}.cif")


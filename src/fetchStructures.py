import requests
import os

#download sequence in FASTA format
pdb_id = input("Digite o ID da sequência de entrada: ").upper()
chain_id = input("Digite a cadeia desejada: ").upper()
fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"

fasta_response = requests.get(fasta_url)
if not fasta_response.ok:
    print("Error obtaining FASTA sequence.")
    exit()

fasta_lines = fasta_response.text.strip().split('\n')

#separate wanted chain
sequence = ''
collecting = False

for line in fasta_lines:
    if line.startswith('>'):
        if f"{pdb_id}" and f"{chain_id}" in line:
            collecting = True
            continue
        else:
            continue
    elif collecting:
        sequence += line.strip()


#search for similar sequencies
query = {
    "query": {
    "type": "terminal",
    "service": "sequence",
    "parameters": {
      "evalue_cutoff": 0.01,
      "identity_cutoff": 0.9,
      "target": "pdb_protein_sequence",
      "sequence_type": "protein",
      "value": sequence
    }
  },
  "request_options": {
    "paginate": {
        "start": 0,
        "rows": 100
    }
  },
  "return_type": "polymer_entity",
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
    download_id = line["identifier"].strip()
    clean_id = download_id.split('_')[0] #remove chain specificatins
    cif_url = f"https://files.rcsb.org/download/{clean_id}.cif"
    cif_response = requests.get(cif_url)

    folder_path = "../data"
    os.makedirs(folder_path, exist_ok=True) #make sure folder exists

    if cif_response.ok:
        file_path = os.path.join(folder_path, f"{clean_id}.cif")
        with open(file_path, 'w') as f:
            f.write(cif_response.text)
        print("Successful download of " + clean_id + ".cif")
    else:
        print("Download error in " + clean_id + ".cif")

import requests
import os

def fetchStructures(folder):
    #access fasta files
    fasta_lines = []
    for file in os.listdir(folder):
        path = os.path.join(folder, file)
        with open(path, 'r') as f:
            fasta_lines = f.readlines()

        #separate wanted chain
        sequence = ''
        collecting = False

        for line in fasta_lines:
            if line.startswith('>'):
                collecting = True
                continue
            elif collecting:
                sequence += line.strip()
                collecting = False


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

            data_folder = file.split(".")[0]
            folder_path = f"../data/{data_folder}"
            os.makedirs(folder_path, exist_ok=True) #make sure folder exists

            if cif_response.ok:
                file_path = os.path.join(folder_path, f"{clean_id}.cif")
                with open(file_path, 'w') as f:
                    f.write(cif_response.text)
                print("Successful download of " + clean_id + ".cif")
            else:
                print("Download error in " + clean_id + ".cif")

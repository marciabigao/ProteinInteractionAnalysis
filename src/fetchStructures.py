import requests

#download sequence in FASTA format
pdb_id = intput("Digite o ID da sequência de entrada: ")
fasta_url = f"https://www.rcsb.org/fasta/entry{pdb_id}/display"

fasta_response = requests.get(fasta_url)
if not fasta_response.ok:
    print("Error obtaining FASTA sequence.")
    exit()

fasta_lines = fasta_response.text.strip().slip('\n')
sequence = ''.join(line.strip() for line in fasta_lines if not line.startswith('>'))

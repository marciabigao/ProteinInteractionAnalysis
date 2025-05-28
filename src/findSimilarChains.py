import requests
import os
from subprocess import run

def findChain(tipo_entrada, variante):
    resultado = []

    pasta = f"../data/{variante}"
    fasta_original = f"../entradas/{tipo_entrada}/{variante}.fasta"
    with open(fasta_original, "r") as f:
        linhas = f.readlines()
        seq_referencia = linhas[1].strip()

    for arquivo in os.listdir(pasta):
        id = arquivo.split('.')[0]
        pdb_id = id.lower()
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        resposta = requests.get(url)

        if not resposta.ok:
            continue #vai para o proximo id se download não foi possível

        #armazenar sequencias das cadeias vindas do fasta
        cadeias = {}
        bloco = ""

        for linha in resposta.text.strip().splitlines():
            if linha.startswith('>'): #caso a linha seja um cabeçalho (ou seja, se eu já tiver mudado para outra cadeia)
                if bloco: #se já existirem informações no bloco (essas informações serão o cabeçalho e a sequência da cadeia anterior)
                    linhas = bloco.split("\n") #processo a cadeia anterior
                    header = linhas[0]
                    seq = linhas[1:]
                    cadeias[header] = "".join(seq) #
                bloco = linha #adiciono o novo cabeçalho
            else:
                bloco += "\n" + linha #adiciona linhas que não são cabeçalhos ao bloco
        
        if bloco: #necessario para garantir que estou processando as informações do último bloco
            linhas = bloco.split("\n")
            header = linhas[0]
            seq = linhas[1:]
            cadeias[header] = "".join(seq)

        #analisa cadeias recuperadas para ver qual é a mais parecida com a entrada
        melhor_cadeia = ""
        melhor_identidade = -1

        #para cada uma das cadeias que eu recuperei, eu preciso criar um fasta temporário para usar o método clustlw
        for header, cadeia_seq in cadeias.items():
            # Extrair o identificador da cadeia individual (ex: A, B, etc.)
            if "|" in header:
                partes = header.split("|")
                cadeia_info = partes[1].strip()
                # Extrai apenas o primeiro identificador (ex: A) mesmo que venha "Chains A, B"
                cadeia_ids = cadeia_info.replace("Chains", "").replace("Chain", "").strip()
                cadeia_ids = cadeia_ids.replace(",", "").split()
                for cadeia_id in cadeia_ids:
                    fasta_temp = f"temp_{pdb_id}_{cadeia_id}.fasta"
                    with open(fasta_temp, "w") as f:
                        f.write(">ref\n" + seq_referencia + "\n")
                        f.write(f">{cadeia_id}\n" + cadeia_seq + "\n")

                    #executo o clustalw a partir de uma linha de comando do programa
                    cmd = ["clustalw", f"-infile={fasta_temp}"] #esse é o comando que seria executado no terminal
                    result = run(cmd, capture_output=True, text=True)

                    if result.returncode != 0:
                        print("Erro ao rodar ClustalW:")
                        print(result.stderr)
                        return

                    #resultado é um arquivo .aln com o mesmo nome do fasta de entrada no mesmo diretorio, então abro e leio seu conteudo
                    aln_file = fasta_temp.replace(".fasta", ".aln")
                    with open(aln_file) as f:
                        aln = f.read()

                    #limpar arquvios temporarios -> melhorar essa organização depois
                    for ext in [".aln", ".dnd", ".fasta"]:
                        path = fasta_temp.replace(".fasta", ext)
                        if os.path.exists(path):
                            os.remove(path)

                    #devo contato também as substituições conservadas (:) e semi-conservadas (.)?
                    identidade = aln.count("*") #um asterisco significa que as posições comparadas são idênticas
                    if identidade > melhor_identidade:
                        melhor_identidade = identidade
                        melhor_cadeia = cadeia_id 
                

        #adicionar o resultado encontrado na lista
        if melhor_cadeia:
            resultado.append((pdb_id, melhor_cadeia))
    
    #adiciona todos os resultados a um arquivo, para serem usados no listInteractions
    with open(f"../results/cadeias_alvo/cadeias_alvo_{variante}.txt", "w") as f:
        for pdb_id, cadeia in resultado:
            pdb_id = pdb_id.upper()
            f.write(f"{pdb_id} {cadeia}\n")
                

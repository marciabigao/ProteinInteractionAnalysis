import os
import math
from collections import defaultdict

#variaveis utilizadas na comparacao das distancias
pasta = "../data"
cadeia_alvo = input("Digite a cadeia alvo da busca: ").strip().upper()
distancia_minima = 6.00
interacoes_arquivos = defaultdict(set) #armazena resultados das buscas

#funcao de calculo da distancia
def distancia_euclidiana(a, b):
    return sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

#abrindo arquivos e registrando seus dados
    for arquivo in os.listdir(pasta):
        caminho = os.path.join(pasta, arquivo)
        with open(caminho, "r") as f:
            linhas = f.readlines()

        #encontrando o bloco _atom_site, que contem as informações dos átomos
        atomos[]
        indices{} #armazena o nome das colunas e os indices de referencia
        in_atom_site = False

        for linha in linhas:
            if linha.strip().startswith("loop_"):
                in_atom_site = False #quando detecto que estou entrando em uma nova tabela, eu paro de ler os átomos que estava lendo (condicao de parada)

            if "_atom_site." in linha:
                in_atom_site = True
                chave = linha.strip().split(".")[1] #pego o identificador da coluna que vem depois de "_atom_site."
                indices[chave] = len(indices) #coloco como identificador o tamanho do dicionario, como um iterador
                continue #faco isso ate ler todas as linhas do cabecalho

            if in_atom_site and linha.strip() and not linha.strastswith("_"):
                instancia = linha.strip().split()

                try:
                    atomos.append({
                        "residuo": instancia[indices["label_comp_id"]], #o que o atomo e em termos biologicos/quimicos (aminoacido, ligante, ion...)
                        "cadeia": instancia[indices["label_asym_id"]],
                        "x": float(partes[indices["Cartn_x"]]),
                        "y": float(partes[indices["Cartn_y"]]),
                        "z": float(partes[indices["Cartn_z"]]),
                    })
                except (IndexError, ValueError, KeyError):
                    continue #continua pegando outras linhas caso de erro em alguma

    atomos_alvo = [a for a in atomos if a["cadeia"] == cadeia_alvo]
    outros_atomos = [a for a in atomos if a["cadeia"] != cadeia_alvo and a["residuo"] not in ("HOH", "WAT")] #atomos que nao sao da cadeia pesquisada nem agua

    #calcula as distancias e verifica interacao
    for a1 in atomos_alvo:
        for a2 in outros_atomos:
            d = distancia((a1["x"], a1["y"], a1["z"]), (a2["x"], a2["y"], a2["z"]))
            if d <= distancia_minima:
                identificador = f"{a2["residuo"]} (Cadeia {a2["cadeia"]})"
                interacoes_arquivos[arquivo].add(identificador)


#adiciona informacoes de interacao a um arquivo
with open("../results/interactionList.txt", "w") as file:
    for arquivo, interacoes in interacoes_arquivos.itens():
        file.write(f"Arquivo {arquivo}\n")

        if not interacoes:
            file.write("Nenhuma interacao")
        else:
            for i in interacoes:
                file.write(f"{i}\n")
    file.write("\n")



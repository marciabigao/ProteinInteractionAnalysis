import os
import math
from collections import defaultdict
import shlex

#dicionarios que mapeiam a letra da cadeia para seu nome
def construir_mapas_cadeias_e_descricoes(linhas):
    mapa_cadeia_para_entidade = {}
    mapa_entidade_para_descricao = {}

    #mapeamento da cadeia para o numero da entidade
    entidade_atual = None
    strands_atuais = None
    for linha in linhas:
        linha = linha.strip()
        if linha.startswith("_entity_poly.entity_id"):
            partes = linha.split()
            if len(partes) > 1:
                entidade_atual = partes[1]
        elif linha.startswith("_entity_poly.pdbx_strand_id"):
            partes = linha.split()
            if len(partes) > 1:
                strands_atuais = partes[1].replace(',', ' ').split()

        if entidade_atual and strands_atuais:
            for cadeia in strands_atuais:
                mapa_cadeia_para_entidade[cadeia.strip()] = entidade_atual
            entidade_atual = None
            strands_atuais = None

    #mapeamento do numero da entidade para sua descrição (nome)
    in_entity = False
    colunas_entity = []
    registros_entity = []

    for linha in linhas:
        linha = linha.strip()
        if linha.startswith("loop_"):
            in_entity = False

        if "_entity.id" in linha:
            in_entity = True
            colunas_entity = [linha]
            registros_entity = []
            continue

        if in_entity:
            if linha.startswith("_entity."):
                colunas_entity.append(linha)
            elif linha:
                registros_entity.append(linha)

    for registro in registros_entity:
        parts = shlex.split(registro)  # RESPEITA aspas simples e duplas!
        if len(parts) < len(colunas_entity):
            continue
        dados = dict(zip(colunas_entity, parts))
        entity_id = dados.get("_entity.id")
        descricao = dados.get("_entity.pdbx_description", "")
        descricao = descricao.strip("'\"")  # Remove aspas externas, se houver

        if entity_id:
            mapa_entidade_para_descricao[entity_id] = descricao

    return mapa_cadeia_para_entidade, mapa_entidade_para_descricao

def listInteractions(variante):
    #variaveis utilizadas na comparacao das distancias
    pasta = f"../data/{variante}"
    distancia_minima = 6.00
    interacoes_arquivos = defaultdict(set) #armazena resultados das buscas

    #encontrar cadeias alvo para cada .cif baixado
    cadeias_alvo = {}
    caminho_cadeias_alvo = f"../results/cadeias_alvo/cadeias_alvo_{variante}.txt"
    with open(caminho_cadeias_alvo, "r") as f:
        linhas_arquivo = f.readlines()
    for linha in linhas_arquivo:
        cif_id = linha.strip().split(" ")[0]
        cadeia_alvo = linha.strip().split(" ")[1]

        cadeias_alvo[cif_id] = cadeia_alvo

    #funcao de calculo da distancia
    def distancia_euclidiana(a, b):
        return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

    #abrindo arquivos e registrando seus dados
    for arquivo in os.listdir(pasta):
        cif_id = arquivo.split(".")[0]
        cadeia_alvo = cadeias_alvo[cif_id]
        caminho = os.path.join(pasta, arquivo)
        with open(caminho, "r") as f:
            linhas = f.readlines()
        mapa_cadeia_para_entidade, mapa_entidade_para_descricao = construir_mapas_cadeias_e_descricoes(linhas) #para identificar nome das cadeias
        #print(mapa_cadeia_para_entidade)
        #print(mapa_entidade_para_descricao)

        #encontrando o bloco _atom_site, que contem as informações dos átomos
        atomos = []
        indices = {} #armazena o nome das colunas e os indices de referencia
        in_atom_site = False

        for linha in linhas:
            if linha.strip().startswith("loop_"):
                in_atom_site = False #quando detecto que estou entrando em uma nova tabela, eu paro de ler os átomos que estava lendo (condicao de parada)

            if "_atom_site." in linha:
                in_atom_site = True
                chave = linha.strip().split(".")[1] #pego o identificador da coluna que vem depois de "_atom_site."
                indices[chave] = len(indices) #coloco como identificador o tamanho do dicionario, como um iterador
                continue #faco isso ate ler todas as linhas do cabecalho

            if in_atom_site and linha.strip() and not linha.startswith("_"):
                instancia = linha.strip().split()

                try:
                    atomos.append({
                        "residuo": instancia[indices["label_comp_id"]], #o que o atomo e em termos biologicos/quimicos (aminoacido, ligante, ion...)
                        "cadeia": instancia[indices["label_asym_id"]],
                        "x": float(instancia[indices["Cartn_x"]]),
                        "y": float(instancia[indices["Cartn_y"]]),
                        "z": float(instancia[indices["Cartn_z"]]),
                    })
                except (IndexError, ValueError, KeyError):
                    continue #continua pegando outras linhas caso de erro em alguma

        residuos_espurios = ("HOH", "WAT", "ACE", "ACT", "BME", "CSD", "CSW", "EDO", "FMT", "GOL", "MSE", "NAG", "NO3", "PO4", "SGM", "SO4", "TPO", "LU8", "U0C", "7IO", "09V", "XK7") #artefatos cristalográficos, água
        atomos_alvo = [a for a in atomos if a["cadeia"] == cadeia_alvo]
        outros_atomos = [a for a in atomos if a["cadeia"] != cadeia_alvo and a["residuo"] not in residuos_espurios] #atomos que nao sao da cadeia pesquisada nem espurios

        print(f"Processando {arquivo} com {len(atomos_alvo)} átomos na cadeia alvo e {len(outros_atomos)} átomos externos")

        #calcula as distancias e verifica interacao
        for a1 in atomos_alvo:
            for a2 in outros_atomos:
                d = distancia_euclidiana((a1["x"], a1["y"], a1["z"]), (a2["x"], a2["y"], a2["z"]))
                if d <= distancia_minima:
                    entity_id = mapa_cadeia_para_entidade.get(a2["cadeia"])
                    descricao = mapa_entidade_para_descricao.get(entity_id)

                    if descricao:
                        identificador = f"{a2['residuo']} (Cadeia {a2['cadeia']}, {descricao})"
                    else:
                        identificador = f"{a2['residuo']} (Cadeia {a2['cadeia']})"

                    interacoes_arquivos[arquivo].add(identificador)


    #adiciona informacoes de interacao a um arquivo
    with open(f"../results/interacoes/interaction_list_{variante}.txt", "w") as file:
        for arquivo, interacoes in interacoes_arquivos.items():
            file.write(f"Arquivo {arquivo}\n")

            if not interacoes:
                file.write("Nenhuma interacao")
            else:
                for i in interacoes:
                    file.write(f"{i}\n")
                file.write("\n")
        file.write("\n")



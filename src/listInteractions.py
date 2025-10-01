import os
import math
from collections import defaultdict
import shlex

def parse_cif_data(linhas):
    """
    analisa as linhas de um arquivo CIF em uma única passagem para extrair:
    1. um mapa de ID da cadeia para ID da entidade (extraído do _atom_site).
    2. um mapa de ID da entidade para sua descrição (extraído do _entity).
    3. uma lista de todos os átomos e suas coordenadas.
    """
    mapa_cadeia_para_entidade = {}
    mapa_entidade_para_descricao = {}
    atomos = []

    colunas_map = {}
    current_loop_type = None

    for linha in linhas:
        linha_strip = linha.strip()

        if not linha_strip:
            continue

        if linha_strip == "loop_":
            current_loop_type = None
            colunas_map.clear()
            continue

        if linha_strip.startswith("_"):
            #bloco _entity é necessário para as descrições
            if "_entity." in linha_strip and "_entity_poly." not in linha_strip:
                current_loop_type = "entity"
                colunas_map[linha_strip] = len(colunas_map)
            #bloco _atom_site é usado para os átomos e para o mapeamento cadeia->entidade
            elif "_atom_site." in linha_strip:
                current_loop_type = "atom_site"
                colunas_map[linha_strip] = len(colunas_map)
            #ignora outros blocos, incluindo o _entity_poly que causava problemas
            elif current_loop_type is not None:
                 colunas_map[linha_strip] = len(colunas_map)
            continue

        if current_loop_type and not linha_strip.startswith("#"):
            try:
                parts = shlex.split(linha_strip)

                if current_loop_type == "entity":
                    if len(parts) == len(colunas_map):
                        entity_id = parts[colunas_map["_entity.id"]]
                        descricao = parts[colunas_map["_entity.pdbx_description"]]
                        mapa_entidade_para_descricao[entity_id] = descricao.strip("'\"")

                elif current_loop_type == "atom_site":
                    # Mapeamento robusto de cadeia para entidade
                    cadeia = parts[colunas_map["_atom_site.auth_asym_id"]]
                    entity_id = parts[colunas_map["_atom_site.label_entity_id"]]
                    if cadeia not in mapa_cadeia_para_entidade:
                        mapa_cadeia_para_entidade[cadeia] = entity_id

                    # Extração dos dados do átomo
                    atomos.append({
                        "residuo": parts[colunas_map["_atom_site.label_comp_id"]],
                        "cadeia": cadeia,
                        "x": float(parts[colunas_map["_atom_site.Cartn_x"]]),
                        "y": float(parts[colunas_map["_atom_site.Cartn_y"]]),
                        "z": float(parts[colunas_map["_atom_site.Cartn_z"]]),
                    })
            except (KeyError, IndexError, ValueError):
                continue
    
    return mapa_cadeia_para_entidade, mapa_entidade_para_descricao, atomos

def listInteractions(variante):
    pasta = f"../data/{variante}"
    distancia_minima = 6.00
    interacoes_arquivos = defaultdict(set)

    caminho_cadeias_alvo = f"../results/cadeias_alvo/cadeias_alvo_{variante}.txt"
    if not os.path.exists(caminho_cadeias_alvo):
        print(f"Erro: Arquivo de cadeias alvo não encontrado em {caminho_cadeias_alvo}")
        return

    cadeias_alvo = {}
    with open(caminho_cadeias_alvo, "r") as f:
        linhas_arquivo = f.readlines()
    for linha in linhas_arquivo:
        parts = linha.strip().split()
        if len(parts) >= 2:
            cif_id, cadeia_alvo = parts[0], parts[1]
            cadeias_alvo[cif_id] = cadeia_alvo

    def distancia_euclidiana(a, b):
        return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

    if not os.path.isdir(pasta):
        print(f"Erro: Pasta de dados não encontrada em {pasta}")
        return

    for arquivo in os.listdir(pasta):
        if not arquivo.endswith((".cif", ".ent")):
            continue

        cif_id = arquivo.split(".")[0]
        cadeia_alvo = cadeias_alvo.get(cif_id)
        if not cadeia_alvo:
            print(f"Aviso: Não foi encontrada cadeia alvo para o arquivo {arquivo}. Pulando...")
            continue
        
        caminho = os.path.join(pasta, arquivo)
        with open(caminho, "r") as f:
            linhas = f.readlines()

        mapa_cadeia_para_entidade, mapa_entidade_para_descricao, atomos = parse_cif_data(linhas)
        
        print(f"--- Processando {arquivo} ---")
        if not mapa_cadeia_para_entidade:
             print("AVISO: Mapa Cadeia -> Entidade está VAZIO.")
        else:
             print("Mapa Cadeia -> Entidade:", mapa_cadeia_para_entidade)
        
        if not mapa_entidade_para_descricao:
             print("AVISO: Mapa Entidade -> Descrição está VAZIO.")
        else:
             print("Mapa Entidade -> Descrição:", mapa_entidade_para_descricao)
        
        residuos_espurios = (
            "HOH", "WAT", "ACE", "ACT", "BME", "CSD", "CSW", "EDO", "FMT", "GOL", 
            "MSE", "NAG", "NO3", "PO4", "SGM", "SO4", "TPO", "LU8", "U0C", "7IO", 
            "09V", "XK7"
        )
        
        atomos_alvo = [a for a in atomos if a["cadeia"] == cadeia_alvo]
        outros_atomos = [a for a in atomos if a["cadeia"] != cadeia_alvo and a["residuo"] not in residuos_espurios]

        print(f"Processando {arquivo} com {len(atomos_alvo)} átomos na cadeia alvo e {len(outros_atomos)} átomos externos")

        for a1 in atomos_alvo:
            for a2 in outros_atomos:
                d = distancia_euclidiana((a1["x"], a1["y"], a1["z"]), (a2["x"], a2["y"], a2["z"]))
                if d <= distancia_minima:
                    entity_id = mapa_cadeia_para_entidade.get(a2["cadeia"])
                    descricao = mapa_entidade_para_descricao.get(entity_id, "Descrição não encontrada")

                    identificador = f"{a2['residuo']} (Cadeia {a2['cadeia']}, {descricao})"
                    interacoes_arquivos[arquivo].add(identificador)

    os.makedirs("../results/interacoes", exist_ok=True)
    with open(f"../results/interacoes/interaction_list_{variante}.txt", "w") as file:
        for arquivo, interacoes in sorted(interacoes_arquivos.items()):
            file.write(f"Arquivo {arquivo}\n")

            if not interacoes:
                file.write("Nenhuma interacao encontrada\n")
            else:
                for i in sorted(list(interacoes)):
                    file.write(f"{i}\n")
            file.write("\n")

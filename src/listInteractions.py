import os
import math
from collections import defaultdict
import shlex

# --- Função parse_cif_data ---
def parse_cif_data(linhas):
    # inicialização
    mapa_cadeia_para_entidade = {}
    mapa_entidade_para_descricao = {}
    atomos = []

    # variáveis de controle
    colunas_map = {}
    coluna_list = []
    current_loop_type = None
    expected_columns = 0
    accumulated_parts = []

    # controle para campos multi-linha (;)
    parsing_multi_line_field = False
    multi_line_buffer = ""
    target_column_index = -1

    # processa cada linha
    for linha_num, linha in enumerate(linhas):
        linha_strip = linha.strip()

        if not linha_strip or linha_strip.startswith("#"):
            continue

        # continua lendo campo multi-linha (;)
        if parsing_multi_line_field:
            if linha_strip.endswith(";"):
                multi_line_buffer += " " + linha_strip[:-1].strip()
                current_field_content = ' '.join(multi_line_buffer.split())
                if target_column_index == len(accumulated_parts):
                     accumulated_parts.append(current_field_content)
                else:
                    accumulated_parts.append(current_field_content)
                parsing_multi_line_field = False
                multi_line_buffer = ""
                target_column_index = -1
            else:
                multi_line_buffer += " " + linha_strip
                continue # pula para próxima linha física

        # início de bloco 'loop_'
        elif linha_strip == "loop_":
            current_loop_type = None
            colunas_map.clear()
            coluna_list = []
            expected_columns = 0
            accumulated_parts = []
            parsing_multi_line_field = False
            multi_line_buffer = ""
            target_column_index = -1
            continue

        # leitura das definições de coluna (_)
        elif linha_strip.startswith("_"):
            if expected_columns > 0 and accumulated_parts:
                 accumulated_parts = [] # limpa partes órfãs se encontrar coluna após dados

            coluna_nome = linha_strip
            if current_loop_type is None:
                if "_entity." in coluna_nome and "_entity_poly." not in coluna_nome:
                    current_loop_type = "entity"
                elif "_atom_site." in coluna_nome:
                    current_loop_type = "atom_site"

            if current_loop_type in ["entity", "atom_site"]:
                 colunas_map[coluna_nome] = len(coluna_list)
                 coluna_list.append(coluna_nome)
                 expected_columns = len(coluna_list)
            continue

        # processamento das linhas de dados (se não estiver lendo campo multi-linha ;)
        elif current_loop_type in ["entity", "atom_site"] and expected_columns > 0:
            try:
                # trata início de campo multi-linha (;) que também pode terminar na mesma linha
                if linha_strip.startswith(";") and not linha_strip.endswith(";"):
                    parsing_multi_line_field = True
                    multi_line_buffer = linha_strip[1:].strip()
                    target_column_index = len(accumulated_parts)
                    continue # Vai para próxima linha para acumular o resto
                elif linha_strip.startswith(";") and linha_strip.endswith(";"):
                    current_parts = [' '.join(linha_strip[1:-1].strip().split())]
                    accumulated_parts.extend(current_parts)
                else:
                    current_parts = shlex.split(linha_strip)
                    accumulated_parts.extend(current_parts)

            except ValueError as e: # erro no shlex.split
                accumulated_parts = []
                continue

        # processa registros completos do acumulador (AGORA APÓS adicionar partes)
        while expected_columns > 0 and len(accumulated_parts) >= expected_columns:
            parts_to_process = []
            process_error = None
            try:
                parts_to_process = accumulated_parts[:expected_columns]

                if len(parts_to_process) != expected_columns:
                    raise IndexError(f"Part count mismatch. Expected {expected_columns}, got {len(parts_to_process)}.")

                # processamento bloco '_entity'
                if current_loop_type == "entity":
                    id_idx = colunas_map.get("_entity.id")
                    desc_idx = colunas_map.get("_entity.pdbx_description")

                    if id_idx is None: raise KeyError("_entity.id column missing")
                    entity_id = parts_to_process[id_idx]

                    if desc_idx is not None:
                         descricao_raw = parts_to_process[desc_idx]
                         descricao = ' '.join(descricao_raw.strip("'\"").split())
                         if descricao not in ('.', '?'):
                            mapa_entidade_para_descricao[entity_id] = descricao

                # processamento bloco '_atom_site'
                elif current_loop_type == "atom_site":
                    chain_idx = colunas_map.get("_atom_site.auth_asym_id")
                    entity_idx = colunas_map.get("_atom_site.label_entity_id")
                    residue_idx = colunas_map.get("_atom_site.label_comp_id")
                    x_idx = colunas_map.get("_atom_site.Cartn_x")
                    y_idx = colunas_map.get("_atom_site.Cartn_y")
                    z_idx = colunas_map.get("_atom_site.Cartn_z")

                    if None in [chain_idx, entity_idx, residue_idx, x_idx, y_idx, z_idx]:
                         raise KeyError("Essential atom_site column definition missing in map")

                    cadeia = parts_to_process[chain_idx] # type: ignore
                    entity_id = parts_to_process[entity_idx] # type: ignore
                    residuo = parts_to_process[residue_idx] # type: ignore
                    x_str = parts_to_process[x_idx] # type: ignore
                    y_str = parts_to_process[y_idx] # type: ignore
                    z_str = parts_to_process[z_idx] # type: ignore

                    x = float(x_str) if x_str not in ('.', '?') else None
                    y = float(y_str) if y_str not in ('.', '?') else None
                    z = float(z_str) if z_str not in ('.', '?') else None

                    if x is not None and y is not None and z is not None:
                        if cadeia != '.' and entity_id != '.' and cadeia not in mapa_cadeia_para_entidade:
                            mapa_cadeia_para_entidade[cadeia] = entity_id

                        atomos.append({
                            "residuo": residuo, "cadeia": cadeia, "x": x, "y": y, "z": z,
                        })

            except (IndexError, ValueError, KeyError) as e:
                 process_error = e
                 pass
            except Exception as e:
                 process_error = e
                 pass

            # remove as partes correspondentes ao registro (tentado) do acumulador
            accumulated_parts = accumulated_parts[expected_columns:]

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
    try:
        with open(caminho_cadeias_alvo, "r") as f:
            linhas_arquivo = f.readlines()
        for linha in linhas_arquivo:
            parts = linha.strip().split()
            if len(parts) >= 2:
                cif_id, cadeia_alvo = parts[0], parts[1]
                cadeias_alvo[cif_id] = cadeia_alvo
        print(f"Cadeias alvo carregadas: {len(cadeias_alvo)} entradas.")
    except Exception as e:
        print(f"Erro ao ler o arquivo de cadeias alvo {caminho_cadeias_alvo}: {e}")
        return

    def distancia_euclidiana(a, b):
        try:
            return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)
        except TypeError:
            return float('inf')

    if not os.path.isdir(pasta):
        print(f"Erro: Pasta de dados não encontrada em {pasta}")
        return

    # itera sobre os arquivos na pasta
    for arquivo in sorted(os.listdir(pasta)): # ordena para processamento consistente
        if not arquivo.endswith((".cif", ".ent")):
            continue

        print(f"\n--- Tentando processar: {arquivo} ---") # indica qual arquivo está sendo tentado

        cif_id = arquivo.split(".")[0]
        cadeia_alvo = cadeias_alvo.get(cif_id)
        if not cadeia_alvo:
            print(f"Aviso: Não foi encontrada cadeia alvo para {arquivo}. Pulando...")
            continue
        print(f"Cadeia alvo para {arquivo}: '{cadeia_alvo}'")

        caminho = os.path.join(pasta, arquivo)
        try:
            with open(caminho, "r") as f:
                linhas = f.readlines()
            print(f"Arquivo {arquivo} lido com sucesso ({len(linhas)} linhas).")
        except Exception as e:
            print(f"Erro ao ler o arquivo {arquivo}: {e}. Pulando...")
            continue

        # chama a função de parsing
        mapa_cadeia_para_entidade, mapa_entidade_para_descricao, atomos = parse_cif_data(linhas)

        # verificações básicas dos resultados do parsing
        if not atomos:
             print(f"AVISO: Nenhum átomo foi extraído de {arquivo} após o parsing. Verifique o arquivo ou o parser.")

        if not mapa_cadeia_para_entidade:
             print("AVISO: Mapa Cadeia -> Entidade está VAZIO.")

        if not mapa_entidade_para_descricao:
             print("AVISO: Mapa Entidade -> Descrição está VAZIO.")

        residuos_espurios = (
            "HOH", "WAT", "ACE", "ACT", "BME", "CSD", "CSW", "EDO", "FMT", "GOL",
            "MSE", "NAG", "NO3", "PO4", "SGM", "SO4", "TPO", "LU8", "U0C", "7IO",
            "09V", "XK7"
        )

        # filtra os átomos
        atomos_alvo = [a for a in atomos if a["cadeia"] == cadeia_alvo]
        outros_atomos = [a for a in atomos if a.get("cadeia") != cadeia_alvo and a.get("residuo") not in residuos_espurios]

        if not atomos_alvo:
             print(f"AVISO: Nenhum átomo encontrado para a cadeia alvo '{cadeia_alvo}' em {arquivo}.")

        print(f"Átomos na cadeia alvo: {len(atomos_alvo)}, Átomos externos (não espúrios): {len(outros_atomos)}")

        # calcula interações
        interacoes_encontradas_neste_arquivo = 0
        for a1 in atomos_alvo:
            if not all(k in a1 and isinstance(a1[k], (int, float)) for k in ["x", "y", "z"]):
                continue
            coord1 = (a1["x"], a1["y"], a1["z"])

            for a2 in outros_atomos:
                if not all(k in a2 and isinstance(a2[k], (int, float)) for k in ["x", "y", "z"]):
                    continue
                coord2 = (a2["x"], a2["y"], a2["z"])

                d = distancia_euclidiana(coord1, coord2)
                if d <= distancia_minima:
                    entity_id = mapa_cadeia_para_entidade.get(a2.get("cadeia"))
                    descricao = mapa_entidade_para_descricao.get(entity_id, "Descrição não encontrada")

                    residuo_a2 = a2.get('residuo', '?')
                    cadeia_a2 = a2.get('cadeia', '?')

                    identificador = f"{residuo_a2} (Cadeia {cadeia_a2}, {descricao})"
                    interacoes_arquivos[arquivo].add(identificador)
                    interacoes_encontradas_neste_arquivo += 1

        if interacoes_encontradas_neste_arquivo > 0:
             print(f"{interacoes_encontradas_neste_arquivo} interações encontradas para {arquivo}.")
        # Se 0 interações foram encontradas, o arquivo será registrado como "Nenhuma interação" na saída.


    # escrita do arquivo de resultados
    os.makedirs("../results/interacoes", exist_ok=True)
    caminho_saida = f"../results/interacoes/interaction_list_{variante}.txt"
    try:
        with open(caminho_saida, "w") as file:
            for arquivo, interacoes in sorted(interacoes_arquivos.items()):
                file.write(f"Arquivo {arquivo}\n")

                if not interacoes:
                    file.write("Nenhuma interacao encontrada\n")
                else:
                    for i in sorted(list(interacoes)):
                        file.write(f"{i}\n")
                file.write("\n")
        print(f"\nResultados salvos em {caminho_saida}")
    except Exception as e:
        print(f"Erro ao escrever no arquivo de saída {caminho_saida}: {e}")

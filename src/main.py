from fetchStructures import fetchStructures
from findSimilarChains import findChain
from listInteractions import listInteractions
import os

folder_drivers = "../entradas/drivers"
folder_non_drivers = "../entradas/non-drivers"

#obter sequencias homologas
fetchStructures(folder_drivers)
fetchStructures(folder_non_drivers)

for file in os.listdir(folder_drivers):
    name = file.split(".")[0]
    findChain("drivers", f"{name}") #encontrar cadeia mais similar à variante de interesse
    listInteractions(f"{name}") #listar residuos com os quais as cadeias de interesse interagem

for file in os.listdir(folder_non_drivers):
    name = file.split(".")[0]
    findChain("non-drivers", f"{name}") #encontrar cadeia mais similar à variante de interesse
    listInteractions(f"{name}") #listar residuos com os quais as cadeias de interesse interagem
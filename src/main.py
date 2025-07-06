from fetchStructures import fetchStructures
from findSimilarChains import findChain
from listInteractions import listInteractions

#obter sequencias homologas
fetchStructures("../entradas/drivers")
fetchStructures("../entradas/non-drivers")

#encontrar cadeia mais similar à variante de interesse
findChain("drivers", "ACVR1_mutated_pG328E")

#listar residuos com os quais as cadeias de interesse interagem
listInteractions("ACVR1_mutated_pG328E")
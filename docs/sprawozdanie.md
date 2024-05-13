# Position SBH Err+ - Sekwencjonowanie z informacją o położeniu, z błędami pozytywnymi
Mateusz Chlebosz 151817, Jakub Aszyk 151841

## 1. Wstęp
Celem laboratorium jest zaimplementowanie algorytmu pozwalającego na sekwencjonowanie z informacją o położeniu, z błędami pozytywnymi. W ramach laboratorium zaimplementowano algorytm Position SBH Err+.


1. Sortujemy po PosL i PosH
2. Hashmapa inkrementalna
3.

Zysk = Ilość odwiedzonych wierzchołków
Wierzchołek = zysk + 1

<!-- Zysk łuku – zależny od odległości pomiędzy wierzchołkami????? -->

Łuki łączą się tylko jeżeli ich pozycje są blisko siebie $h_i+1 ≥ l_{i+1}$

np. PosL=1 PosH=10 może się połączyć z -> PosL=11 PosH=20
ale nie z PosL=15 PosH=25

Koszt = (Łuk) różnica w dodatkowych Zasadach azotowych nukleotydów np ACTG -> TGCA koszt = 2
Koszt <= Długości sekwencji dna

Max (Zysk)


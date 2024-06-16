# Position SBH Err+ - Sekwencjonowanie z informacją o położeniu, z błędami pozytywnymi
Mateusz Chlebosz 151817, Jakub Aszyk 151841

## 1. Opis problemu

Celem laboratorium jest zaimplementowanie algorytmu pozwalającego na sekwencjonowanie z informacją o położeniu, z błędami pozytywnymi. W ramach laboratorium zaimplementowano algorytm Position SBH Err+.

Sekwencjonowanie to proces określania sekwencji nukleotydów w DNA. Jest to proces określenia kolejności nukleotydów w DNA. Sekwencjonowanie DNA jest kluczowym narzędziem w dziedzinie biologii molekularnej, pozwala na zrozumienie struktury genów, ich funkcji oraz ewolucji.

SBH (Sequencing by Hybridization) to technika sekwencjonowania DNA, która polega na hybrydyzacji krótkich sond DNA z DNA, które chcemy sekwencjonować. W wyniku hybrydyzacji otrzymujemy informację o sekwencji DNA.

Hybrydyzacja to proces łączenia dwóch nici DNA, które są komplementarne do siebie. W wyniku hybrydyzacji otrzymujemy informację o sekwencji DNA.

PSBH (Position SBH) to technika sekwencjonowania DNA, która dodatkowo zawiera informację o położeniu sond DNA na DNA, które chcemy sekwencjonować.

Oligonukleotydy to krótkie sondy DNA, które są używane w technikach sekwencjonowania DNA. Oligonukleotydy są używane do hybrydyzacji z DNA, które chcemy sekwencjonować.


### Problem sekwencjonowania z informacją o położeniu

Oprócz sekwencji DNA, na wejściu otrzymujemy informację o położeniu sond DNA na DNA, które chcemy sekwencjonować. Pomaga to w zrekonstruowaniu sekwencji DNA.

PosH, PosL - pozycje na DNA

Jak to wygląda, fragment danych wejściowych, rekonstrukcja kolejności

### Błędy pozytywne

Błędy pozytywne to sytuacja, w której w wyniku sekwencjonowania otrzymujemy więcej informacji niż jest w rzeczywistości. W wyniku błędów pozytywnych otrzymujemy dodatkowe zasady azotowe nukleotydów. Błędy pozytywne mogą powstać w wyniku błędów pomiarowych. Musimy więc uwzględnić błędy pozytywne w naszym algorytmie. i założyć możliwość pominięcia oligonukleotydów.


## 2. Opis algorytmu dokładnego

### Algorytm Position SBH Err+ - STS

Selective Traveling Salesman
Każdy wierzchołek ma zysk a łuk koszt
Całkowity koszt ścieżki nie przekroczy ustalonej wartości

Możemy tylko przechodzić do łuków, które są blisko siebie

$h_i+1 ≥ l_{j}$ i $l_i <= h_{j}$

dodatkowo ze względu na odległość między wierzchołkami

$h_i+1 ≥ l_{j}$ i $l_i + diff(j,i) <= h_{j} + 1$


gdzie diff(j,i) to różnica w dodatkowych Zasadach azotowych nukleotydów np
ACTG -> TGCA koszt = 2


### Implementacja

#### Dokładna

Selected Traveling Salesman Comivoyeur

#### Przyblizona

Algorytm Zachłanny

### Podsumowanie
Analizę uzyskanych wyników


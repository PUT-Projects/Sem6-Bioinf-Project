# Position SBH Err+ - Sekwencjonowanie z informacją o położeniu, z błędami pozytywnymi
Mateusz Chlebosz 151817, Jakub Aszyk 151841

## Opis problemu
Problem dotyczy analizy sekwencji DNA w celu znalezienia optymalnej ścieżki między różnymi segmentami DNA. Sekwencje te są reprezentowane w formacie XML, a celem jest stworzenie algorytmu, który efektywnie znajdzie najlepszą ścieżkę przejścia przez te segmenty, minimalizując różnice między kolejnymi segmentami.

## Sformalizowanie problemu
Problem można sformalizować jako problem znajdowania najkrótszej ścieżki (pathfinding) w grafie, gdzie węzłami są segmenty DNA, a krawędziami są różnice między kolejnymi segmentami. Celem jest maksymalizacja liczby odwiedzonych węzłów, przy jednoczesnym minimalizowaniu całkowitego kosztu przejścia, który nie może przekroczyć określonej długości.


## Opis algorytmu dokładnego

### Algorytm Position SBH Err+ - STS

Użyty został Selective Traveling Salesman gdzie
Każdy wierzchołek ma zysk a łuk koszt
Całkowity koszt ścieżki nie może przekroczyć ustalonej wartości
Możemy tylko przechodzić do łuków, które są blisko siebie

$h_i+1 ≥ l_{j}$ i $l_i <= h_{j}$

dodatkowo ze względu na odległość między wierzchołkami

$h_i+1 ≥ l_{j}$ i $l_i + diff(j,i) <= h_{j} + 1$


gdzie diff(j,i) to różnica w dodatkowych Zasadach azotowych nukleotydów np
ACTG -> TGCA koszt = 2

### Kod
```python
class PathFinder:
    def __init__(self, graph, diff_matrix, start_id, length):
        self.__graph = graph
        self.__diff_matrix = diff_matrix
        self.__length = length
        self.__start_id = start_id
        self.best_path = []
        self.max_vertices = 0

    def run(self):
        # Start the search from the specified start_id
        visited = [False] * len(self.__graph)
        self.__backtrack(self.__start_id, [], visited, 0)
        return self.best_path

    def __backtrack(self, current_vertex, current_path, visited, current_cost):
        # Include current vertex in path
        visited[current_vertex] = True
        current_path.append(current_vertex)

        # Check if the current cost equals the desired length
        if current_cost == self.__length:
            # Check if the current path is longer (in terms of vertices) than the previously found paths
            if len(current_path) > self.max_vertices:
                self.max_vertices = len(current_path)
                self.best_path = current_path.copy()
            # Since we found a perfect cost match, we can return here for optimization
            visited[current_vertex] = False
            current_path.pop()
            return

        # Explore adjacent vertices
        for neighbor in self.__graph[current_vertex]:
            if not visited[neighbor]:
                edge_cost = self.__diff_matrix[neighbor][current_vertex]
                if current_cost + edge_cost <= self.__length:  # Only consider if within total cost
                    self.__backtrack(neighbor, current_path, visited, current_cost + edge_cost)

        # Backtrack
        visited[current_vertex] = False
        current_path.pop()

```

## Opis algorytmu przybliżonego
Algorytm wielomianowy opiera się na heurystycznym podejściu do problemu, gdzie zamiast przeszukiwać wszystkie możliwe ścieżki, algorytm korzysta z metody zachłannej, aby w każdej iteracji wybierać najkorzystniejszą opcję. Takie podejście nie gwarantuje znalezienia optymalnego rozwiązania, ale znacznie przyspiesza czas wykonania.

Przybliżony algorytm (solverHeuristick.py) działa na zasadzie przeszukiwania grafu, gdzie każdy węzeł reprezentuje segment DNA, a krawędzie reprezentują różnice między segmentami. Algorytm ten:

Parsuje dane wejściowe z pliku XML (dnaParser.py).
Tworzy graf, gdzie węzłami są segmenty DNA, a krawędziami różnice między segmentami (obliczane w offset.py).
Wykorzystuje strategię zachłanną do przeszukiwania grafu, starając się maksymalizować liczbę odwiedzonych węzłów, jednocześnie minimalizując koszt przejścia.

### Kod
```python
class PathFinderGreedy:
    def __init__(self, graph, diff_matrix, start_id, length):
        self.__graph = graph
        self.__diff_matrix = diff_matrix
        self.__length = length
        self.__start_id = start_id
        self.best_path = []

    def run(self):
        visited = [False] * len(self.__graph)
        self.__greedy_search(self.__start_id, visited)
        return self.best_path

    def __greedy_search(self, current_vertex, visited):
        current_path = []
        current_cost = 0
        while True:
            visited[current_vertex] = True
            current_path.append(current_vertex)
            self.best_path = current_path.copy()

            min_cost = float('inf')
            next_vertex = None

            for neighbor in self.__graph[current_vertex]:
                if not visited[neighbor]:
                    edge_cost = self.__diff_matrix[neighbor][current_vertex]
                    if current_cost + edge_cost <= self.__length and edge_cost < min_cost:
                        min_cost = edge_cost
                        next_vertex = neighbor

            if next_vertex is None:
                break

            current_vertex = next_vertex
            current_cost += min_cost
```


## Analiza wyników

Algorytm został przetestowany na danych wejściowych z witryny internatowej z danymi wejściowymi. Wyniki pokazują, że algorytm zachłanny jest w stanie efektywnie znaleźć ścieżkę między segmentami DNA nawet dla 1000 oligonukleotydów o długości 10. Algorytm zachłanny znacznie spowalnia swoją pracę już przy 30 oligonukleotydach.

### Platforma testowa
- Procesor: AMD Ryzen 7 7840HS
- Pamięć RAM: 32GB 6400MHz


### Wyniki
$n$ - długość sekwencji DNA, Zakres $[16, 65536]$

$k$ - długość sondy DNA, Zakres $[4, 10]$

$sqpep$ - Liczba błędów pozytywnych procentowo, wartość ≤ $25\%$

$posep$ - Szerokość przedziału ufności dla informacji o położeniu, dla starego trybu (position=='old') liczba błędów w informacji o położeniach, wartość ≤ 50%

Jeżeli nie zostało napisane inaczej używane są wartości domyślne:

Dla algorytmu dokładnego:
$$ n = 25, k = 8, sqpep = 10\%, posep = 10\% $$

Dla algorytmu przybliżonego:
$$ n = 100,  k = 8, sqpep = 10\%, posep = 10\% $$

### Dla zmiennej n

#### Algorytm Dokładny

| n   | czas [ms] | optimum |
| --- | --------- | ------- |
| 20  | 0.0       | True    |
| 25  | 5.5       | True    |
| 30  | 74.53     | True    |
| 35  | 1723.9    | True    |


#### Algorytm Przybliżony

| n   | czas [ms] | optimum |
| --- | --------- | ------- |
| 25  | 1.01      | True    |
| 20  | 1.01      | True    |
| 30  | 1.0       | True    |
| 40  | 1.0       | True    |
| 50  | 2.02      | True    |
| 100 | 10.98     | True    |
| 500 | 301.29    | True    |



### Dla zmiennej k


#### Algorytm Dokładny


| k   | czas [ms] | optimum |
| --- | --------- | ------- |
| 4   | 95.7      | True    |
| 5   | 29.54     | True    |
| 6   | 17.76     | True    |
| 7   | 12.51     | True    |
| 8   | 5.45      | True    |
| 9   | 3.0       | True    |
| 10  | 2.01      | True    |


#### Algorytm Przybliżony

| k   | czas [ms] | optimum |
| --- | --------- | ------- |
| 6   | 7.99      | True    |
| 7   | 10.0      | True    |
| 8   | 10.0      | True    |
| 9   | 11.0      | True    |
| 10  | 11.0      | True    |

### Dla zmiennej sqpep


#### Algorytm Dokładny

| sqpep [%] | time [ms] | optimum |
| --------- | --------- | ------- |
| 5         | 6.54      | True    |
| 10        | 3.0       | True    |
| 15        | 3.96      | True    |
| 20        | 5.03      | True    |
| 25        | 5.0       | True    |


#### Algorytm Przybliżony

| sqpep [%] | time [ms] | optimum |
| --------- | --------- | ------- |
| 5         | 10.0      | True    |
| 10        | 10.49     | True    |
| 15        | 10.97     | True    |
| 20        | 12.61     | True    |
| 25        | 12.91     | True    |


### Dla zmiennej posep

#### Algorytm Dokładny

| posep [%] | time [ms] | optimum |
| --------- | --------- | ------- |
| 5         | 0.0       | True    |
| 10        | 5.89      | True    |
| 15        | 4.51      | True    |
| 20        | 10.18     | True    |
| 25        | 15.6      | True    |
| 30        | 19.59     | True    |
| 40        | 33.6      | True    |
| 50        | 63.36     | True    |




#### Algorytm Przybliżony

| posep [%] | time [ms] | optimum |
| --------- | --------- | ------- |
| 5         | 9.61      | True    |
| 10        | 10.0      | True    |
| 15        | 11.0      | True    |
| 20        | 11.0      | True    |
| 25        | 10.0      | True    |
| 30        | 11.0      | True    |
| 40        | 11.0      | True    |
| 50        | 12.03     | True    |



## Wnioski

Algorytm dokładny jest w stanie znaleźć optymalną ścieżkę między segmentami DNA, ale jest znacznie wolniejszy niż algorytm przybliżony. Algorytm przybliżony, mimo że nie gwarantuje optymalnego rozwiązania, jest w stanie znaleźć satysfakcjonujące rozwiązanie w krótszym czasie. Dla większych danych wejściowych zaleca się użycie algorytmu przybliżonego ze względu na krótszy czas wykonania.

Co ciekawe algorytm przybliżony wykonuje się szybciej dla większego k niż dla mniejszego k. Dla algorytmu dokładnego jest odwrotnie, im większe k tym dłuższy czas wykonania.
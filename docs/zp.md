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

### Dla zmiennej n

$$ k = 6, sqpep = 10\%, posep = 10\% $$

| n    | Dokładny             | Przybliżony            |
| ---- | -------------------- | ---------------------- |
| 20   | 0.0239961s           | 0.005000114440917969 s |
| 25   | 0.09800052642822266s | 0.0059967041015625 s   |
| 30   | 3.79111s             | 0.026000 s             |
| 40   | -                    | 0.0100 s               |
| 50   | -                    | 0.11 s                 |
| 100  | -                    | 0.11 s                 |
| 500  | -                    | 0.11 s                 |
| 1000 | -                    | 1.3 s                  |
| 2000 | -                    | 1.3 s                  |
| 5000 | -                    | 1.3 s                  |

### Dla zmiennej k

$$ n = 25, sqpep = 5, posep = 5 $$


### Dla zmiennej sqpep

$$ n = 25, k = 6, posep = 5 $$

| sqpep [%] | Dokładny             | Przybliżony          |
| --------- | -------------------- | -------------------- |
| 5         | 0.09800052642822266s | 0.0059967041015625 s |
| 10        | 0.09800052642822266s | 0.0059967041015625 s |
| 15        | 0.09800052642822266s | 0.0059967041015625 s |
| 20        | 0.09800052642822266s | 0.0059967041015625 s |
| 25        | 0.09800052642822266s | 0.0059967041015625 s |

### Dla zmiennej posep

$$ n = 25, k = 6, sqpep = 5 $$

| posep [%] | Dokładny             | Przybliżony          |
| --------- | -------------------- | -------------------- |
| 5         | 0.09800052642822266s | 0.0059967041015625 s |
| 10        | 0.09800052642822266s | 0.0059967041015625 s |
| 15        | 0.09800052642822266s | 0.0059967041015625 s |
| 20        | 0.09800052642822266s | 0.0059967041015625 s |
| 25        | 0.09800052642822266s | 0.0059967041015625 s |
| 30        | 0.09800052642822266s | 0.0059967041015625 s |
| 40        | 0.09800052642822266s | 0.0059967041015625 s |
| 50        | 0.09800052642822266s | 0.0059967041015625 s |
###

| n    | k   | sqpe | pose | Dokładny             | Przybliżony            |
| ---- | --- | ---- | ---- | -------------------- | ---------------------- |
| 20   | 5   | 4    | 4    | 0.0239961s           | 0.005000114440917969 s |
| 25   | 6   | 5    | 5    | 0.09800052642822266s | 0.0059967041015625 s   |
| 30   | 6   | 5    | 5    | 3.79111s             | 0.026000 s             |
| 40   | 6   | 5    | 5    | -                    | 0.0100 s               |
| 50   | 8   | 5    | 5    | -                    | 0.11 s                 |
| 100  | 8   | 5    | 5    |
| 1000 | 10  | 5    | 5    | -                    | 1.3 s                  |




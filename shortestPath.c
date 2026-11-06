
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_VERTICES 500
#define MAX_EDGES (MAX_VERTICES * (MAX_VERTICES - 1))
#define INF (LLONG_MAX / 4)

typedef long long Distance;

typedef struct {
    int from;
    int to;
} Edge;

static Distance graph[MAX_VERTICES][MAX_VERTICES];
static Distance floydSerialResult[MAX_VERTICES][MAX_VERTICES];
static Distance floydOpenMPResult[MAX_VERTICES][MAX_VERTICES];
static Distance dijkstraSerialResult[MAX_VERTICES][MAX_VERTICES];
static Distance dijkstraOpenMPResult[MAX_VERTICES][MAX_VERTICES];

static Edge possibleEdges[MAX_EDGES];

/* Compile with OpenMP enabled, for example:
 * gcc -O3 -fopenmp shortestPath.c -o shortest_paths.exe
 */
static double wallTime(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}

void generateRandomGraph(int vertexCount, int edgeCount) {
    int possibleEdgeCount = 0;

    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            graph[i][j] = (i == j) ? 0 : INF;
        }
    }

    for (int from = 0; from < vertexCount; from++) {
        for (int to = 0; to < vertexCount; to++) {
            if (from != to) {
                possibleEdges[possibleEdgeCount].from = from;
                possibleEdges[possibleEdgeCount].to = to;
                possibleEdgeCount++;
            }
        }
    }

    for (int i = possibleEdgeCount - 1; i > 0; i--) {
        int randomIndex = rand() % (i + 1);

        Edge temporary = possibleEdges[i];
        possibleEdges[i] = possibleEdges[randomIndex];
        possibleEdges[randomIndex] = temporary;
    }

    for (int i = 0; i < edgeCount; i++) {
        int from = possibleEdges[i].from;
        int to = possibleEdges[i].to;
        int weight = 1 + rand() % 100;

        graph[from][to] = weight;
    }
}

void floydWarshallBasic(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            distance[i][j] = graph[i][j];
        }
    }

    for (int k = 0; k < vertexCount; k++) {
        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < vertexCount; j++) {
                if (distance[i][k] != INF &&
                    distance[k][j] != INF) {

                    Distance candidate =
                        distance[i][k] + distance[k][j];

                    if (candidate < distance[i][j]) {
                        distance[i][j] = candidate;
                    }
                }
            }
        }
    }
}

void floydWarshallOpenMP(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            distance[i][j] = graph[i][j];
        }
    }

    for (int k = 0; k < vertexCount; k++) {
        /* k is sequential, but all rows are independent for a fixed k. */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < vertexCount; i++) {
            Distance distanceToK = distance[i][k];

            if (distanceToK == INF) {
                continue;
            }

            for (int j = 0; j < vertexCount; j++) {
                if (distance[k][j] != INF) {
                    Distance candidate = distanceToK + distance[k][j];

                    if (candidate < distance[i][j]) {
                        distance[i][j] = candidate;
                    }
                }
            }
        }
    }
}


static void dijkstraSingleSource(
    int vertexCount,
    int source,
    Distance result[MAX_VERTICES]
) {
    Distance currentDistance[MAX_VERTICES];
    int visited[MAX_VERTICES];

    for (int vertex = 0; vertex < vertexCount; vertex++) {
        currentDistance[vertex] = INF;
        visited[vertex] = 0;
    }

    currentDistance[source] = 0;

    for (int iteration = 0; iteration < vertexCount; iteration++) {
        int currentVertex = -1;
        Distance smallestDistance = INF;

        for (int vertex = 0; vertex < vertexCount; vertex++) {
            if (!visited[vertex] &&
                currentDistance[vertex] < smallestDistance) {

                smallestDistance = currentDistance[vertex];
                currentVertex = vertex;
            }
        }

        if (currentVertex == -1) {
            break;
        }

        visited[currentVertex] = 1;

        for (int neighbor = 0; neighbor < vertexCount; neighbor++) {
            if (graph[currentVertex][neighbor] == INF) {
                continue;
            }

            Distance candidate = currentDistance[currentVertex] +
                                 graph[currentVertex][neighbor];

            if (candidate < currentDistance[neighbor]) {
                currentDistance[neighbor] = candidate;
            }
        }
    }

    for (int vertex = 0; vertex < vertexCount; vertex++) {
        result[vertex] = currentDistance[vertex];
    }
}




void dijkstraAllPairs(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    for (int source = 0; source < vertexCount; source++) {
        dijkstraSingleSource(vertexCount, source, distance[source]);
    }
}

void dijkstraAllPairsOpenMP(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    /* Every source has an independent Dijkstra calculation. */
    #pragma omp parallel for schedule(static)
    for (int source = 0; source < vertexCount; source++) {
        dijkstraSingleSource(vertexCount, source, distance[source]);
    }
}


int matricesEqual(
    int vertexCount,
    Distance first[MAX_VERTICES][MAX_VERTICES],
    Distance second[MAX_VERTICES][MAX_VERTICES]
) {
    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            if (first[i][j] != second[i][j]) {
                return 0;
            }
        }
    }

    return 1;
}

void printDistances(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    printf("\nAll-pairs shortest-distance matrix:\n");

    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            if (distance[i][j] == INF) {
                printf("%8s", "INF");
            } else {
                printf("%8lld", distance[i][j]);
            }
        }
        printf("\n");
    }
}

void printGraph(int vertexCount) {
    printf("\nGenerated graph:\n");
    printf("From -> To (Weight)\n");

    for (int from = 0; from < vertexCount; from++) {
        for (int to = 0; to < vertexCount; to++) {
            if (from != to && graph[from][to] != INF) {
                printf("%d -> %d (weight: %lld)\n",
                       from, to, graph[from][to]);
            }
        }
    }
}

int main(void) {
    int vertexCount;
    int edgeCount;
    int printGraphBool = 0;
    srand((unsigned int)time(NULL));

    printf("Enter number of vertices (1-%d): ", MAX_VERTICES);

    if (scanf("%d", &vertexCount) != 1 ||
        vertexCount < 1 ||
        vertexCount > MAX_VERTICES) {

        printf("Invalid number of vertices.\n");
        return 1;
    }

    int maximumEdgeCount = vertexCount * (vertexCount - 1);

    printf("Enter number of directed edges (0-%d): ",
           maximumEdgeCount);

    if (scanf("%d", &edgeCount) != 1 ||
        edgeCount < 0 ||
        edgeCount > maximumEdgeCount) {

        printf("Invalid number of edges.\n");
        return 1;
    }

    generateRandomGraph(vertexCount, edgeCount);

    printf("\nRandom graph generated successfully.\n");
    printf("Vertices: %d\n", vertexCount);
    printf("Directed edges: %d\n", edgeCount);
    printf("Edge weights: 1-100\n");

    printf("Do you want to print the graph? (1 for yes, 0 for no): ");
    scanf("%d", &printGraphBool);
    if(printGraphBool) {
        printGraph(vertexCount);
    }
    

    

#ifdef _OPENMP
    int availableProcessors = omp_get_num_procs();
    int threadCount;

    printf("Enter number of OpenMP threads (1-%d): ",
           availableProcessors);

    if (scanf("%d", &threadCount) != 1 ||
        threadCount < 1 ||
        threadCount > availableProcessors) {

        printf("Invalid number of OpenMP threads.\n");
        return 1;
    }

    omp_set_num_threads(threadCount);
    printf("OpenMP threads selected: %d\n", threadCount);
#else
    printf("\nOpenMP is not enabled; only serial execution is available.\n");
#endif

    double start = wallTime();
    floydWarshallBasic(vertexCount, floydSerialResult);
    double floydSerialTime = wallTime() - start;

    start = wallTime();
    floydWarshallOpenMP(vertexCount, floydOpenMPResult);
    double floydOpenMPTime = wallTime() - start;

    start = wallTime();
    dijkstraAllPairs(vertexCount, dijkstraSerialResult);
    double dijkstraSerialTime = wallTime() - start;

    start = wallTime();
    dijkstraAllPairsOpenMP(vertexCount, dijkstraOpenMPResult);
    double dijkstraOpenMPTime = wallTime() - start;

    printf("\n--- Performance comparison ---\n");
    printf("Serial Floyd-Warshall:    %.6f seconds\n",
           floydSerialTime);
    printf("OpenMP Floyd-Warshall:    %.6f seconds\n",
           floydOpenMPTime);
    printf("Serial repeated Dijkstra: %.6f seconds\n",
           dijkstraSerialTime);
    printf("OpenMP repeated Dijkstra: %.6f seconds\n",
           dijkstraOpenMPTime);

    printf("\n--- Correctness checks ---\n");
    printf("Serial FW vs. OpenMP FW:  %s\n",
           matricesEqual(vertexCount, floydSerialResult,
                         floydOpenMPResult) ? "PASS" : "FAIL");
    printf("Serial Dijkstra vs. OpenMP Dijkstra: %s\n",
           matricesEqual(vertexCount, dijkstraSerialResult,
                         dijkstraOpenMPResult) ? "PASS" : "FAIL");
    printf("Serial FW vs. Serial Dijkstra: %s\n",
           matricesEqual(vertexCount, floydSerialResult,
                         dijkstraSerialResult) ? "PASS" : "FAIL");
    printf("OpenMP FW vs. OpenMP Dijkstra: %s\n",
           matricesEqual(vertexCount, floydOpenMPResult,
                         dijkstraOpenMPResult) ? "PASS" : "FAIL");
    return 0;
}

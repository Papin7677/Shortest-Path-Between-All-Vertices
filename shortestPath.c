
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define MAX_VERTICES 500
#define MAX_EDGES (MAX_VERTICES * (MAX_VERTICES - 1))
#define INF (LLONG_MAX / 4)

typedef long long Distance;

typedef struct {
    int from;
    int to;
} Edge;

static Distance graph[MAX_VERTICES][MAX_VERTICES];
static Distance basicResult[MAX_VERTICES][MAX_VERTICES];
static Distance dijkstraResult[MAX_VERTICES][MAX_VERTICES];

static Edge possibleEdges[MAX_EDGES];

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




void dijkstraAllPairs(
    int vertexCount,
    Distance distance[MAX_VERTICES][MAX_VERTICES]
) {
    Distance currentDistance[MAX_VERTICES];
    int visited[MAX_VERTICES];

    for (int source = 0; source < vertexCount; source++) {
        for (int vertex = 0; vertex < vertexCount; vertex++) {
            currentDistance[vertex] = INF;
            visited[vertex] = 0;
        }

        currentDistance[source] = 0;

        for (int iteration = 0;
             iteration < vertexCount;
             iteration++) {

            int currentVertex = -1;
            Distance smallestDistance = INF;

            for (int vertex = 0;
                 vertex < vertexCount;
                 vertex++) {

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

            for (int neighbor = 0;
                 neighbor < vertexCount;
                 neighbor++) {

                if (graph[currentVertex][neighbor] == INF) {
                    continue;
                }

                Distance candidate =
                    currentDistance[currentVertex] +
                    graph[currentVertex][neighbor];

                if (candidate < currentDistance[neighbor]) {
                    currentDistance[neighbor] = candidate;
                }
            }
        }

        for (int vertex = 0; vertex < vertexCount; vertex++) {
            distance[source][vertex] = currentDistance[vertex];
        }
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
    

    

    clock_t start;
    clock_t finish;

    start = clock();
    floydWarshallBasic(vertexCount, basicResult);
    finish = clock();
    double basicTime =
        (double)(finish - start) / CLOCKS_PER_SEC;

    start = clock();
    

    start = clock();
    dijkstraAllPairs(vertexCount, dijkstraResult);
    finish = clock();
    double dijkstraTime =
        (double)(finish - start) / CLOCKS_PER_SEC;

    printf("\n--- Performance comparison ---\n");
    printf("Basic Floyd-Warshall:     %.6f seconds\n",
           basicTime);
 
    printf("Repeated Dijkstra:        %.6f seconds\n",
           dijkstraTime);

    printf("\n--- Correctness checks ---\n");

    printf("Basic FW vs. Dijkstra:     %s\n",
           matricesEqual(vertexCount, basicResult,
                         dijkstraResult) ? "PASS" : "FAIL");
    return 0;
}

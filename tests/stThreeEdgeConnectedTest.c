/*
 * stThreeEdgeConnectedTest.c
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stThreeEdgeConnected.h"
#include <stdlib.h>

/*
 * Gets the components of the graph with the given edges, as pairs of vertices, laid out in
 * compressed sparse row form with each vertex's edges in the order they are given, and checks
 * the output holds each vertex exactly once.
 */
static void getComponents(CuTest *testCase, int n, int edgeNumber, const int *edges, int **members, int **starts, int *componentNumber) {
    int *offsets = st_calloc(n + 1, sizeof(int));
    for (int i = 0; i < 2 * edgeNumber; i++) {
        offsets[edges[i] + 1]++;
    }
    for (int v = 0; v < n; v++) {
        offsets[v + 1] += offsets[v];
    }
    int *fill = st_malloc((n + 1) * sizeof(int));
    for (int v = 0; v <= n; v++) {
        fill[v] = offsets[v];
    }
    int *adj = st_malloc((2 * edgeNumber + 1) * sizeof(int));
    for (int i = 0; i < edgeNumber; i++) {
        adj[fill[edges[2 * i]]++] = edges[2 * i + 1];
        adj[fill[edges[2 * i + 1]]++] = edges[2 * i];
    }
    stThreeEdgeConnected_getComponents(n, offsets, adj, members, starts, componentNumber);
    CuAssertIntEquals(testCase, 0, (*starts)[0]);
    CuAssertIntEquals(testCase, n, (*starts)[*componentNumber]);
    int *seen = st_calloc(n + 1, sizeof(int));
    for (int i = 0; i < *componentNumber; i++) {
        CuAssertTrue(testCase, (*starts)[i] < (*starts)[i + 1]);
        for (int j = (*starts)[i]; j < (*starts)[i + 1]; j++) {
            CuAssertTrue(testCase, (*members)[j] >= 0 && (*members)[j] < n);
            CuAssertIntEquals(testCase, 0, seen[(*members)[j]]++);
        }
    }
    free(seen);
    free(offsets);
    free(fill);
    free(adj);
}

static int getComponentNumber(CuTest *testCase, int n, int edgeNumber, const int *edges) {
    int *members, *starts, componentNumber;
    getComponents(testCase, n, edgeNumber, edges, &members, &starts, &componentNumber);
    free(members);
    free(starts);
    return componentNumber;
}

static void testStThreeEdgeConnected_examples(CuTest *testCase) {
    CuAssertIntEquals(testCase, 0, getComponentNumber(testCase, 0, 0, NULL));
    CuAssertIntEquals(testCase, 3, getComponentNumber(testCase, 3, 0, NULL));
    int selfLoops[] = { 0, 0, 0, 0 };
    CuAssertIntEquals(testCase, 1, getComponentNumber(testCase, 1, 2, selfLoops));
    //Two vertices are joined in a component by three parallel edges but not by two
    int parallel[] = { 0, 1, 1, 0, 0, 1 };
    CuAssertIntEquals(testCase, 2, getComponentNumber(testCase, 2, 2, parallel));
    CuAssertIntEquals(testCase, 1, getComponentNumber(testCase, 2, 3, parallel));
    //Removing any two edges of a cycle cuts it
    int triangle[] = { 0, 1, 1, 2, 2, 0 };
    CuAssertIntEquals(testCase, 3, getComponentNumber(testCase, 3, 3, triangle));
    int k4[] = { 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3 };
    CuAssertIntEquals(testCase, 1, getComponentNumber(testCase, 4, 6, k4));
    //Two K4s joined by a bridge, with a pendant vertex
    int bridged[] = { 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 3, 4, 4, 5, 4, 6, 4, 7, 5, 6, 5, 7, 6, 7, 7, 8 };
    CuAssertIntEquals(testCase, 3, getComponentNumber(testCase, 9, 14, bridged));
}

static void checkOrder(CuTest *testCase, int n, int edgeNumber, const int *edges, int componentNumber, const int *members, const int *starts) {
    int *m, *s, k;
    getComponents(testCase, n, edgeNumber, edges, &m, &s, &k);
    CuAssertIntEquals(testCase, componentNumber, k);
    for (int i = 0; i < n; i++) {
        CuAssertIntEquals(testCase, members[i], m[i]);
    }
    for (int i = 0; i <= k; i++) {
        CuAssertIntEquals(testCase, starts[i], s[i]);
    }
    free(m);
    free(s);
}

/*
 * The order of the components, and of the vertices in each, decides how stCactusGraph_collapseToCactus
 * merges nodes, which the cactus graphs that callers get out depend on. These are the orders the
 * previous implementation gave.
 */
static void testStThreeEdgeConnected_order(CuTest *testCase) {
    int edges1[] = { 6, 7, 1, 1, 4, 6, 3, 5, 4, 6, 0, 1, 6, 8, 2, 5, 3, 4, 0, 1, 4, 1, 5, 7, 5, 4 };
    int members1[] = { 8, 7, 2, 3, 4, 6, 5, 1, 0 };
    int starts1[] = { 0, 1, 2, 3, 4, 7, 8, 9 };
    checkOrder(testCase, 9, 13, edges1, 7, members1, starts1);

    int edges2[] = { 3, 3, 0, 1, 6, 8, 2, 2, 7, 0, 0, 2, 1, 3, 4, 4, 1, 7, 8, 1, 8, 1, 1, 2, 4, 8, 5, 6, 2, 4 };
    int members2[] = { 7, 3, 5, 6, 4, 0, 1, 8, 2 };
    int starts2[] = { 0, 1, 2, 3, 4, 5, 9 };
    checkOrder(testCase, 9, 15, edges2, 6, members2, starts2);

    int edges3[] = { 1, 0, 4, 4, 4, 5, 7, 4, 6, 0, 2, 3, 7, 0, 0, 1, 7, 1, 1, 6, 1, 4, 2, 4, 2, 3, 4, 6 };
    int members3[] = { 3, 2, 5, 0, 6, 7, 4, 1 };
    int starts3[] = { 0, 1, 2, 3, 8 };
    checkOrder(testCase, 8, 14, edges3, 4, members3, starts3);

    //Several connected components, with the vertices of each interleaved
    int edges4[] = { 2, 0, 8, 8, 8, 10, 14, 8, 12, 0, 4, 6, 14, 0, 0, 2, 14, 2, 2, 12, 2, 8, 4, 8, 4, 6, 8, 12,
                     1, 3, 1, 5, 1, 7, 3, 5, 3, 7, 5, 7, 7, 9, 9, 9, 11, 5, 11, 3 };
    int members4[] = { 6, 4, 10, 0, 12, 14, 8, 2, 9, 11, 1, 3, 5, 7, 13, 15 };
    int starts4[] = { 0, 1, 2, 3, 8, 9, 10, 14, 15, 16 };
    checkOrder(testCase, 16, 24, edges4, 9, members4, starts4);
}

static int findRoot(int *parent, int v) {
    while (parent[v] != v) {
        v = parent[v] = parent[parent[v]];
    }
    return v;
}

/*
 * Two vertices are in the same 3-edge-connected component if no cut of at most two edges separates
 * them, so check every such cut.
 */
static bool *getSameComponentByBruteForce(int n, int edgeNumber, const int *edges) {
    bool *same = st_malloc(n * n * sizeof(bool));
    for (int i = 0; i < n * n; i++) {
        same[i] = 1;
    }
    int *parent = st_malloc(n * sizeof(int));
    for (int cut1 = -1; cut1 < edgeNumber; cut1++) {
        for (int cut2 = cut1 + 1; cut2 <= edgeNumber; cut2++) { //cut2 == edgeNumber cuts just cut1
            for (int v = 0; v < n; v++) {
                parent[v] = v;
            }
            for (int i = 0; i < edgeNumber; i++) {
                if (i != cut1 && i != cut2) {
                    parent[findRoot(parent, edges[2 * i])] = findRoot(parent, edges[2 * i + 1]);
                }
            }
            for (int u = 0; u < n; u++) {
                for (int v = 0; v < n; v++) {
                    if (findRoot(parent, u) != findRoot(parent, v)) {
                        same[u * n + v] = 0;
                    }
                }
            }
        }
    }
    free(parent);
    return same;
}

static void testStThreeEdgeConnected_random(CuTest *testCase) {
    for (int test = 0; test < 1000; test++) {
        int n = st_randomInt(1, 12);
        int edgeNumber = st_randomInt(0, 3 * n);
        int *edges = st_malloc((2 * edgeNumber + 1) * sizeof(int));
        for (int i = 0; i < edgeNumber; i++) {
            edges[2 * i] = st_randomInt(0, n);
            //Mostly edges to the next vertex, to make cycles
            edges[2 * i + 1] = st_random() < 0.5 ? (edges[2 * i] + 1) % n : st_randomInt(0, n);
        }
        int *members, *starts, componentNumber;
        getComponents(testCase, n, edgeNumber, edges, &members, &starts, &componentNumber);
        int *component = st_malloc(n * sizeof(int));
        for (int i = 0; i < componentNumber; i++) {
            for (int j = starts[i]; j < starts[i + 1]; j++) {
                component[members[j]] = i;
            }
        }
        bool *same = getSameComponentByBruteForce(n, edgeNumber, edges);
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                CuAssertIntEquals(testCase, same[u * n + v], component[u] == component[v]);
            }
        }
        free(same);
        free(component);
        free(members);
        free(starts);
        free(edges);
    }
}

CuSuite* stThreeEdgeConnectedTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStThreeEdgeConnected_examples);
    SUITE_ADD_TEST(suite, testStThreeEdgeConnected_order);
    SUITE_ADD_TEST(suite, testStThreeEdgeConnected_random);
    return suite;
}

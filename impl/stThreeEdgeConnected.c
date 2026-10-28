/*
 * stThreeEdgeConnected.c
 *
 * The 3-edge-connected components of an undirected multigraph, by the absorb-eject algorithm of
 * Norouzi and Tsin, "A simple 3-edge connected component algorithm revisited", Information
 * Processing Letters 114 (2014), which revises Tsin, "A simple 3-edge-connected component
 * algorithm", Theory of Computing Systems 40 (2007).
 *
 * This is a C translation of the independent implementation of that algorithm in vg
 * (https://github.com/vgteam/vg, src/algorithms/three_edge_connected_components.cpp), which is
 * released under the MIT license:
 *
 * Copyright (c) 2014 Erik Garrison
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "sonLib.h"
#include "stThreeEdgeConnected.h"

/*
 * The algorithm in brief. A depth first search walks the graph, and on coming back up each tree
 * edge it absorbs into the parent the vertices that must share its 3-edge-connected component,
 * having first ejected the child if its effective degree (its degree in the graph with every
 * absorb so far applied) is 2, as the child's component is then cut off by two edges and complete.
 * None of the graph surgery is done: it suffices to track each vertex's effective degree, its low
 * point, and a path from it down to the vertex whose back edge gave that low point. Paths share
 * tails, so a path is stored as each vertex's successor on it, with a flag to take a vertex that
 * was absorbed off its own path. Self loops and bridges, which the paper excludes, are hidden from
 * the degree counts, and nothing is absorbed across a bridge.
 *
 * As in the paper, each vertex that has not been absorbed keeps the list of the vertices absorbed
 * into it, and a component is output as that list when it is complete: when the vertex is ejected
 * or cut off by a bridge, or is the root of the search. The lists are circular, so absorbing one
 * into another is a swap of two successors. The edges of each vertex are looked at from the last
 * to the first, and with the circular lists that puts the components, and the vertices in each,
 * in the same order as the implementation this replaced did, which the cactus graph built from
 * them depends on.
 */

typedef struct _threeEdgeState {
    int *dfsCounter; //1 + the number of vertices visited before this one, or 0 if it is unvisited
    int *dfsExit; //the DFS counter after the vertex's subtree is done
    int *lowPoint; //the lowest DFS counter that the vertex's subtree has a back edge to
    int *effectiveDegree;
    int *pathTail; //the vertex's successor on its path, or -1 where the path ends
    char *isOnPath;
    int *nextMember; //the vertex's successor in the circular list of its component
    int *members, *starts, components, membersLength; //the output
} ThreeEdgeState;

/*
 * Absorbs into one component the vertices on the path from start that are on it, stopping before
 * pastEnd, or at the end of the path if pastEnd is -1. They are absorbed into into, or into the
 * first of them if into is -1, and into's effective degree is updated as if the edge joining each
 * one to it was contracted.
 */
static void absorbPath(ThreeEdgeState *s, int into, int start, int pastEnd) {
    bool merged = false; //only read by the assert
    (void) merged;
    for (int v = start; v != pastEnd; v = s->pathTail[v]) {
        if (v == -1) {
            //Only a path without edges can end before reaching pastEnd
            assert(!merged);
            break;
        }
        if (s->isOnPath[v]) {
            if (into == -1) {
                into = v;
            } else {
                merged = true;
                s->effectiveDegree[into] += s->effectiveDegree[v] - 2;
                int t = s->nextMember[into];
                s->nextMember[into] = s->nextMember[v];
                s->nextMember[v] = t;
            }
        }
    }
}

//Outputs the component of the vertices absorbed into v, v first
static void outputComponent(ThreeEdgeState *s, int v) {
    s->starts[s->components++] = s->membersLength;
    int u = v;
    do {
        s->members[s->membersLength++] = u;
        u = s->nextMember[u];
    } while (u != v);
}

typedef struct _threeEdgeFrame {
    int vertex;
    int edge; //index in adj of the edge being looked at, counting down
    bool sawParentTreeEdge; //the first edge back to the parent is the tree edge, any more are back edges
    bool recursing; //whether the search went down the edge, and has to finish it on return
} ThreeEdgeFrame;

static void visit(ThreeEdgeState *s, int v, int *dfsCount) {
    s->dfsCounter[v] = s->lowPoint[v] = (*dfsCount)++;
    s->pathTail[v] = -1;
    s->isOnPath[v] = 1;
    s->nextMember[v] = v;
}

void stThreeEdgeConnected_getComponents(int n, const int *offsets, const int *adj, int **membersOut, int **startsOut, int *nComponentsOut) {
    int m = n > 0 ? n : 1;
    ThreeEdgeState s;
    s.dfsCounter = st_calloc(m, sizeof(int));
    s.dfsExit = st_malloc(m * sizeof(int));
    s.lowPoint = st_malloc(m * sizeof(int));
    s.effectiveDegree = st_calloc(m, sizeof(int));
    s.pathTail = st_malloc(m * sizeof(int));
    s.isOnPath = st_malloc(m * sizeof(char));
    s.nextMember = st_malloc(m * sizeof(int));
    s.members = st_malloc(m * sizeof(int));
    s.starts = st_malloc((n + 1) * sizeof(int));
    s.components = 0;
    s.membersLength = 0;
    //An explicit stack, as the search can go as deep as there are vertices
    ThreeEdgeFrame *stack = st_malloc(m * sizeof(ThreeEdgeFrame));

    int dfsCount = 1;
    for (int root = 0; root < n; root++) {
        if (s.dfsCounter[root] != 0) {
            continue;
        }
        visit(&s, root, &dfsCount);
        stack[0] = (ThreeEdgeFrame) { root, offsets[root + 1] - 1, 0, 0 };
        int depth = 1;
        while (depth > 0) {
            ThreeEdgeFrame *frame = &stack[depth - 1];
            int v = frame->vertex;
            if (frame->edge < offsets[v]) {
                s.dfsExit[v] = dfsCount;
                depth--;
                continue;
            }
            int w = adj[frame->edge];
            if (!frame->recursing) {
                s.effectiveDegree[v]++;
                if (s.dfsCounter[w] == 0) {
                    frame->recursing = 1;
                    visit(&s, w, &dfsCount);
                    stack[depth++] = (ThreeEdgeFrame) { w, offsets[w + 1] - 1, 0, 0 };
                    continue;
                }
                if (depth > 1 && w == stack[depth - 2].vertex && !frame->sawParentTreeEdge) {
                    frame->sawParentTreeEdge = 1;
                } else if (s.dfsCounter[w] < s.dfsCounter[v]) {
                    //An outgoing back edge (step 1.2 of the paper). If it reaches lower than our low
                    //point, absorb our path and start a new one, from us, for the new low point.
                    if (s.dfsCounter[w] < s.lowPoint[v]) {
                        absorbPath(&s, -1, v, -1);
                        s.lowPoint[v] = s.dfsCounter[w];
                        s.isOnPath[v] = 1;
                        s.pathTail[v] = -1;
                    }
                } else if (s.dfsCounter[v] < s.dfsCounter[w]) {
                    //An incoming back edge (step 1.3), from a descendant w that closes a cycle
                    //through us. Absorb our path as far as w. w may have been absorbed off the path,
                    //in which case go as far as the last vertex on it that is an ancestor of w,
                    //which is what w was absorbed into.
                    s.effectiveDegree[v] -= 2;
                    int last = v;
                    int next = s.pathTail[last];
                    while (next != -1 && s.dfsCounter[next] <= s.dfsCounter[w] && s.dfsExit[next] >= s.dfsExit[w]) {
                        last = next;
                        next = s.pathTail[last];
                    }
                    absorbPath(&s, -1, v, s.pathTail[last]);
                    //Our path continues from where the absorbing stopped
                    s.pathTail[v] = s.pathTail[last];
                } else {
                    s.effectiveDegree[v]--; //a self loop
                }
                frame->edge--;
            } else {
                //Back from the search down the tree edge to w (step 1.1)
                if (s.lowPoint[w] == s.dfsCounter[w]) {
                    //Nothing under w has a back edge above it, so the tree edge is a bridge
                    s.effectiveDegree[w]--;
                    s.effectiveDegree[v]--;
                    outputComponent(&s, w);
                } else {
                    if (s.effectiveDegree[w] == 2) {
                        s.isOnPath[w] = 0; //eject w
                        outputComponent(&s, w);
                    }
                    assert(s.effectiveDegree[w] != 1); //as the bridges are hidden
                    if (s.lowPoint[v] <= s.lowPoint[w]) {
                        absorbPath(&s, v, w, -1);
                    } else {
                        s.lowPoint[v] = s.lowPoint[w];
                        absorbPath(&s, -1, v, -1);
                        s.isOnPath[v] = 1;
                        s.pathTail[v] = w;
                    }
                }
                frame->recursing = 0;
                frame->edge--;
            }
        }
        outputComponent(&s, root);
    }
    assert(s.membersLength == n);
    s.starts[s.components] = n;

    free(stack);
    free(s.dfsCounter);
    free(s.dfsExit);
    free(s.lowPoint);
    free(s.effectiveDegree);
    free(s.pathTail);
    free(s.isOnPath);
    free(s.nextMember);
    *membersOut = s.members;
    *startsOut = s.starts;
    *nComponentsOut = s.components;
}

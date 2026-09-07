/*
 * stCactusGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "stCactusGraphs.h"
#include "3_Absorb3edge2x.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

/*
 * The nodes of a graph are kept on an intrusive doubly linked list, in the order the old
 * object-to-node hash iterated them: creation order, with a merged node moving to the tail (the hash
 * removed and re-inserted it). That order decides how the graph is traversed and so what the caller
 * gets out, so it is preserved exactly. Nodes and edge ends come from chunk arenas owned by the graph:
 * a merged-away node is simply unlinked, and destruction frees the chunks.
 */
struct _stCactusNode {
    stCactusEdgeEnd *head;
    stCactusEdgeEnd *tail;
    void *nodeObject;
    stCactusNode *prevNode, *nextNode;
    stCactusEdgeEnd *onChainEdgeEnd; //scratch for stCactusGraph_markCycles
    int64_t index; //scratch: dense id for the 3-edge-connected computation
    uint64_t mark; //scratch: the epoch the node was last visited in by a component search
};

struct _stCactusEdgeEnd {
    stCactusEdgeEnd *otherEdgeEnd;
    stCactusNode *node;
    stCactusEdgeEnd *nEdgeEnd;
    void *endObject;
    stCactusEdgeEnd *link;
    bool linkOrientation;
    bool isChainEnd;
};

#define ST_CACTUS_CHUNK_SIZE 4096

struct _stCactusGraph {
    stCactusNode *firstNode, *lastNode;
    int64_t nodeNumber;
    uint64_t markEpoch;
    stList *nodeChunks; //arrays of ST_CACTUS_CHUNK_SIZE nodes
    int64_t nodeChunkUsed;
    stList *edgeEndChunks; //arrays of ST_CACTUS_CHUNK_SIZE edge ends
    int64_t edgeEndChunkUsed;
    void (*destructNodeObjectFn)(void *);
    void (*destructEdgeEndObjectFn)(void *);
};

//Node functions

static void stCactusGraph_appendNode(stCactusGraph *graph, stCactusNode *node) {
    node->prevNode = graph->lastNode;
    node->nextNode = NULL;
    if (graph->lastNode == NULL) {
        graph->firstNode = node;
    } else {
        graph->lastNode->nextNode = node;
    }
    graph->lastNode = node;
    graph->nodeNumber++;
}

static void stCactusGraph_unlinkNode(stCactusGraph *graph, stCactusNode *node) {
    if (node->prevNode == NULL) {
        assert(graph->firstNode == node);
        graph->firstNode = node->nextNode;
    } else {
        node->prevNode->nextNode = node->nextNode;
    }
    if (node->nextNode == NULL) {
        assert(graph->lastNode == node);
        graph->lastNode = node->prevNode;
    } else {
        node->nextNode->prevNode = node->prevNode;
    }
    node->prevNode = node->nextNode = NULL;
    graph->nodeNumber--;
}

stCactusNode *stCactusNode_construct(stCactusGraph *graph, void *nodeObject) {
    if (graph->nodeChunkUsed == ST_CACTUS_CHUNK_SIZE || stList_length(graph->nodeChunks) == 0) {
        stList_append(graph->nodeChunks, st_calloc(ST_CACTUS_CHUNK_SIZE, sizeof(stCactusNode)));
        graph->nodeChunkUsed = 0;
    }
    stCactusNode *node = &((stCactusNode *) stList_peek(graph->nodeChunks))[graph->nodeChunkUsed++];
    node->nodeObject = nodeObject;
    stCactusGraph_appendNode(graph, node);
    return node;
}

void *stCactusNode_getObject(stCactusNode *node) {
    return node->nodeObject;
}

stCactusNodeEdgeEndIt stCactusNode_getEdgeEndIt(stCactusNode *node) {
    stCactusNodeEdgeEndIt it;
    it.edgeEnd = node->head;
    return it;
}

stCactusEdgeEnd *stCactusNodeEdgeEndIt_getNext(stCactusNodeEdgeEndIt *it) {
    stCactusEdgeEnd *edgeEnd = it->edgeEnd;
    if (edgeEnd != NULL) {
        it->edgeEnd = it->edgeEnd->nEdgeEnd;
    }
    return edgeEnd;
}

stCactusEdgeEnd *stCactusNode_getFirstEdgeEnd(stCactusNode *node) {
    return node->head;
}

int64_t stCactusGraph_getNodeNumber(stCactusGraph *graph) {
    return graph->nodeNumber;
}

//Private node functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edge, void(*destructEdgeEndObjectFn)(void *));

static void stCactusNode_destruct(stCactusNode *node, void(*destructNodeObjectFn)(void *), void(*destructEdgeEndObjectFn)(void *)) {
    stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&it)) != NULL) {
        stCactusEdgeEnd_destruct(edgeEnd, destructEdgeEndObjectFn);
    }
    if (destructNodeObjectFn != NULL) {
        destructNodeObjectFn(node->nodeObject);
    }
    //the node's memory belongs to the graph's arena
}

static void stCactusNode_mergeNodes(stCactusGraph *graph, stCactusNode *node1, stCactusNode *node2,
        void *(*mergeNodeObjects)(void *, void *)) {
    if (node1 == node2) {
        return;
    }
    //The surviving node moves to the end of the node order, as the hash re-insert used to put it
    stCactusGraph_unlinkNode(graph, node2);
    stCactusGraph_unlinkNode(graph, node1);
    node1->nodeObject = mergeNodeObjects(node1->nodeObject, node2->nodeObject);
    stCactusGraph_appendNode(graph, node1);
    if (node2->head != NULL) {
        if (node1->head == NULL) {
            node1->head = node2->head;
        } else {
            assert(node1->tail->nEdgeEnd == NULL);
            node1->tail->nEdgeEnd = node2->head;
        }
        stCactusEdgeEnd *edgeEnd = node2->head;
        while (edgeEnd != NULL) {
            edgeEnd->node = node1;
            edgeEnd = edgeEnd->nEdgeEnd;
        }
        node1->tail = node2->tail;
    }
    //node2's memory belongs to the graph's arena
}

//Edge functions

static void connectUpEdgeEnd(stCactusEdgeEnd *edgeEnd, stCactusNode *node, stCactusEdgeEnd *otherEdgeEnd, void *endObject) {
    edgeEnd->node = node;
    edgeEnd->otherEdgeEnd = otherEdgeEnd;
    edgeEnd->endObject = endObject;
    if (node->head == NULL) {
        node->head = edgeEnd;
    } else {
        node->tail->nEdgeEnd = edgeEnd;
    }
    node->tail = edgeEnd;
}

static stCactusEdgeEnd *stCactusGraph_allocateEdgeEnd(stCactusGraph *graph) {
    if (graph->edgeEndChunkUsed == ST_CACTUS_CHUNK_SIZE || stList_length(graph->edgeEndChunks) == 0) {
        stList_append(graph->edgeEndChunks, st_calloc(ST_CACTUS_CHUNK_SIZE, sizeof(stCactusEdgeEnd)));
        graph->edgeEndChunkUsed = 0;
    }
    return &((stCactusEdgeEnd *) stList_peek(graph->edgeEndChunks))[graph->edgeEndChunkUsed++];
}

stCactusEdgeEnd *stCactusEdgeEnd_construct(stCactusGraph *graph, stCactusNode *node1, stCactusNode *node2, void *edgeEndObject1,
        void *edgeEndObject2) {
    stCactusEdgeEnd *edgeEnd1 = stCactusGraph_allocateEdgeEnd(graph);
    stCactusEdgeEnd *edgeEnd2 = stCactusGraph_allocateEdgeEnd(graph);

    connectUpEdgeEnd(edgeEnd1, node1, edgeEnd2, edgeEndObject1);
    connectUpEdgeEnd(edgeEnd2, node2, edgeEnd1, edgeEndObject2);
    return edgeEnd1;
}

void *stCactusEdgeEnd_getObject(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->endObject;
}

stCactusNode *stCactusEdgeEnd_getNode(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->node;
}

stCactusNode *stCactusEdgeEnd_getOtherNode(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->otherEdgeEnd->node;
}

stCactusEdgeEnd *stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->otherEdgeEnd;
}

stCactusEdgeEnd *stCactusEdgeEnd_getLink(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->link;
}

bool stCactusEdgeEnd_getLinkOrientation(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->linkOrientation;
}

bool stCactusEdgeEnd_isChainEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->isChainEnd;
}

stCactusEdgeEnd *stCactusEdgeEnd_getNextEdgeEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->nEdgeEnd;
}

//Private edge functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edgeEnd, void(*destructEdgeEndObjectFn)(void *)) {
    if (destructEdgeEndObjectFn != NULL) {
        destructEdgeEndObjectFn(edgeEnd->endObject);
    }
    //the edge end's memory belongs to the graph's arena
}

static void stCactusEdgeEnd_setLink(stCactusEdgeEnd *edgeEnd, stCactusEdgeEnd *otherEdgeEnd) {
    edgeEnd->link = otherEdgeEnd;
}

static void stCactusEdgeEnd_setLinkOrientation(stCactusEdgeEnd *edgeEnd, bool orientation) {
    edgeEnd->linkOrientation = orientation;
}

static void stCactusEdgeEnd_setIsChainEnd(stCactusEdgeEnd *edgeEnd, bool isChainEnd) {
    edgeEnd->isChainEnd = isChainEnd;
}

//Graph functions

stCactusGraph *stCactusGraph_construct2(void(*destructNodeObjectFn)(void *), void(*destructEdgeEndObjectFn)(void *)) {
    stCactusGraph *cactusGraph = st_calloc(1, sizeof(stCactusGraph));
    cactusGraph->nodeChunks = stList_construct3(0, free);
    cactusGraph->edgeEndChunks = stList_construct3(0, free);
    cactusGraph->destructNodeObjectFn = destructNodeObjectFn;
    cactusGraph->destructEdgeEndObjectFn = destructEdgeEndObjectFn;
    return cactusGraph;
}

stCactusGraph *stCactusGraph_construct() {
    return stCactusGraph_construct2(NULL, NULL);
}

void stCactusGraph_destruct(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(graph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNode_destruct(node, graph->destructNodeObjectFn, graph->destructEdgeEndObjectFn);
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    stList_destruct(graph->nodeChunks);
    stList_destruct(graph->edgeEndChunks);
    free(graph);
}

stCactusNode *stCactusGraph_getNode(stCactusGraph *graph, void *nodeObject) {
    //A linear scan: only the tests look nodes up by object
    for (stCactusNode *node = graph->firstNode; node != NULL; node = node->nextNode) {
        if (node->nodeObject == nodeObject) {
            return node;
        }
    }
    return NULL;
}

/*
 * A stack of pointers on a plain array. The depth first searches below push and pop at every edge
 * end, and the bounds checks and calls a list costs per operation were a measurable share of every
 * cactus graph build. Pushes and pops are in the same order as with the list they replace.
 */
typedef struct _ptrStack {
    void **items;
    int64_t length, capacity;
} PtrStack;

static void ptrStack_init(PtrStack *stack) {
    stack->capacity = 1024;
    stack->length = 0;
    stack->items = st_malloc(stack->capacity * sizeof(void *));
}

static inline void ptrStack_push(PtrStack *stack, void *item) {
    if (stack->length == stack->capacity) {
        stack->capacity *= 2;
        stack->items = st_realloc(stack->items, stack->capacity * sizeof(void *));
    }
    stack->items[stack->length++] = item;
}

static inline void *ptrStack_pop(PtrStack *stack) {
    assert(stack->length > 0);
    return stack->items[--stack->length];
}

static inline void *ptrStack_peek(PtrStack *stack) {
    assert(stack->length > 0);
    return stack->items[stack->length - 1];
}

static void ptrStack_destruct(PtrStack *stack) {
    free(stack->items);
}

void stCactusGraph_collapseToCactus(stCactusGraph *graph, void *(*mergeNodeObjects)(void *, void *), stCactusNode *startNode) {
    //computeThreeEdgeConnectedComponents indexes vertices with int, so bail out here rather
    //than truncate the node count into an int and work on a graph that is not the input
    int64_t nodeNumber = stCactusGraph_getNodeNumber(graph);
    if (nodeNumber >= INT32_MAX) {
        st_errAbort("Cactus graph has %" PRIi64 " nodes, exceeding the %" PRIi64 " node limit of the "
                    "3-edge-connected component algorithm", nodeNumber, (int64_t)INT32_MAX - 1);
    }

    //Number the nodes in iteration order and lay the graph out in compressed sparse row form,
    //edge ends in their list order, for the three edge connected code
    stCactusNode **nodes = st_malloc((nodeNumber > 0 ? nodeNumber : 1) * sizeof(stCactusNode *));
    int *offsets = st_malloc((nodeNumber + 1) * sizeof(int));
    int64_t edgeEndNumber = 0;
    int64_t nodeIdCounter = 0;
    for (stCactusNode *node = graph->firstNode; node != NULL; node = node->nextNode) {
        node->index = nodeIdCounter;
        nodes[nodeIdCounter] = node;
        offsets[nodeIdCounter++] = edgeEndNumber;
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        while (stCactusNodeEdgeEndIt_getNext(&edgeEndIt) != NULL) {
            edgeEndNumber++;
        }
    }
    offsets[nodeNumber] = edgeEndNumber;
    if (edgeEndNumber >= INT32_MAX) {
        st_errAbort("Cactus graph has %" PRIi64 " edge ends, exceeding the limit of the 3-edge-connected component algorithm", edgeEndNumber);
    }
    int *adj = st_malloc((edgeEndNumber > 0 ? edgeEndNumber : 1) * sizeof(int));
    int64_t k = 0;
    for (stCactusNode *node = graph->firstNode; node != NULL; node = node->nextNode) {
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            adj[k++] = stCactusEdgeEnd_getOtherNode(edgeEnd)->index;
        }
    }
    assert(k == edgeEndNumber);
    int *members, *starts, componentNumber;
    computeThreeEdgeConnectedComponentsCSR(nodeNumber, offsets, adj, &members, &starts, &componentNumber);
    //Now do the merging
    for (int64_t i = 0; i < componentNumber; i++) {
        stCactusNode *node = nodes[members[starts[i]]];
        for (int64_t j = starts[i] + 1; j < starts[i + 1]; j++) {
            stCactusNode *otherNode = nodes[members[j]];
            assert(node != otherNode);
            if (otherNode == startNode) { //This prevents the start node from being destructed.
                otherNode = node;
                node = startNode;
            }
            stCactusNode_mergeNodes(graph, node, otherNode, mergeNodeObjects);
        }
    }
    //Cleanup
    free(members);
    free(starts);
    free(nodes);
    free(offsets);
    free(adj);

    //Mark the cycles, one connected component at a time, starting each from its first node in
    //iteration order (or from startNode, for the component that holds it)
    PtrStack stack;
    ptrStack_init(&stack);
    uint64_t firstEpoch = graph->markEpoch + 1;
    for (stCactusNode *root = graph->firstNode; root != NULL; root = root->nextNode) {
        if (root->mark >= firstEpoch) {
            continue; //already in a component found in this pass
        }
        uint64_t epoch = ++graph->markEpoch;
        root->mark = epoch;
        ptrStack_push(&stack, root);
        while (stack.length > 0) {
            stCactusNode *node = ptrStack_pop(&stack);
            stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd;
            while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
                stCactusNode *otherNode = stCactusEdgeEnd_getOtherNode(edgeEnd);
                if (otherNode->mark != epoch) {
                    otherNode->mark = epoch;
                    ptrStack_push(&stack, otherNode);
                }
            }
        }
        if (startNode != NULL && startNode->mark == epoch) {
            stCactusGraph_markCycles(graph, startNode);
        } else {
            stCactusGraph_markCycles(graph, root);
        }
    }
    ptrStack_destruct(&stack);
}

stCactusGraphNodeIt *stCactusGraphNodeIterator_construct(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = st_malloc(sizeof(stCactusGraphNodeIt));
    nodeIt->next = graph->firstNode;
    nodeIt->graph = graph;
    return nodeIt;
}

stCactusNode *stCactusGraphNodeIterator_getNext(stCactusGraphNodeIt *nodeIt) {
    stCactusNode *node = nodeIt->next;
    if (node != NULL) {
        nodeIt->next = node->nextNode; //read now, so the caller may destroy the returned node
    }
    return node;
}

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIt *nodeIt) {
    free(nodeIt);
}

void stCactusGraph_unmarkCycles(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(graph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stCactusEdgeEnd_setIsChainEnd(edgeEnd, 0);
            stCactusEdgeEnd_setLink(edgeEnd, NULL);
            stCactusEdgeEnd_setLinkOrientation(edgeEnd, 0);
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
}

static void makeChain(stCactusEdgeEnd *edgeEnd, stCactusEdgeEnd *edgeEnd2, PtrStack *chainPath) {
    assert(edgeEnd != edgeEnd2);
    stCactusEdgeEnd_setLink(edgeEnd, edgeEnd2);
    stCactusEdgeEnd_setLink(edgeEnd2, edgeEnd);
    stCactusEdgeEnd_setIsChainEnd(edgeEnd, 1);
    stCactusEdgeEnd_setIsChainEnd(edgeEnd2, 1);
    stCactusEdgeEnd_setLinkOrientation(edgeEnd2, 1);
    for (int64_t j = chainPath->length - 1; j >= 0; j -= 2) {
        assert(j >= 1);
        stCactusEdgeEnd *edgeEnd3 = chainPath->items[j];
        stCactusEdgeEnd *edgeEnd4 = chainPath->items[j - 1];
        assert(edgeEnd3 != edgeEnd2);
        assert(edgeEnd4 != edgeEnd2);
        assert(edgeEnd4 != edgeEnd);
        if (edgeEnd3 == edgeEnd) {
            break;
        }
        stCactusEdgeEnd_setLink(edgeEnd3, edgeEnd4);
        stCactusEdgeEnd_setLink(edgeEnd4, edgeEnd3);
        stCactusEdgeEnd_setLinkOrientation(edgeEnd4, 1);
    }
}

/*
 * The nodes on the current path each remember the edge end the path left them by, in onChainEdgeEnd;
 * it is NULL for every node off the path, and every node is off the path again on return.
 */
static void stCactusGraph_markCyclesP(stCactusEdgeEnd *edgeEnd, PtrStack *chainPath) {
    ptrStack_push(chainPath, edgeEnd);
    ptrStack_push(chainPath, NULL);
    while (chainPath->length > 0) {
        stCactusEdgeEnd *edgeEnd2 = ptrStack_pop(chainPath);
        edgeEnd = ptrStack_peek(chainPath);
        stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
        if (edgeEnd2 == NULL) {
            edgeEnd2 = stCactusNode_getFirstEdgeEnd(node);
        } else {
            edgeEnd2 = stCactusEdgeEnd_getNextEdgeEnd(edgeEnd2);
        }
        if (edgeEnd2 == NULL) {
            ptrStack_pop(chainPath);
            node->onChainEdgeEnd = NULL;
        } else {
            ptrStack_push(chainPath, edgeEnd2);
            if (edgeEnd2 != edgeEnd && stCactusEdgeEnd_getLink(edgeEnd2) == NULL) {
                node->onChainEdgeEnd = edgeEnd2;
                stCactusEdgeEnd *edgeEnd3 = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2);
                stCactusEdgeEnd *edgeEnd4;
                if ((edgeEnd4 = stCactusEdgeEnd_getNode(edgeEnd3)->onChainEdgeEnd) != NULL) { //We've traversed a cycle
                    makeChain(edgeEnd4, edgeEnd3, chainPath);
                } else {
                    ptrStack_push(chainPath, edgeEnd3);
                    ptrStack_push(chainPath, NULL);
                }
            }
        }
    }
}

void stCactusGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode) {
    PtrStack chainPath;
    ptrStack_init(&chainPath);
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd);
            if (stCactusEdgeEnd_getNode(edgeEnd2) == startNode) {
                makeChain(edgeEnd, edgeEnd2, &chainPath);
            } else {
                startNode->onChainEdgeEnd = edgeEnd;
                stCactusGraph_markCyclesP(edgeEnd2, &chainPath);
                startNode->onChainEdgeEnd = NULL;
            }
        }
    }
    ptrStack_destruct(&chainPath);
#ifndef NDEBUG
    for (stCactusNode *node = graph->firstNode; node != NULL; node = node->nextNode) {
        assert(node->onChainEdgeEnd == NULL);
    }
#endif
}

static void stCactusGraph_collapseBridgesP(stCactusNode *parentNode, stCactusEdgeEnd *edgeEnd, stList *nodesToMerge) {
    PtrStack stack;
    ptrStack_init(&stack);
    ptrStack_push(&stack, parentNode);
    ptrStack_push(&stack, edgeEnd);
    while (stack.length > 0) {
        edgeEnd = ptrStack_pop(&stack);
        parentNode = ptrStack_pop(&stack);
        stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
        if (stCactusEdgeEnd_getLink(edgeEnd) == NULL) { //Is a bridge
            //Establish if this is a leaf
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd2;
            int64_t bridges = 0;
            while ((edgeEnd2 = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                if (edgeEnd2 != edgeEnd) {
                    if (stCactusEdgeEnd_getLink(edgeEnd2) == NULL) {
                        bridges++;
                    }
                    if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                        ptrStack_push(&stack, parentNode);
                        ptrStack_push(&stack, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                    }
                }
            }
            if (bridges != 1) {
                stList_append(nodesToMerge, parentNode);
                stList_append(nodesToMerge, node);
            }
        } else if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd2;
            while ((edgeEnd2 = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                if (edgeEnd2 != edgeEnd && !stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                    ptrStack_push(&stack, node);
                    ptrStack_push(&stack, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                }
            }
        }
    }
    ptrStack_destruct(&stack);
}

void mergeNodes(stCactusGraph *graph, stCactusNode *startNode, stList *nodesToMerge, void *(*mergeNodeObjects)(void *, void *)) {
    for (int64_t i = 0; i < stList_length(nodesToMerge); i += 2) {
        stCactusNode *parentNode = stList_get(nodesToMerge, i);
        stCactusNode *nodeToMerge = stList_get(nodesToMerge, i + 1);
        stCactusNode_mergeNodes(graph, parentNode, nodeToMerge, mergeNodeObjects);
    }
    stList_destruct(nodesToMerge);
    stCactusGraph_unmarkCycles(graph);
    stCactusGraph_markCycles(graph, startNode);
}

void stCactusGraph_collapseBridges(stCactusGraph *graph, stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *)) {
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    stList *bridgesToMerge = stList_construct();
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stCactusGraph_collapseBridgesP(startNode, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd), bridgesToMerge);
        }
    }
    mergeNodes(graph, startNode, bridgesToMerge, mergeNodeObjects); //Now merge bridges
}

int64_t stCactusEdgeEnd_getChainLength(stCactusEdgeEnd *edgeEnd) {
    stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd));
    int64_t chainLength = 1;
    while (edgeEnd != edgeEnd2) {
        chainLength++;
        edgeEnd2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
    }
    return chainLength;
}

int64_t stCactusNode_getTotalEdgeLengthOfFlower(stCactusNode *node) {
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    int64_t totalLength = 0;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (stCactusEdgeEnd_isChainEnd(edgeEnd) && stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            totalLength += stCactusEdgeEnd_getChainLength(edgeEnd);
        }
    }
    return totalLength;
}

int64_t stCactusNode_getChainNumber(stCactusNode *node) {
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    int64_t totalNumber = 0;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (stCactusEdgeEnd_isChainEnd(edgeEnd) && stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            totalNumber++;
        }
    }
    return totalNumber;
}

static void stCactusGraph_collapseLongChainsP(stCactusNode *node, stList *nodesToMerge, int64_t longChain) {
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (stCactusEdgeEnd_isChainEnd(edgeEnd) && stCactusEdgeEnd_getLinkOrientation(edgeEnd) && stCactusEdgeEnd_getChainLength(edgeEnd)
                > longChain) {
            stCactusNode *node2 = stCactusEdgeEnd_getNode(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd));
            stCactusNode *node3 = stCactusEdgeEnd_getNode(stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd)));
            if (node2 != node3) {
                assert(node2 != node);
                assert(node3 != node);
                stList_append(nodesToMerge, node3);
                stList_append(nodesToMerge, node2); //subtle, but this ordering ensures that when
                //we remark the nodes with orientations they end up with different orientations between the edges that are
                //really linked but now merged in the same node.
            } else {
                assert(longChain <= 1);
            }
        }
    }
}

stSet *stCactusGraph_collapseLongChainsOfBigFlowers(stCactusGraph *graph, stCactusNode *startNode, int64_t chainLengthForBigFlower,
        int64_t longChain, void *(*mergeNodeObjects)(void *, void *), bool recursive) {
    int64_t nodesMerged = 0;
    stSet *bigNodes = stSet_construct();
    do {
        stCactusGraphNodeIt *it = stCactusGraphNodeIterator_construct(graph);
        stCactusNode *node;
        stList *chainNodesToMerge = stList_construct();
        while ((node = stCactusGraphNodeIterator_getNext(it))) {
            if (stCactusNode_getTotalEdgeLengthOfFlower(node) > chainLengthForBigFlower) {
                stCactusGraph_collapseLongChainsP(node, chainNodesToMerge, longChain);
                stSet_insert(bigNodes, node);
            }
        }
        stCactusGraphNodeIterator_destruct(it);
        if(stList_length(chainNodesToMerge) > 0) { //We go around because the result of the merge may require more merges.
            nodesMerged += stList_length(chainNodesToMerge)/2;
            mergeNodes(graph, startNode, chainNodesToMerge, mergeNodeObjects); //Now merge bridges
        }
        else {
            stList_destruct(chainNodesToMerge);
            break;
        }
    } while(recursive);
    return bigNodes;
}

/*
 * Following functions are used for breaking chains by edge-ends that can-not be part of a link, i.e. because they have an incident effective self loops.
 */

static void stCactusEdgeEnd_mergeIncidentNodesInMarkedCactusGraph(stCactusGraph *graph,
        stCactusEdgeEnd *leftEdgeEnd, stCactusEdgeEnd *rightEdgeEnd,
        void *(*mergeNodeObjects)(void *, void *)) {
    /*
     * Merges together the nodes incident with the given edge ends, if distinct, updating the
     * chain links in the process.
     */
    if(stCactusEdgeEnd_getNode(leftEdgeEnd) != stCactusEdgeEnd_getNode(rightEdgeEnd)) {
        stCactusNode_mergeNodes(graph, stCactusEdgeEnd_getNode(leftEdgeEnd), stCactusEdgeEnd_getNode(rightEdgeEnd), mergeNodeObjects);
        assert(leftEdgeEnd->link != rightEdgeEnd);
        assert(rightEdgeEnd->link != leftEdgeEnd);
        assert(rightEdgeEnd->link != leftEdgeEnd->link);
        assert(stCactusEdgeEnd_getNode(leftEdgeEnd) == stCactusEdgeEnd_getNode(rightEdgeEnd));
        assert(stCactusEdgeEnd_getNode(leftEdgeEnd->link) == stCactusEdgeEnd_getNode(rightEdgeEnd->link));
        assert(stCactusEdgeEnd_getNode(leftEdgeEnd) == stCactusEdgeEnd_getNode(rightEdgeEnd->link));

        //Deal with the edge ends connected to their links
        leftEdgeEnd->link->link = rightEdgeEnd->link;
        rightEdgeEnd->link->link = leftEdgeEnd->link;
        leftEdgeEnd->link->isChainEnd = leftEdgeEnd->link->isChainEnd || rightEdgeEnd->link->isChainEnd;
        rightEdgeEnd->link->isChainEnd = leftEdgeEnd->link->isChainEnd;

        assert(leftEdgeEnd->link->link == rightEdgeEnd->link);
        assert(rightEdgeEnd->link->link == leftEdgeEnd->link);
        assert(leftEdgeEnd->link->linkOrientation != rightEdgeEnd->link->linkOrientation);
        assert(stCactusEdgeEnd_getNode(leftEdgeEnd->link) == stCactusEdgeEnd_getNode(rightEdgeEnd->link));

        //Deal with making the chain circular
        leftEdgeEnd->link = rightEdgeEnd;
        rightEdgeEnd->link = leftEdgeEnd;
        leftEdgeEnd->isChainEnd = 1;
        rightEdgeEnd->isChainEnd = 1;

        assert(leftEdgeEnd->link == rightEdgeEnd);
        assert(rightEdgeEnd->link == leftEdgeEnd);
        assert(leftEdgeEnd->linkOrientation != rightEdgeEnd->linkOrientation);
    }
}

static bool getNestedEdgeEndsNotInChain(stCactusEdgeEnd *edgeEnd,
        stCactusEdgeEnd **leftEdgeEnd, stCactusEdgeEnd **rightEdgeEnd,
        bool (*endIsNotInChain)(stCactusEdgeEnd *, void *extraArg), void *extraArg) {
    /*
     * Finds a pair of edge-ends for which (1) each end either returns non-zero for endIsNotInChain or is an edge-end at the end of chain
     * and (2) the ends are in the chain that linked by a chain of cactus edges and links, starting on both ends with cactus edges.
     * If no such pair exists, except for the very ends of the chain, then returns 0, else returns 1.
     */
    assert(stCactusEdgeEnd_getLinkOrientation(edgeEnd));
    *leftEdgeEnd = edgeEnd;
    while(1) {
        *rightEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd);
        if(stCactusEdgeEnd_isChainEnd(*rightEdgeEnd)) { //looped to beginning of chain
            return stCactusEdgeEnd_getNode(*leftEdgeEnd) != stCactusEdgeEnd_getNode(*rightEdgeEnd) &&
                    (stCactusEdgeEnd_isChainEnd(*leftEdgeEnd) || endIsNotInChain(*leftEdgeEnd, extraArg)); //We have distinct node to merge and the left end is a chain.
        }
        if(endIsNotInChain(*rightEdgeEnd, extraArg)) {
            return 1;
        }
        edgeEnd = stCactusEdgeEnd_getLink(*rightEdgeEnd);
        if(endIsNotInChain(edgeEnd, extraArg)) {
            *leftEdgeEnd = edgeEnd;
        }
    }
    return 0;
}

static stCactusNode *stCactusGraph_breakChainByEndsNotInChain(stCactusGraph *graph,
        stCactusEdgeEnd *edgeEnd, void *(*mergeNodeObjects)(void *, void *),
        bool (*endIsNotInChain)(stCactusEdgeEnd *, void *extraArg), void *extraArg) {
    /*
     * The edgeEnd argument represents the end of a chain. Method breaks up the chain so that the resulting set of chains
     * contains no internal edge end for which endIsNotInChain is true. This is deterministic and results in the resulting set of chains being maximal.
     */
    assert(stCactusEdgeEnd_isChainEnd(edgeEnd));
    assert(stCactusEdgeEnd_getLink(edgeEnd) != NULL);
    assert(stCactusEdgeEnd_getLinkOrientation(edgeEnd));

    //get location of first pair of nested parentheses
    stCactusEdgeEnd *leftEdgeEnd, *rightEdgeEnd;
    bool hasPair = getNestedEdgeEndsNotInChain(edgeEnd, &leftEdgeEnd, &rightEdgeEnd, endIsNotInChain, extraArg);

    //while pair of parentheses exists
    while(hasPair) {
        assert(stCactusEdgeEnd_getNode(leftEdgeEnd) != stCactusEdgeEnd_getNode(rightEdgeEnd));
        assert(stCactusEdgeEnd_getLinkOrientation(leftEdgeEnd));
        assert(!stCactusEdgeEnd_getLinkOrientation(rightEdgeEnd));
        //if we are going to brake the chain into two by merging an internal node with "node", then we add
        if(stCactusEdgeEnd_getNode(leftEdgeEnd) == stCactusEdgeEnd_getNode(edgeEnd)) {
            assert(stCactusEdgeEnd_isChainEnd(leftEdgeEnd));
            assert(!stCactusEdgeEnd_isChainEnd(rightEdgeEnd));
            stCactusEdgeEnd *nextChainEdgeEnd = stCactusEdgeEnd_getLink(rightEdgeEnd);
            stCactusEdgeEnd_mergeIncidentNodesInMarkedCactusGraph(graph, leftEdgeEnd, rightEdgeEnd, mergeNodeObjects);
            stCactusGraph_breakChainByEndsNotInChain(graph, nextChainEdgeEnd, mergeNodeObjects, endIsNotInChain, extraArg);
            hasPair = getNestedEdgeEndsNotInChain(edgeEnd, &leftEdgeEnd, &rightEdgeEnd, endIsNotInChain, extraArg);
        }
        else if(stCactusEdgeEnd_getNode(rightEdgeEnd) == stCactusEdgeEnd_getNode(edgeEnd)) {
            assert(!stCactusEdgeEnd_isChainEnd(leftEdgeEnd));
            assert(stCactusEdgeEnd_isChainEnd(rightEdgeEnd));
            stCactusEdgeEnd_mergeIncidentNodesInMarkedCactusGraph(graph, leftEdgeEnd, rightEdgeEnd, mergeNodeObjects);
            stCactusGraph_breakChainByEndsNotInChain(graph, leftEdgeEnd, mergeNodeObjects, endIsNotInChain, extraArg);
            hasPair = getNestedEdgeEndsNotInChain(edgeEnd, &leftEdgeEnd, &rightEdgeEnd, endIsNotInChain, extraArg);
        }
        else { //Two internal nodes
            assert(!stCactusEdgeEnd_isChainEnd(leftEdgeEnd));
            assert(!stCactusEdgeEnd_isChainEnd(rightEdgeEnd));
            stCactusEdgeEnd *nextChainEdgeEnd = stCactusEdgeEnd_getLink(rightEdgeEnd);
            stCactusEdgeEnd_mergeIncidentNodesInMarkedCactusGraph(graph, leftEdgeEnd, rightEdgeEnd, mergeNodeObjects);
            hasPair = getNestedEdgeEndsNotInChain(nextChainEdgeEnd, &leftEdgeEnd, &rightEdgeEnd, endIsNotInChain, extraArg);
            if(!hasPair) {
                hasPair = getNestedEdgeEndsNotInChain(edgeEnd, &leftEdgeEnd, &rightEdgeEnd, endIsNotInChain, extraArg);
            }
        }
        assert(stCactusEdgeEnd_isChainEnd(edgeEnd));
    }
    return stCactusEdgeEnd_getNode(edgeEnd);
}

stList *stCactusNode_getChains(stCactusNode *startNode) {
    /*
     * Collate a static list of chains attached to the node.
     */
    stList *chainStack = stList_construct();
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (stCactusEdgeEnd_isChainEnd(edgeEnd) && stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stList_append(chainStack, edgeEnd);
        }
    }
    return chainStack;
}

stCactusNode *stCactusGraph_breakChainsByEndsNotInChains(stCactusGraph *graph,
        stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *),
        bool (*endIsNotInChain)(stCactusEdgeEnd *, void *extraArg), void *extraArg) {
    /*
     * Removes all ends that return non-zero from endIsNotInChain from within internal links within chains, breaking them down deterministically and correctly, in that
     * all resulting chains are maximal given this constraint. Returns the startNode, which may change due to merges.
     */
    //First call function recursively
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        assert(stCactusEdgeEnd_getNode(edgeEnd) == startNode);
        if (stCactusEdgeEnd_isChainEnd(edgeEnd) && stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd);
            while(stCactusEdgeEnd_getNode(edgeEnd2) != startNode) {
                stCactusGraph_breakChainsByEndsNotInChains(graph, stCactusEdgeEnd_getNode(edgeEnd2), mergeNodeObjects, endIsNotInChain, extraArg);
                edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd2));
            }
            assert(stCactusEdgeEnd_getLink(edgeEnd2) == edgeEnd);
            assert(stCactusEdgeEnd_isChainEnd(edgeEnd2));
            assert(!stCactusEdgeEnd_getLinkOrientation(edgeEnd2));
            assert(stCactusEdgeEnd_isChainEnd(edgeEnd));
            assert(stCactusEdgeEnd_getLinkOrientation(edgeEnd));
        }
    }
    //Now process the chains we have for this node
    stList *chainStack = stCactusNode_getChains(startNode);
    while(stList_length(chainStack) > 0) {
        startNode = stCactusGraph_breakChainByEndsNotInChain(graph, stList_pop(chainStack), mergeNodeObjects, endIsNotInChain, extraArg);
    }
    stList_destruct(chainStack);
    return startNode;
}

static void buildComponent(stSet *component, stCactusNode *cactusNode, bool ignoreBridgeEdges) {
    PtrStack stack; // Nodes to be visited
    ptrStack_init(&stack);
    ptrStack_push(&stack, cactusNode);

    // While there exists nodes to visit
    while(stack.length > 0) {
        stCactusNode *cactusNode2 = ptrStack_pop(&stack);

        // If not already in the component
        if(stSet_search(component, cactusNode2) == NULL) {

            // Add to the component
            stSet_insert(component, cactusNode2);

            // Recursively add connected nodes
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(cactusNode2);
            stCactusEdgeEnd *edgeEnd;
            while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                assert(cactusNode2 == stCactusEdgeEnd_getNode(edgeEnd));

                // If not ignoreBridgeEdges or is not a bridge edge
                if(!ignoreBridgeEdges || stCactusEdgeEnd_getLink(edgeEnd) != NULL) {
                    // Add to nodes to visit
                    ptrStack_push(&stack, stCactusEdgeEnd_getNode(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd)));
                }
            }
        }
    }

    ptrStack_destruct(&stack);
}

stList *stCactusGraph_getComponents(stCactusGraph *cactusGraph,
                                    bool ignoreBridgeEdges) {
    stSet *seen = stSet_construct(); // Set of nodes already visited
    stList *components = stList_construct3(0, (void (*)(void *))stSet_destruct); // Components

    // For each node in the cactus graph
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *cactusNode;
    while((cactusNode = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {

        // Check if the node has not already been seen, otherwise it is already in
        // a component
        if(stSet_search(seen, cactusNode) == NULL) {

            // Make a new component
            stSet *component = stSet_construct();
            stList_append(components, component);

            // Recursively build the component
            buildComponent(component, cactusNode, ignoreBridgeEdges);

            // Add set of nodes to seen
            stSet_insertAll(seen, component);

            assert(stSet_search(seen, cactusNode) == cactusNode);
        }
    }

    // Cleanup
    stCactusGraphNodeIterator_destruct(nodeIt);
    stSet_destruct(seen);

    return components;
}

stHash *stCactusGraphComponents_getNodesToComponentsMap(stList *components) {
    /*
     * For set of components, as computed by stCactusGraph_getComponents returns
     * map (hash) of cactus nodes to their respective component (each a stSet).
     */
    stHash *nodesToComponents = stHash_construct();
    for(int64_t i=0; i<stList_length(components); i++) {
        stSet *component = stList_get(components, i);
        stSetIterator *componentIt = stSet_getIterator(component);
        stCactusNode *cactusNode;
        while((cactusNode = stSet_getNext(componentIt)) != NULL) {
            stHash_insert(nodesToComponents, cactusNode, component);
        }
        stSet_destructIterator(componentIt);
    }
    return nodesToComponents;
}

/*
 * Bridge graphs
 */

stBridgeNode *stBridgeNode_construct() {
    stBridgeNode *bridgeNode = st_malloc(sizeof(stBridgeNode));
    bridgeNode->connectedNodes = stList_construct();
    bridgeNode->bridgeEnds = stSet_construct();
    bridgeNode->cactusNodes = stSet_construct();
    bridgeNode->connectedCactusNodes = stSet_construct(); // Cactus nodes on a path between pair of bridges in cactus tree

    return bridgeNode;
}

void stBridgeNode_destruct(stBridgeNode *bridgeNode) {
    stList_destruct(bridgeNode->connectedNodes);
    stSet_destruct(bridgeNode->bridgeEnds);
    stSet_destruct(bridgeNode->cactusNodes);
    stSet_destruct(bridgeNode->connectedCactusNodes);
    free(bridgeNode);
}

static void buildBridgeComponent(stBridgeNode *bridgeNode, stCactusEdgeEnd *edgeEnd) {
    // Stack object containing pairs of cactus edge ends
    stList *stack = stList_construct();

    // Initialise
    stList_append(stack, edgeEnd);
    stList_append(stack, edgeEnd);

    // Map of edge ends to booleans indicating if the node incident with the edge end
    // is either incident with a bridge edge or on a path between two bridge edge ends
    stHash *edgeEndsToConnectedStatus = stHash_construct2(NULL, free);

    while(stList_length(stack) > 0) {
        // Get the edge end and its predecessor edge end
        edgeEnd = stList_pop(stack);
        stCactusEdgeEnd *pEdgeEnd = stList_pop(stack);

        // Get the incident node
        stCactusNode *cactusNode = stCactusEdgeEnd_getNode(edgeEnd);

        // Boolean flag to indicate if incident with a bridge edge end or on a path between
        // two bridge edge ends
        bool *onConnectedPath;

        // First time through this will be true
        if((onConnectedPath = stHash_search(edgeEndsToConnectedStatus, edgeEnd)) == NULL) {

            // Initialise onConnectedPath
            onConnectedPath = st_malloc(sizeof(bool));
            stCactusEdgeEnd *linkEdgeEnd = stCactusEdgeEnd_getLink(edgeEnd);
            assert(linkEdgeEnd == NULL || stCactusEdgeEnd_getNode(linkEdgeEnd) == cactusNode);
            // If has not link then is an edge end of a bridge and so node is on a connected path
            onConnectedPath[0] = linkEdgeEnd == NULL;

            // Add back to the stack to process after all its successors
            stList_append(stack, pEdgeEnd);
            stList_append(stack, edgeEnd);
            stHash_insert(edgeEndsToConnectedStatus, edgeEnd, onConnectedPath);

            // For each edge incident with cactusNode
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(cactusNode);
            stCactusEdgeEnd *edgeEnd2;
            while ((edgeEnd2 = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                assert(cactusNode == stCactusEdgeEnd_getNode(edgeEnd2));

                // If not the incoming chain
                if(edgeEnd2 != edgeEnd && edgeEnd2 != linkEdgeEnd) {

                    // If part of a chain
                    if(stCactusEdgeEnd_getLink(edgeEnd2) != NULL) {

                        // If positive orientation (to avoid going around twice)
                        if(stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {

                            // For each other node in the chain
                            stCactusEdgeEnd *edgeEnd3 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                            while (edgeEnd3 != edgeEnd2) {
                                assert(stCactusEdgeEnd_getNode(edgeEnd3) != cactusNode);

                                // Add to the stack
                                stList_append(stack, edgeEnd);
                                stList_append(stack, edgeEnd3);

                                edgeEnd3 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd3));
                            }
                        }
                    }
                    // Else is a bridge edge end, so on a connected path
                    else {
                        onConnectedPath[0] = 1;
                    }
                }
            }
        }
        else {
            // If on a connected path add to the bridge node's connected nodes.
            if(onConnectedPath[0]) {
                stSet_insert(bridgeNode->connectedCactusNodes, cactusNode);

                // Search for the predecessor and update its connected path variable to
                // indicate it is also one a connected path
                onConnectedPath = stHash_search(edgeEndsToConnectedStatus, pEdgeEnd);
                assert(onConnectedPath != NULL);
                onConnectedPath[0] = 1;
            }
        }
    }

    // Cleanup
    stList_destruct(stack);
    stHash_destruct(edgeEndsToConnectedStatus);
}

static void getBridges(stBridgeNode *bridgeNode) {
    /*
     * Compute set of bridge edge ends to bridgeNode->bridgeEnds.
     */

    // For each cactus node in bridgeNode->cactusNodes
    stSetIterator *nodeIt = stSet_getIterator(bridgeNode->cactusNodes);
    stCactusNode *cactusNode;
    while((cactusNode = stSet_getNext(nodeIt)) != NULL) {

        // For each incident edge end
        stCactusNodeEdgeEndIt edgeIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *edgeEnd;
        while((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIt)) != NULL) {
            if(stCactusEdgeEnd_getLink(edgeEnd) == NULL) {

                // Add to the set
                stSet_insert(bridgeNode->bridgeEnds, edgeEnd);
            }
        }
    }
    stSet_destructIterator(nodeIt);
}

void stBridgeNode_print(stBridgeNode *bridgeNode, FILE *fileHandle) {
    // Print nodes and edges indicating bridges and chain edges.
    fprintf(fileHandle, "Bridge node graph: %" PRIi64 " nodes, %" PRIi64 " bridges, %" PRIi64 " bridge connecting nodes\n",
            stSet_size(bridgeNode->cactusNodes),
            stSet_size(bridgeNode->bridgeEnds),
            stSet_size(bridgeNode->connectedCactusNodes));
    stSetIterator *nodeIt = stSet_getIterator(bridgeNode->cactusNodes);
    stCactusNode *cactusNode;
    while((cactusNode = stSet_getNext(nodeIt)) != NULL) {
        fprintf(fileHandle, "Node: %" PRIi64 " in connected nodes: %i, connected to:\n",
                (int64_t)cactusNode, stSet_search(bridgeNode->connectedCactusNodes, cactusNode) != NULL);
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *edgeEnd;
        while((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            assert(stCactusEdgeEnd_getNode(edgeEnd) == cactusNode);
            fprintf(fileHandle, "\tEdge: %" PRIi64 "\tis_bridge: %i\n",
                    (int64_t)stCactusEdgeEnd_getNode(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd)),
                    stCactusEdgeEnd_getLink(edgeEnd) == NULL);
        }
    }
    stSet_destructIterator(nodeIt);
}

stBridgeGraph *stBridgeGraph_getBridgeGraph(stCactusNode *cactusNode) {
    // Create the bridge graph container
    stBridgeGraph *bridgeGraph = st_malloc(sizeof(stBridgeGraph));
    bridgeGraph->bridgeNodes = stList_construct3(0, (void (*)(void *))stBridgeNode_destruct);

    // Get cactus nodes
    stSet *component = stSet_construct();
    buildComponent(component, cactusNode, 0);

    stSet *seen = stSet_construct(); // A record of which nodes are in a bridge node;

    stHash *bridgeEndsToBridgeNodes = stHash_construct();

    // For each cactus node
    stSetIterator *componentIt = stSet_getIterator(component);
    stCactusNode *cactusNode2;
    while((cactusNode2 = stSet_getNext(componentIt)) != NULL) {

        // If not already in a bridge graph node
        if(stSet_search(seen, cactusNode2) == NULL) {

            // Make a bridge node
            stBridgeNode *bridgeNode = stBridgeNode_construct();

            // Add to bridge nodes
            stList_append(bridgeGraph->bridgeNodes, bridgeNode);

            // Get nodes in the bridge component
            buildComponent(bridgeNode->cactusNodes, cactusNode2, 1);

            // Add nodes in bridge component to seen
            stSet_insertAll(seen, bridgeNode->cactusNodes);
            assert(stSet_search(bridgeNode->cactusNodes, cactusNode2) == cactusNode2);

            // Iterate to find bridge edges
            getBridges(bridgeNode);

            // For each bridge end recursively search to find connected nodes
            // and add the bridge end to the bridgeEndsToBridgeNodes map
            stSetIterator *bridgeEndIt = stSet_getIterator(bridgeNode->bridgeEnds);
            stCactusEdgeEnd *bridgeEnd;
            while((bridgeEnd = stSet_getNext(bridgeEndIt)) != NULL) {
                buildBridgeComponent(bridgeNode, bridgeEnd);
                stHash_insert(bridgeEndsToBridgeNodes, bridgeEnd, bridgeNode);
            }
            stSet_destructIterator(bridgeEndIt);
        }
    }
    assert(stSet_size(component) == stSet_size(seen));
    stSet_destructIterator(componentIt);

    // Link the bridge nodes together

    // For each bridge node
    for(int64_t i=0; i<stList_length(bridgeGraph->bridgeNodes); i++) {
        stBridgeNode *bridgeNode = stList_get(bridgeGraph->bridgeNodes, i);

        // For each bridge end in the node
        stSetIterator *bridgeEndIt = stSet_getIterator(bridgeNode->bridgeEnds);
        stCactusEdgeEnd *bridgeEnd;
        while((bridgeEnd = stSet_getNext(bridgeEndIt)) != NULL) {

            // Find the connected bridge node
            stBridgeNode *bridgeNode2 = stHash_search(bridgeEndsToBridgeNodes,
                    stCactusEdgeEnd_getOtherEdgeEnd(bridgeEnd));
            assert(bridgeNode2 != NULL);

            // Attach to the bridge node
            stList_append(bridgeNode->connectedNodes, bridgeNode2);
        }
        stSet_destructIterator(bridgeEndIt);
    }

    // Cleanup
    stSet_destruct(seen);
    stSet_destruct(component);
    stHash_destruct(bridgeEndsToBridgeNodes);

    return bridgeGraph;
}

void stBridgeGraph_destruct(stBridgeGraph *bridgeGraph) {
    stList_destruct(bridgeGraph->bridgeNodes);
    free(bridgeGraph);
}

/*
 * Core functions for computing ultrabubbles for a cactus graph.
 */

stUltraBubble *stUltraBubble_construct(stList *parentChain,
        stCactusEdgeEnd *edgeEnd1, stCactusEdgeEnd *edgeEnd2) {
    stUltraBubble *ultraBubble = st_malloc(sizeof(stUltraBubble));
    if(parentChain != NULL) {
        stList_append(parentChain, ultraBubble);
    }
    ultraBubble->chains = stList_construct3(0,
            (void (*)(void *))stList_destruct);
    ultraBubble->edgeEnd1 = edgeEnd1;
    ultraBubble->edgeEnd2 = edgeEnd2;
    return ultraBubble;
}

void stUltraBubble_destruct(stUltraBubble *ultraBubble) {
    stList_destruct(ultraBubble->chains);
    free(ultraBubble);
}

/*
 * Functions for printing ultra bubble decomposition.
 */

void stUltraBubble_print2(stUltraBubble *ultraBubble, FILE *fileHandle, const char *parentPrefix, const char *parentChainPrefix);

void stUltraBubble_printChains2(stList *ultraBubbleChains, FILE *fileHandle, const char *parentPrefix, const char *parentChainPrefix) {
    for(int64_t i=0; i<stList_length(ultraBubbleChains); i++) {
        stList *chain = stList_get(ultraBubbleChains, i);
        char *chainPrefix = stString_print("%s-%" PRIi64 "", parentChainPrefix, i);
        fprintf(fileHandle, "%s\tChain: %s\tLength: %" PRIi64 "\n", parentPrefix, chainPrefix, stList_length(chain));
        for(int64_t j=0; j<stList_length(chain); j++) {
            char *prefix = stString_print("%s%s", parentPrefix, "\t\t");
            stUltraBubble_print2(stList_get(chain, j), fileHandle, prefix, chainPrefix);
            free(prefix);
        }
        free(chainPrefix);
    }
}

void stUltraBubble_print2(stUltraBubble *ultraBubble, FILE *fileHandle, const char *parentPrefix, const char *parentChainPrefix) {
    fprintf(fileHandle, "%sBubble,\tchild-chain-number: %" PRIi64 "\tSide1:(%" PRIi64 ",%" PRIi64 ")\tSide2:(%" PRIi64 ",%" PRIi64 ")\n", parentPrefix,
            stList_length(ultraBubble->chains), (int64_t)stCactusEdgeEnd_getOtherEdgeEnd(ultraBubble->edgeEnd1), (int64_t)ultraBubble->edgeEnd1,
            (int64_t)ultraBubble->edgeEnd2, (int64_t)stCactusEdgeEnd_getOtherEdgeEnd(ultraBubble->edgeEnd2));
    stUltraBubble_printChains2(ultraBubble->chains, fileHandle, parentPrefix, parentChainPrefix);
}

void stUltraBubble_print(stUltraBubble *ultraBubble, FILE *fileHandle) {
    stUltraBubble_print2(ultraBubble, fileHandle, "", "C");
}

void stUltraBubble_printChains(stList *ultraBubbleChains, FILE *fileHandle) {
    fprintf(fileHandle, "Top level chains, total: %" PRIi64 "\n", stList_length(ultraBubbleChains));
    stUltraBubble_printChains2(ultraBubbleChains, fileHandle, "", "C");
}

/*
 * Functions for building ultrabubble decomposition
 */

static stList *addNewNestedChain(stList *chains) {
    // Adds an ultrabubble chain to a list of ultrabubble chains.
    stList *chain = stList_construct3(0, (void (*)(void *))stUltraBubble_destruct);
    stList_append(chains, chain);
    return chain;
}

static void recursiveEnumerateChainBubbles(stList *parentUltraBubbleChain, stCactusNode *parentNode, stCactusEdgeEnd *parentChain) {
    // Stack to keep track of recursion
    stList *stack = stList_construct();
    stList_append(stack, parentUltraBubbleChain);
    stList_append(stack, parentNode);
    stList_append(stack, parentChain);

    while(stList_length(stack) > 0) {
        // Pop the stack
        parentChain = stList_pop(stack);
        parentNode = stList_pop(stack);
        parentUltraBubbleChain = stList_pop(stack);

        // Checks
        assert(stCactusEdgeEnd_getNode(parentChain) == parentNode);
        assert(stCactusEdgeEnd_getNode(stCactusEdgeEnd_getLink(parentChain)) == parentNode);
        assert(stCactusEdgeEnd_getLinkOrientation(parentChain));
        assert(!stCactusEdgeEnd_getLinkOrientation(stCactusEdgeEnd_getLink(parentChain)));

        // Create a bubble and append it to the parent's ultrabubble chain
        stUltraBubble *ultraBubble = stUltraBubble_construct(parentUltraBubbleChain,
                parentChain, stCactusEdgeEnd_getLink(parentChain));

        // For each child chain incident with parentNode and not equal to parentChain
        stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(parentNode);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
            assert(parentNode == stCactusEdgeEnd_getNode(edgeEnd));
            assert(stCactusEdgeEnd_getLink(edgeEnd) != NULL); // Can not be a bridge, by definition
            if (stCactusEdgeEnd_getLinkOrientation(edgeEnd) && edgeEnd != parentChain) {

                stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd));

                if(edgeEnd != edgeEnd2) {
                    // Make a new ultrabubble chain
                    stList *childUltraBubbleChain = addNewNestedChain(ultraBubble->chains);

                    // For each node in the chain not equal to parentNode
                    do {

                        // Add to the stack
                        stList_append(stack, childUltraBubbleChain);
                        stList_append(stack, stCactusEdgeEnd_getNode(edgeEnd2));
                        stList_append(stack, edgeEnd2);

                        edgeEnd2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                    } while (edgeEnd != edgeEnd2);
                }
            }
        }
    }

    // Cleanup
    stList_destruct(stack);
}

void bridgeGraphToUltraBubbles(stBridgeNode *pBridgeNode, stBridgeNode *bridgeNode, stSet *seen,
        stList *topLevelChains, stSet *cactusGraphComponent, stCactusNode *startNode) {
    // Stack for recursion
    stList *stack = stList_construct();
    stList_append(stack, pBridgeNode);
    stList_append(stack, bridgeNode);

    // Loop while there are bridge nodes to consider
    while(stList_length(stack) > 0) {
        bridgeNode = stList_pop(stack);
        pBridgeNode = stList_pop(stack);

        // This check ensures we do not loop in the case that the bridge graph contains a cycle
        if(stSet_search(seen, bridgeNode) != NULL) {
            continue;
        }
        stSet_insert(seen, bridgeNode);

        // If the bridge node has incident edges
        if(stSet_size(bridgeNode->bridgeEnds) > 0) {
            stList *workingLevelChains = topLevelChains;

            // If incident edges equals 2
            if(stSet_size(bridgeNode->bridgeEnds) == 2) {

                // Continue the existing top level chain if previous node was link in the chain, else start a new one
                stList *topLevelChain = (pBridgeNode != NULL && stSet_size(pBridgeNode->bridgeEnds) == 2) ?
                        stList_peek(topLevelChains) : addNewNestedChain(topLevelChains);

                // Create a bubble and add it to the top level chain
                stList *bridgeEnds = stSet_getList(bridgeNode->bridgeEnds);

                if(stList_length(topLevelChain) > 0) { // If continuing chain orient
                    // the ends correctly in the chain
                    stUltraBubble *pUltraBubble = stList_peek(topLevelChain);
                    stCactusEdgeEnd *oppEnd = stCactusEdgeEnd_getOtherEdgeEnd(pUltraBubble->edgeEnd2);

                    // This statement is true if the start of the chain was oriented backwards
                    // with respect to this link
                    if(!stList_contains(bridgeEnds, oppEnd)) {
                        assert(stList_length(topLevelChain) == 1);
                        oppEnd = stCactusEdgeEnd_getOtherEdgeEnd(pUltraBubble->edgeEnd1);
                        assert(stList_contains(bridgeEnds, oppEnd));
                        stCactusEdgeEnd *edgeEnd = pUltraBubble->edgeEnd1;
                        pUltraBubble->edgeEnd1 = pUltraBubble->edgeEnd2;
                        pUltraBubble->edgeEnd2 = edgeEnd;
                    }

                    // Reverse the orientation if backwards with respect to the previous
                    // link
                    if(oppEnd == stList_get(bridgeEnds, 1)) {
                        stList_reverse(bridgeEnds);
                    }
                }

                stUltraBubble *ultraBubble = stUltraBubble_construct(topLevelChain,
                            stList_get(bridgeEnds, 0), stList_get(bridgeEnds, 1));

                stList_destruct(bridgeEnds);

                // Set working chains to be those nested in the ultraBubble
                workingLevelChains = ultraBubble->chains;
            }

            // A set to track which nodes we've constructed ultrabubbles for
            stSet *visitedCactusNodes = stSet_construct();

            // For each node in cactusNodes
            stSetIterator *cactusNodesIt = stSet_getIterator(bridgeNode->connectedCactusNodes);
            stCactusNode *cactusNode;
            while((cactusNode = stSet_getNext(cactusNodesIt)) != NULL) {

                // For each chain incident with the cactusNode
                stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(cactusNode);
                stCactusEdgeEnd *edgeEnd;
                while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                    assert(cactusNode == stCactusEdgeEnd_getNode(edgeEnd));
                    if (stCactusEdgeEnd_getLink(edgeEnd) != NULL && stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {

                        // Make new chain
                        stList *ultraBubbleChain = addNewNestedChain(workingLevelChains);

                        // For each node in the chain not in cactusNodes and not yet visited
                        stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd));;
                        while (edgeEnd != edgeEnd2) {
                            stCactusNode *cactusNode2 = stCactusEdgeEnd_getNode(edgeEnd2);
                            assert(cactusNode2 != cactusNode);
                            if(stSet_search(bridgeNode->connectedCactusNodes, cactusNode2) == NULL) {
                                if(stSet_search(visitedCactusNodes, cactusNode2) == NULL) {

                                    // Call recursivelyEnumerateChainBubbles
                                    recursiveEnumerateChainBubbles(ultraBubbleChain, cactusNode2, edgeEnd2);

                                    // Add to visited
                                    stSet_insert(visitedCactusNodes, cactusNode2);
                                }
                            }
                            else {
                                // As we encountered a node in connected nodes we break the chain if it is already started
                                if(stList_length(ultraBubbleChain) > 0) {
                                    ultraBubbleChain = addNewNestedChain(workingLevelChains);
                                }
                            }
                            edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd2));
                        }

                        // If new chain has zero length delete
                        if(stList_length(ultraBubbleChain) == 0) {
                            assert(stList_peek(workingLevelChains) == ultraBubbleChain);
                            stList_pop(workingLevelChains);
                            stList_destruct(ultraBubbleChain);
                        }
                    }
                }
            }
            stSet_destructIterator(cactusNodesIt);
            stSet_destruct(visitedCactusNodes);

            // Recursively walk the remainder of the bridge graph to construct chains for the other nodes
            for(int64_t i=0; i<stList_length(bridgeNode->connectedNodes); i++) {
                stBridgeNode *nBridgeNode = stList_get(bridgeNode->connectedNodes, i);
                if(nBridgeNode != pBridgeNode) {
                    assert(nBridgeNode != bridgeNode);
                    stList_append(stack, bridgeNode);
                    stList_append(stack, nBridgeNode);
                }
            }
        }

        // Else the bridgeNode is in a circular component
        else {
            // Pick the start node, otherwise the highest degree node in the cactus component
            // Mostly this is just to make things deterministic.
            stSetIterator *componentIt = stSet_getIterator(cactusGraphComponent);
            stCactusNode *cactusNode = stSet_getNext(componentIt);
            int64_t chainNumber = stCactusNode_getChainNumber(cactusNode);
            stCactusNode *cactusNode2;
            while((cactusNode2 = stSet_getNext(componentIt)) != NULL) {
                if(cactusNode2 == startNode || (stCactusNode_getChainNumber(cactusNode2) > chainNumber && cactusNode != startNode)) {
                    cactusNode = cactusNode2;
                    chainNumber = stCactusNode_getChainNumber(cactusNode2);
                }
            }
            stSet_destructIterator(componentIt);

            // For each chain incident with the cactusNode
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(cactusNode);
            stCactusEdgeEnd *edgeEnd;
            while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                assert(cactusNode == stCactusEdgeEnd_getNode(edgeEnd));
                assert(stCactusEdgeEnd_getLink(edgeEnd) != NULL); // Can not be a bridge, by definition
                if (stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {

                    // For each node in the chain not equal to cactusNode
                    stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd));

                    if(edgeEnd != edgeEnd2) { // If not a self loop
                        // Make new chain
                        stList *ultraBubbleChain = addNewNestedChain(topLevelChains);

                        do {
                            assert(stCactusEdgeEnd_getNode(edgeEnd2) != cactusNode);

                            // Call recursivelyEnumerateChainBubbles
                            recursiveEnumerateChainBubbles(ultraBubbleChain, stCactusEdgeEnd_getNode(edgeEnd2), edgeEnd2);

                            edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(edgeEnd2));
                        } while (edgeEnd != edgeEnd2);

                        assert(stList_length(ultraBubbleChain) > 0);
                    }
                }
            }
        }
    }

    // Clean up
    stList_destruct(stack);
}

stBridgeNode *getStartingNode(stBridgeGraph *bridgeGraph) {
    /*
     * Picks a non-chain node (if it exists). A chain node is a node with two incident edges.
     * This ensures that all chains are complete, or if cyclic are broken arbitrarily.
     */
    for(int64_t i=0; i<stList_length(bridgeGraph->bridgeNodes); i++) {
        stBridgeNode *bridgeNode = stList_get(bridgeGraph->bridgeNodes, i);
        if(stList_length(bridgeNode->connectedNodes) != 2) {
            return bridgeNode;
        }
    }
    return stList_peek(bridgeGraph->bridgeNodes);
}

stList *stCactusGraph_getUltraBubbles(stCactusGraph *graph, stCactusNode *startNode) {
    /*
     * Constructs the set of nested ultrabubbles for a cactus graph, returning the
     * set of top level chains as a list. A bubble is top level if it is not
     * contained in any other. A chain is a sequence of bubbles, such that each successive
     * bubble in the chain starts with the opposite side of a side in the previous bubble.
     *
     * The startNode is used as a tip for the ultra-bubble construction, it may be NULL.
     */

    // Create the list to contain the top level chains
    stList *topLevelChains = stList_construct3(0, (void (*)(void *))stList_destruct);

    // Get cactus graph components
    stList *cactusGraphComponents = stCactusGraph_getComponents(graph, 0);

    // For each component in the cactus graph
    for(int64_t i=0; i<stList_length(cactusGraphComponents); i++) {
        stSet *cactusGraphComponent = stList_get(cactusGraphComponents, i);

        // Determine the bridge graph for the component
        stBridgeGraph *bridgeGraph = stBridgeGraph_getBridgeGraph(stSet_peek(cactusGraphComponent));
        assert(stList_length(bridgeGraph->bridgeNodes) > 0);

        // Get a starting node for recursion
        stBridgeNode *bridgeNode = getStartingNode(bridgeGraph);

        // Recursively build the ultra bubble chains for the bridge graph
        stSet *seen = stSet_construct();
        bridgeGraphToUltraBubbles(NULL, bridgeNode, seen, topLevelChains, cactusGraphComponent, startNode);
        assert(stSet_size(seen) == stList_length(bridgeGraph->bridgeNodes));
        stSet_destruct(seen);

        // Cleanup the bridge graph
        stBridgeGraph_destruct(bridgeGraph);
    }

    // Cleanup the cactus graph components
    stList_destruct(cactusGraphComponents);

    return topLevelChains;
}

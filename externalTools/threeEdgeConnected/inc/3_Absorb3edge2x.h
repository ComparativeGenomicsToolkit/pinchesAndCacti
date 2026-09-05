/*
 * Copyright (C) 2011 by Yung H Tsin and Nima Norouzi
 *
 * The authors of this code have given their permission
 * to distribute and redistribute this code within this software product,
 * however they require that any further incorporation or copying
 * be permitted by request only, hence the BSD/MIT license which
 * applies to much of the other code in this software product does not apply to this file.
 */

#ifndef ABSORB_3_EDGE_2X_H_
#define ABSORB_3_EDGE_2X_H_

/*
 * Function takes an graph represented as an adjacency list. Each vertex is represented
 * by a list in the argument "vertices", it's index in the list is its identifier.
 * Each vertice's list contain 1 length stIntTuples that indicate what the vertex is
 * connected to. The return value is a list of lists of nodes,
 * also represented using stLists and stIntTuples.
 */
stList *computeThreeEdgeConnectedComponents(stList *vertices);

/*
 * The same computation on a graph in compressed sparse row form: vertex v's neighbours are
 * adj[offsets[v]] .. adj[offsets[v+1]-1], with n vertices numbered from 0. On return *membersOut
 * holds the n vertices grouped by component and *startsOut the start of each of the *nComponentsOut
 * components in it (with a sentinel at the end); both are malloc'd and belong to the caller. The
 * components come out in the same order, and their members in the same order, as the list interface.
 */
void computeThreeEdgeConnectedComponentsCSR(int n, const int *offsets, const int *adj, int **membersOut, int **startsOut, int *nComponentsOut);

#endif

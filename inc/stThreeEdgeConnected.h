/*
 * stThreeEdgeConnected.h
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_THREE_EDGE_CONNECTED_H_
#define ST_THREE_EDGE_CONNECTED_H_

#ifdef __cplusplus
extern "C"{
#endif

/*
 * Computes the 3-edge-connected components of an undirected multigraph given in compressed sparse
 * row form: vertex v's neighbours are adj[offsets[v]] .. adj[offsets[v+1]-1], with n vertices
 * numbered from 0. Each edge appears once in the list of each of its ends, so a self loop appears
 * twice in its vertex's list. On return *membersOut holds the n vertices grouped by component and
 * *startsOut the start of each of the *nComponentsOut components in it (with a sentinel at the
 * end); both are malloc'd and belong to the caller.
 */
void stThreeEdgeConnected_getComponents(int n, const int *offsets, const int *adj, int **membersOut, int **startsOut, int *nComponentsOut);

#ifdef __cplusplus
}
#endif
#endif

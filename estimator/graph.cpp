#ifndef GRAPH_H
#define GRAPH_H

#include <cstdlib>
#include <vector>

namespace estimator {

	class Edge {
	public:
		Edge(uint32_t from, uint32_t to, float cost) : from(from), to(to), cost(cost) {}

		uint32_t from, to;
		float cost;
	};

	class Graph {
	public:
		Graph(uint32_t size) {
			n = size;
			edges = (std::vector<Edge*>**) malloc(sizeof(std::vector<Edge>*) * size);
			for (int i = 0; i < n; i++) {
				edges[i] = NULL;
			}
		}
		~Graph() {
			for (int i = 0; i < n; i++) {
				if (edges[i] != NULL) {
					for (int j = 0; j < edges[i]->size(); j++) {
						delete edges[i]->at(j);
					}
					delete edges[i];
				}
			}
			free(edges);
		}
		void addEdge(uint32_t from, uint32_t to, float cost) {
			if (edges[from] == NULL) edges[from] = new std::vector<Edge*>();
			if (edges[to] == NULL) edges[to] = new std::vector<Edge*>();
			edges[from]->push_back(new Edge(from, to, cost));
			edges[to]->push_back(new Edge(to, from, cost));
		}
		std::vector<Edge*>* getNeighbors(uint32_t index) {
			return edges[index];
		}
	private:
		uint32_t n;
		std::vector<Edge*>** edges;
	};

}

#endif

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <float.h>
#include <string.h>
#include <set>
#include <math.h>
#include <algorithm>

#include "octree/Octree.hpp"
#include "estimator/graph.cpp"
#include "estimator/point.cpp"

std::vector<estimator::Point> points;
std::vector<estimator::Point> normals;

estimator::Point planeFromPoints(std::vector<uint32_t> neighbors) {
	const size_t n = neighbors.size();
	if (n < 3) throw "at least 3 points required";

	estimator::Point sum(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < neighbors.size(); i++) {
		estimator::Point p = points.at(neighbors.at(i));
		sum.x += p.x;
		sum.y += p.y;
		sum.z += p.z;
	}
	estimator::Point centroid(sum.x / n, sum.y / n, sum.z / n);

	float xx = 0.0, xy = 0.0, xz = 0.0, yy = 0.0, yz = 0.0, zz = 0.0;
	for (int i = 0; i < neighbors.size(); i++) {
		estimator::Point p = points.at(neighbors.at(i));
		estimator::Point r(p.x - centroid.x, p.y - centroid.y, p.z - centroid.z);
		xx += r.x * r.x;
		xy += r.x * r.y;
		xz += r.x * r.z;
		yy += r.y * r.y;
		yz += r.y * r.z;
		zz += r.z * r.z;
	}

	const float detX = (yy * zz) - (yz * yz);
	const float detY = (xx * zz) - (xz * xz);
	const float detZ = (xx * yy) - (xy * xy);

	float detMax = std::max(std::max(detX, detY), detZ);
	if (detMax <= 0) throw "points dont span a plane";

	estimator::Point dir(0.0, 0.0, 0.0);
	if (detMax == detX) {
		dir.x = detX;
		dir.y = (xz * yz) - (xy * zz);
		dir.z = (xy * yz) - (xz * yy);
	}
	else if (detMax == detY) {
		dir.x = (xz * yz) - (xy * zz);
		dir.y = detY;
		dir.z = (xy * xz) - (yz * xx);
	}
	else {
		dir.x = (xy * yz) - (xz * yy);
		dir.y = (xy * xz) - (yz * xx);
		dir.z = detZ;
	}
	dir.normalize();

	return dir;
}

void propagateOrientation(estimator::Point p, uint32_t nextIndex) {
	estimator::Point nextPoint = normals.at(nextIndex);
	estimator::Point invertedPoint(nextPoint.x * -1, nextPoint.y * -1, nextPoint.z * -1);

	if (acos(invertedPoint.dot(p)) < acos(nextPoint.dot(p))) {
		normals.at(nextIndex).flipSign();
	}
}

bool pointSort(estimator::Point a, estimator::Point b) {
	return a.z < b.z;
}

int main(int argc, char *argv[]) {

	// Read points in format X Y Z per line from STDIN
	float p[3];
	while (scanf("%f %f %f", &p[0], &p[1], &p[2]) == 3) {
		points.push_back(estimator::Point(p[0], p[1], p[2]));
	}

	// Sort the points so that the highest Z value is at the end
	std::sort(points.begin(), points.end(), pointSort);

	// Build an octree to search on the points
	unibn::Octree<estimator::Point> octree;
	unibn::OctreeParams params;
	octree.initialize(points);

	// Estimate a tangent plane for each point to serve as a normal
	for (uint32_t i = 0; i < points.size(); ++i) {
		estimator::Point point = points.at(i);
		std::vector<uint32_t> results;
		octree.radiusNeighbors<unibn::L2Distance<estimator::Point> >(point, 5.0f, results);
		estimator::Point normal = planeFromPoints(results);
		normals.push_back(normal);
	}

	// Build Riemmanian Graph from Euclidean Minimum Spanning Tree
	estimator::Graph* riemann = new estimator::Graph(points.size());
	uint32_t curIndex = floor(points.size() / 2);
	// Visited points
	std::set<uint32_t> emstSet;
	emstSet.insert(curIndex);
	// Potential edges to add to graph
	std::vector<estimator::Edge*>* pEdges = new std::vector<estimator::Edge*>();
	// Visit every node and create n-1 edges total
	while (emstSet.size() != points.size()) {
		float minDist = FLT_MAX;
		uint32_t minIndex = -1;
		float searchRadius = 1.0f;
		std::vector<estimator::Edge*> localEdges;
		std::vector<uint32_t> neighbors;
		// Search until we find atleast 1 unvisited neighbor point
		while (localEdges.size() < 1) {
			octree.radiusNeighbors<unibn::L2Distance<estimator::Point> >(points.at(curIndex), searchRadius, neighbors);
			for (int i = 0; i < neighbors.size(); i++) {
				uint32_t idx = neighbors.at(i);
				std::set<uint32_t>::iterator find = emstSet.find(idx);
				bool visited = find != emstSet.end();
				if (idx != curIndex && !visited) {
					float dist = std::sqrt(unibn::L2Distance<estimator::Point>::compute(points.at(curIndex), points.at(idx)));
					estimator::Edge* candidateEdge = new estimator::Edge(curIndex, idx, dist);
					localEdges.push_back(candidateEdge);
				}
			}
			// Increase search radius until we get a result
			searchRadius += 0.5f;
			neighbors.clear();
		}

		// Add new edges to potential edges
		pEdges->insert(pEdges->end(), localEdges.begin(), localEdges.end());
		// Find potential edge with lost cost
		for (auto i = pEdges->begin(); i != pEdges->end();) {
			estimator::Edge* edge = *i;
			std::set<uint32_t>::iterator find = emstSet.find(edge->to);
			bool visited = find != emstSet.end();
			if (visited) {
				free(*i);
				i = pEdges->erase(i);
			} else {
				if (edge->cost < minDist) {
					minDist = edge->cost;
					minIndex = i - pEdges->begin();
				}
				++i;
			}
		}
		estimator::Edge* newEdge = pEdges->at(minIndex);
		
		// Add new edge to riemannian graph using cosine similarity of normals as cost
		riemann->addEdge(newEdge->from, newEdge->to, acos(normals.at(newEdge->from).dot(normals.at(newEdge->to))));
		pEdges->erase(pEdges->begin() + minIndex);

		curIndex = newEdge->to;
		emstSet.insert(curIndex);
		free(newEdge);
	}
	delete pEdges;

	// Add additional edges to reduce tree sparsity
	for (uint32_t i = 0; i < points.size(); i++) {
		std::vector<uint32_t> neighbors;
		float searchRadius = 0.5f;
		while (neighbors.size() < 6) {
			octree.radiusNeighbors<unibn::L2Distance<estimator::Point> >(points.at(i), searchRadius, neighbors);
			searchRadius += 0.1f;
		}
		for (int j = 0; j < neighbors.size(); j++) {
			uint32_t idx = neighbors.at(j);
			if (idx != i) {
				riemann->addEdge(i, idx, acos(normals.at(i).dot(normals.at(idx))));
			}
		}
	}

	// Sign the largest Z point with Z+
	estimator::Point zPlus(0.0f, 0.0f, 1.0f);
	propagateOrientation(zPlus, normals.size()-1);
	// Traverse the MST of the riemmanian graph using normal cosine similarity as cost and propagate the sign
	std::set<uint32_t> mstSet;
	mstSet.insert(points.size() - 1);
	pEdges = new std::vector<estimator::Edge*>();
	std::vector<estimator::Edge*>* seedNeighbors = riemann->getNeighbors(points.size() - 1);
	pEdges->insert(pEdges->end(), seedNeighbors->begin(), seedNeighbors->end());
	while (mstSet.size() != points.size()) {
		uint32_t from = -1;
		uint32_t to = -1;
		float minCost = FLT_MAX;

		for (auto i = pEdges->begin(); i != pEdges->end();) {
			estimator::Edge* edge = *i;
			std::set<uint32_t>::iterator find = mstSet.find(edge->to);
			bool visited = find != mstSet.end();
			if (visited) {
				i = pEdges->erase(i);
			} else {
				if (edge->cost < minCost) {
					from = edge->from;
					to = edge->to;
					minCost = edge->cost;
				}
				++i;
			}
		}

		propagateOrientation(normals.at(from), to);
		mstSet.insert(to);
		std::vector<estimator::Edge*>* neighbors = riemann->getNeighbors(to);
		pEdges->insert(pEdges->end(), neighbors->begin(), neighbors->end());
	}

	// Output results to STDOUT
	for (uint32_t i = 0; i < points.size(); ++i) {
		estimator::Point point = points.at(i);
		estimator::Point normal = normals.at(i);
		std::cout << point.x << " " << point.y << " " << point.z << " ";
		std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;
	}

	return 0;
}
//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_DIJKSTRA_010506
#define GEODESIC_ALGORITHM_DIJKSTRA_010506

#include "geodesic_algorithm_graph_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic{

class DijkstraNode
{
	typedef DijkstraNode* node_pointer;
public: 
	DijkstraNode(){};
	~DijkstraNode(){};

	double& distance_from_source(){return m_distance;};
	node_pointer& previous(){return m_previous;};
	unsigned& source_index(){return m_source_index;};
	vertex_pointer& vertex(){return m_vertex;};

	void clear()
	{
		m_distance = GEODESIC_INF;
		m_previous = NULL;
	}

	bool operator()(node_pointer const s1, node_pointer const s2) const
	{
		return s1->distance_from_source() != s2->distance_from_source() ?
			   s1->distance_from_source() < s2->distance_from_source() :
			   s1->vertex()->id() < s2->vertex()->id();
	};

	double distance(SurfacePoint* p)
	{
		return m_vertex->distance(p);
	}

	SurfacePoint surface_point()
	{
		return SurfacePoint(m_vertex);
	}

private: 
	double m_distance;					//distance to the closest source
	unsigned m_source_index;			//closest source index
	node_pointer m_previous;			//previous node in the geodesic path
	vertex_pointer m_vertex;			//correspoding vertex
};

class GeodesicAlgorithmDijkstra: public GeodesicAlgorithmGraphBase<DijkstraNode>
{
public:
	typedef DijkstraNode Node;
	typedef Node* node_pointer;

	GeodesicAlgorithmDijkstra(geodesic::Mesh* mesh):
		GeodesicAlgorithmGraphBase<Node>(mesh)
	{
		m_type = DIJKSTRA;

		m_nodes.resize(mesh->vertices().size());
		for(unsigned i=0; i<m_nodes.size(); ++i)
		{
			m_nodes[i].vertex() = &m_mesh->vertices()[i];
		}
	};

	~GeodesicAlgorithmDijkstra(){};

protected:

	void list_nodes_visible_from_source(MeshElementBase* p, 
										std::vector<node_pointer>& storage);		//list all nodes that belong to this mesh element

	void list_nodes_visible_from_node(node_pointer node,			//list all nodes that belong to this mesh element
									  std::vector<node_pointer>& storage,
									  std::vector<double>& distances, 
									  double threshold_distance);	//list only the nodes whose current distance is larger than the threshold
};

void GeodesicAlgorithmDijkstra::list_nodes_visible_from_source(MeshElementBase* p,
															   std::vector<node_pointer>& storage)
{
	assert(p->type() != UNDEFINED_POINT);

	if(p->type() == FACE)
	{
		face_pointer f = static_cast<face_pointer>(p);
		for(unsigned i=0; i<3; ++i)
		{
			vertex_pointer v = f->adjacent_vertices()[i];
			storage.push_back(&m_nodes[node_index(v)]);
		}
	}
	else if(p->type() == EDGE)
	{
		edge_pointer e = static_cast<edge_pointer>(p);
		for(unsigned i=0; i<2; ++i)
		{
			vertex_pointer v = e->adjacent_vertices()[i];
			storage.push_back(&m_nodes[node_index(v)]);
		}
	}
	else			//VERTEX
	{
		vertex_pointer v = static_cast<vertex_pointer>(p);
		storage.push_back(&m_nodes[node_index(v)]);
	}
}

inline void GeodesicAlgorithmDijkstra::list_nodes_visible_from_node(node_pointer node, //list all nodes that belong to this mesh element
															 std::vector<node_pointer>& storage,
															 std::vector<double>& distances,
															 double threshold_distance)
{
	vertex_pointer v = node->vertex();
	assert(storage.size() == distances.size());

	for(unsigned i=0; i<v->adjacent_edges().size(); ++i)
	{
		edge_pointer e = v->adjacent_edges()[i];
		vertex_pointer new_v = e->opposite_vertex(v);
		node_pointer new_node = &m_nodes[node_index(new_v)];

		if(new_node->distance_from_source() > threshold_distance + e->length())
		{
			storage.push_back(new_node);
			distances.push_back(e->length());
		}
	}
}

}		//geodesic

#endif //GEODESIC_ALGORITHM_DIJKSTRA_010506

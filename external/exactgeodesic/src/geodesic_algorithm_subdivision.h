//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_SUBDIVISION_122806
#define GEODESIC_ALGORITHM_SUBDIVISION_122806

#include "geodesic_algorithm_graph_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic{

class SubdivisionNode: public SurfacePoint
{
	typedef SubdivisionNode* node_pointer;
public: 
	SubdivisionNode(){};

	template <class Pointer>
	SubdivisionNode(Pointer p):
		SurfacePoint(p),
		m_previous(NULL),
		m_distance(0.0)
	{};

	template <class Pointer, class Parameter>
	SubdivisionNode(Pointer p, Parameter param):
		SurfacePoint(p, param),
		m_previous(NULL),
		m_distance(0.0)
	{};

	~SubdivisionNode(){};

	double& distance_from_source(){return m_distance;};
	node_pointer& previous(){return m_previous;};
	unsigned& source_index(){return m_source_index;};

	void clear()
	{
		m_distance = GEODESIC_INF;
		m_previous = NULL;
	}

	bool operator()(node_pointer const s1, node_pointer const s2) const
	{
		if(s1 == s2)
		{
			return false;
		}
		if(s1->distance_from_source() != s2->distance_from_source())
		{
			return s1->distance_from_source() < s2->distance_from_source();
		}
/*		if(s1->type() != s2->type())
		{
			return s1->type() < s2->type();
		}
		if(s1->base_element()->id() != s2->base_element()->id())
		{
		    return s1->base_element()->id() < s2->base_element()->id();
		} */
		if(s1->x() != s2->x())		//two nodes cannot be located in the same space 
		{
			return s1->x() < s2->x();
		}
		if(s1->y() != s2->y())		
		{
			return s1->y() < s2->y();
		}
		if(s1->z() != s2->z())		
		{
			return s1->z() < s2->z();
		}

		assert(0);
		return true;
	};

	SurfacePoint& surface_point(){return static_cast<SurfacePoint&>(*this);};

private: 
	double m_distance;					//distance to the closest source
	unsigned m_source_index;			//closest source index
	node_pointer m_previous;			//previous node in the geodesic path
};

class GeodesicAlgorithmSubdivision: public GeodesicAlgorithmGraphBase<SubdivisionNode>
{
	typedef SubdivisionNode Node;
public:
    GeodesicAlgorithmSubdivision() {}
	GeodesicAlgorithmSubdivision(geodesic::Mesh* mesh, 
								 unsigned subdivision_level):
		GeodesicAlgorithmGraphBase<Node>(mesh)
	{
		m_type = SUBDIVISION;

		m_nodes.reserve(mesh->vertices().size());
		for(unsigned i=0; i<mesh->vertices().size(); ++i)
		{
			vertex_pointer v = &mesh->vertices()[i];
			m_nodes.push_back(Node(v));		//!!
		}

		set_subdivision_level(subdivision_level);
	};

	~GeodesicAlgorithmSubdivision(){};

	unsigned subdivision_level(){return m_subdivision_level;};

	void set_subdivision_level(unsigned subdivision_level)
	{
		m_subdivision_level = subdivision_level;

		m_nodes.resize(m_mesh->vertices().size());
		m_nodes.reserve(m_mesh->vertices().size() + 
						m_mesh->edges().size()*subdivision_level);

		for(unsigned i=0; i<m_mesh->edges().size(); ++i)
		{
			edge_pointer e = &m_mesh->edges()[i];
			for(unsigned i=0; i<subdivision_level; ++i)
			{
				double offset = (double)(i+1)/(double)(subdivision_level+1);
				m_nodes.push_back(Node(e, offset));
			}
		}
	};

protected:
	void list_nodes_visible_from_source(MeshElementBase* p, 
										std::vector<node_pointer>& storage);		//list all nodes that belong to this mesh element

	void list_nodes_visible_from_node(node_pointer node,			//list all nodes that belong to this mesh element
									  std::vector<node_pointer>& storage,
									  std::vector<double>& distances, 
									  double threshold_distance);	//list only the nodes whose current distance is larger than the threshold
	
	unsigned node_indexx(edge_pointer e)
	{
		return e->id()*m_subdivision_level + m_mesh->vertices().size();
	};

private:
	void list_nodes(MeshElementBase* p,		//list nodes that belong to this mesh element
					std::vector<node_pointer>& storage,
					double threshold_distance = -1.0);				//list only the nodes whose current distance is larger than the threshold

	unsigned m_subdivision_level;	//when level is equal to 1, this algorithm corresponds to the Dijkstra algorithm
};

inline void GeodesicAlgorithmSubdivision::list_nodes(MeshElementBase* p,
											         std::vector<node_pointer>& storage,
													 double threshold_distance)
{
	assert(p->type() != UNDEFINED_POINT);

	if(p->type() == VERTEX)
	{
		vertex_pointer v = static_cast<vertex_pointer>(p);
		node_pointer node = &m_nodes[node_index(v)];
		if(node->distance_from_source() > threshold_distance)
		{
			storage.push_back(node);
		}
	}
	else if(p->type() == EDGE)
	{
		edge_pointer e = static_cast<edge_pointer>(p);
		unsigned node_index = node_indexx(e);
		for(unsigned i=0; i<m_subdivision_level; ++i)
		{
			node_pointer node = &m_nodes[node_index++];
			if(node->distance_from_source() > threshold_distance)
			{
				storage.push_back(node);
			}
		}
	}
	//FACE has no nodes
}

inline void GeodesicAlgorithmSubdivision::list_nodes_visible_from_source(MeshElementBase* p,
																  std::vector<node_pointer>& storage)
{
	assert(p->type() != UNDEFINED_POINT);

	if(p->type() == FACE)
	{
		face_pointer f = static_cast<face_pointer>(p);
		for(unsigned i=0; i<3; ++i)
		{
			list_nodes(f->adjacent_vertices()[i],storage);
			list_nodes(f->adjacent_edges()[i],storage);
		}
	}
	else if(p->type() == EDGE)
	{
		list_nodes(p,storage);
		list_nodes(p->adjacent_vertices()[0],storage);
		list_nodes(p->adjacent_vertices()[1],storage);
	}
	else			//VERTEX
	{
		list_nodes(p,storage);
	}
}

inline void GeodesicAlgorithmSubdivision::list_nodes_visible_from_node(node_pointer node, //list all nodes that belong to this mesh element
																std::vector<node_pointer>& storage,
																std::vector<double>& distances,
																double threshold_distance)
{
	MeshElementBase* p = node->base_element();
	assert(p->type() != UNDEFINED_POINT);
	assert(storage.size() == distances.size());

	if(p->type() == VERTEX)
	{
		vertex_pointer v = static_cast<vertex_pointer>(p);

		for(unsigned i=0; i<v->adjacent_edges().size(); ++i)
		{
			edge_pointer e = v->adjacent_edges()[i];
			vertex_pointer v_opposite = e->opposite_vertex(v);
			list_nodes(e, storage, threshold_distance);
			list_nodes(v_opposite, storage, threshold_distance);
		}
		for(unsigned i=0; i<v->adjacent_faces().size(); ++i)
		{
			face_pointer f = v->adjacent_faces()[i];
			edge_pointer e = f->opposite_edge(v);
			list_nodes(e, storage, threshold_distance);
		}
	}
	else if(p->type() == EDGE)
	{
		edge_pointer e = static_cast<edge_pointer>(p);

		vertex_pointer v0 = e->adjacent_vertices()[0];
		vertex_pointer v1 = e->adjacent_vertices()[1];
		list_nodes(v0, storage, threshold_distance);
		list_nodes(v1, storage, threshold_distance);

		for(unsigned i=0; i<e->adjacent_faces().size(); ++i)
		{
			face_pointer f = e->adjacent_faces()[i];

			list_nodes(f->next_edge(e,v0), storage, threshold_distance);
			list_nodes(f->next_edge(e,v1), storage, threshold_distance);
			list_nodes(f->opposite_vertex(e), storage, threshold_distance);
		}
	}
	else 
	{
		assert(0);
	}

	unsigned index = distances.size();
	distances.resize(storage.size());
	for(; index<storage.size(); ++index)
	{	
		distances[index] = node->distance(&storage[index]->surface_point());
	}
}

}		//geodesic

#endif //GEODESIC_ALGORITHM_SUBDIVISION_122806

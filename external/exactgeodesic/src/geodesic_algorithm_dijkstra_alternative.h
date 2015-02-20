//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_DIJKSTRA_ALTERNATIVE_010506
#define GEODESIC_ALGORITHM_DIJKSTRA_ALTERNATIVE_010506

#include "geodesic_algorithm_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic{

class DijkstraNode1
{
	typedef DijkstraNode1* node_pointer;
public: 
	DijkstraNode1(){};
	~DijkstraNode1(){};

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

private: 
	double m_distance;					//distance to the closest source
	unsigned m_source_index;			//closest source index
	node_pointer m_previous;			//previous node in the geodesic path
	vertex_pointer m_vertex;			//correspoding vertex
};

class GeodesicAlgorithmDijkstraAlternative: public GeodesicAlgorithmBase
{
public:
	typedef DijkstraNode1 Node;
	typedef Node* node_pointer;

	GeodesicAlgorithmDijkstraAlternative(geodesic::Mesh* mesh = NULL):
		GeodesicAlgorithmBase(mesh),
		m_nodes(mesh->vertices().size())
	{
		m_type = DIJKSTRA;
		for(unsigned i=0; i<m_nodes.size(); ++i)
		{
			m_nodes[i].vertex() = &m_mesh->vertices()[i];
		}
	};

	~GeodesicAlgorithmDijkstraAlternative(){};

	virtual void propagate(std::vector<SurfacePoint>& sources,
   						   double max_propagation_distance  = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
						   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

	virtual void trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
							std::vector<SurfacePoint>& path);

	virtual unsigned best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
									double& best_source_distance); 
private:

	void set_sources(std::vector<SurfacePoint>& sources)
	{
		m_sources = sources;
	}

	bool check_stop_conditions();
	
	std::vector<Node> m_nodes;	//nodes of the subdivision graph located on the vertices

	std::vector<SurfacePoint> m_sources;		//for simplicity, we keep sources as they are

	typedef std::set<node_pointer, Node> queue_type;
	queue_type m_queue;
};

inline unsigned GeodesicAlgorithmDijkstraAlternative::best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
																	   double& best_source_distance)
{
	if(point.type() == VERTEX)		
	{
		vertex_pointer v = (vertex_pointer)point.base_element();
		node_pointer node = &m_nodes[v->id()];
		best_source_distance = node->distance_from_source();
		return node->source_index();
	}
	else
	{
		std::vector<vertex_pointer> possible_vertices;		//initialize vertices directly visible from sources
		m_mesh->closest_vertices(&point, &possible_vertices);

		best_source_distance = GEODESIC_INF;
		vertex_pointer min_vertex = NULL;
		for(unsigned i=0; i<possible_vertices.size(); ++i)
		{
			vertex_pointer v = possible_vertices[i];

			double distance_from_source = m_nodes[v->id()].distance_from_source();
			double distance_from_dest = v->distance(&point);
			if(distance_from_source + distance_from_dest < best_source_distance)
			{
				best_source_distance = distance_from_source + distance_from_dest;
				min_vertex = v;
			}
		}
		assert(min_vertex);
		node_pointer node = &m_nodes[min_vertex->id()];
		return node->source_index();
	}
}

inline void GeodesicAlgorithmDijkstraAlternative::trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
														   std::vector<SurfacePoint>& path) 
{
	path.clear();
	path.push_back(destination);

	if(destination.type() != VERTEX)		//snap to the closest vertex
	{
		std::vector<vertex_pointer> possible_vertices;		//initialize vertices directly visible from sources
		m_mesh->closest_vertices(&destination, &possible_vertices);

		double min_distance = GEODESIC_INF;
		vertex_pointer min_vertex = NULL;
		for(unsigned i=0; i<possible_vertices.size(); ++i)
		{
			vertex_pointer v = possible_vertices[i];

			double distance_from_source = m_nodes[v->id()].distance_from_source();
			double distance_from_dest = v->distance(&destination);
			if(distance_from_source + distance_from_dest < min_distance)
			{
				min_distance = distance_from_source + distance_from_dest;
				min_vertex = v;
			}
		}
		assert(min_vertex);
		path.push_back(SurfacePoint(min_vertex));
	}

	node_pointer node = &m_nodes[path.back().base_element()->id()];
	while(node->previous())			//follow the path
	{
		node = node->previous();
		path.push_back(SurfacePoint(node->vertex()));
	}

	SurfacePoint& source = m_sources[node->source_index()];		//add source to the path if it is not already there
	if(source.type() != VERTEX ||
	   source.base_element()->id() != node->vertex()->id())
	{
		path.push_back(source);
	}
}

inline void GeodesicAlgorithmDijkstraAlternative::propagate(std::vector<SurfacePoint>& sources,
   															double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
															std::vector<SurfacePoint>* stop_points)   //or after ensuring that all the stop_points are covered
{
	set_stop_conditions(stop_points, max_propagation_distance);
	set_sources(sources);
	
	m_queue.clear();									//clear everything
	for(unsigned i=0; i<m_nodes.size(); ++i)
	{
		m_nodes[i].clear();
	}

	clock_t start = clock();

	std::vector<vertex_pointer> visible_vertices;		//initialize vertices directly visible from sources
	for(unsigned i=0; i<m_sources.size(); ++i)
	{
		SurfacePoint* s = &m_sources[i];
		m_mesh->closest_vertices(s, &visible_vertices);
		for(unsigned j=0; j<visible_vertices.size(); ++j)
		{
			vertex_pointer v = visible_vertices[j];
			double distance = s->distance(v);

			node_pointer node = &m_nodes[v->id()]; 
			if(distance < node->distance_from_source())
			{
				node->distance_from_source() = distance;
				node->source_index() = i;
				node->previous() = NULL;
			}
		}
		visible_vertices.clear();
	}

	for(unsigned i=0; i<m_nodes.size(); ++i)		//initialize the queue
	{
		if(m_nodes[i].distance_from_source() < GEODESIC_INF)
		{
			m_queue.insert(&m_nodes[i]);
		}
	}

	unsigned counter = 0;
	
	while(!m_queue.empty())					//main cycle
	{
		if(counter++ % 10 == 0)		//check if we covered all required vertices
		{
			if (check_stop_conditions())
			{
				break;
			}
			//std::cout << counter << " " << m_queue.size() << " " << (*m_queue.begin())->distance_from_source()<< std::endl;
		}

		node_pointer min_node = *m_queue.begin();
		m_queue.erase(m_queue.begin());

		vertex_pointer v = min_node->vertex();

		for(unsigned i=0; i<v->adjacent_edges().size(); ++i)		//update all the adgecent vertices
		{
			edge_pointer e = v->adjacent_edges()[i];
			vertex_pointer next_v = e->opposite_vertex(v);
			node_pointer next_node = &m_nodes[next_v->id()];
			//double current_distance = 
			if(next_node->distance_from_source() > min_node->distance_from_source() + e->length())
			{
				if(next_node->distance_from_source() < GEODESIC_INF)		//remove it from the queue
				{
					queue_type::iterator iter = m_queue.find(next_node);
					assert(iter != m_queue.end());
					m_queue.erase(iter);
				}
				next_node->distance_from_source() = min_node->distance_from_source() + e->length();
				next_node->source_index() = min_node->source_index();
				next_node->previous() = min_node;
				m_queue.insert(next_node);
			}
		}
	} 
	//std::cout << std::endl;

	clock_t finish = clock();
	m_time_consumed = (static_cast<double>(finish)-static_cast<double>(start))/CLOCKS_PER_SEC;
}

inline bool GeodesicAlgorithmDijkstraAlternative::check_stop_conditions()
{
	double queue_min_distance = (*m_queue.begin())->distance_from_source();

	if(queue_min_distance > m_max_propagation_distance)
	{
		return true;
	}

	for (unsigned index = 0; index < m_stop_vertices.size(); ++index)
	{
		vertex_pointer v = m_stop_vertices[index].first;
		Node& node = m_nodes[v->id()];
		if(queue_min_distance > node.distance_from_source() + m_stop_vertices[index].second)
		{
			return true;
		}
	}
	return false;
}

}		//geodesic

#endif //GEODESIC_ALGORITHM_DIJKSTRA_ALTERNATIVE_010506

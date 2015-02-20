#ifndef GEODESIC_ALGORITHM_GRAPH_BASE_010907
#define GEODESIC_ALGORITHM_GRAPH_BASE_010907

#include "geodesic_algorithm_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic{

template<class Node>
class GeodesicAlgorithmGraphBase: public GeodesicAlgorithmBase
{
public:
	typedef Node* node_pointer;

	GeodesicAlgorithmGraphBase(geodesic::Mesh* mesh = 0):
		GeodesicAlgorithmBase(mesh)
	{};

	~GeodesicAlgorithmGraphBase(){};

	void propagate(std::vector<SurfacePoint>& sources,
   				   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
				   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

	void trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
		std::vector<SurfacePoint>& path);

	unsigned best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
							double& best_source_distance); 

	void print_statistics()
	{
		GeodesicAlgorithmBase::print_statistics();

		double memory = m_nodes.size()*sizeof(Node);
		std::cout << "uses about " << memory/1e6 << "Mb of memory" <<std::endl;
	}

protected:
	unsigned node_index(vertex_pointer v)		//gives index of the node that corresponds to this vertex
	{
		return v->id();
	};

	void set_sources(std::vector<SurfacePoint>& sources)
	{
		m_sources = sources;
	}

	node_pointer best_first_node(SurfacePoint& point, double& best_total_distance)
	{
		node_pointer best_node = NULL;	
		if(point.type() == VERTEX)		
		{
			vertex_pointer v = (vertex_pointer)point.base_element();
			best_node = &m_nodes[node_index(v)];
			best_total_distance = best_node->distance_from_source(); 
		}
		else
		{
			std::vector<node_pointer> possible_nodes;
			list_nodes_visible_from_source(point.base_element(), possible_nodes);

			best_total_distance = GEODESIC_INF;
			for(unsigned i=0; i<possible_nodes.size(); ++i)
			{
				node_pointer node = possible_nodes[i];

				double distance_from_dest = node->distance(&point);
				if(node->distance_from_source() + distance_from_dest < best_total_distance)
				{
					best_total_distance = node->distance_from_source() + distance_from_dest;
					best_node = node;
				}
			}
		}

		//assert(best_node);
		//assert(best_total_distance<GEODESIC_INF);
		if(best_total_distance > m_propagation_distance_stopped)		//result is unreliable
		{
			best_total_distance = GEODESIC_INF;
			return NULL;
		}
		else
		{
			return best_node;
		}
	};	//quickly find what node will be the next one in geodesic path

	bool check_stop_conditions();		//check when propagation should stop

	virtual void list_nodes_visible_from_source(MeshElementBase* p, 
												std::vector<node_pointer>& storage) = 0;		//list all nodes that belong to this mesh element

	virtual void list_nodes_visible_from_node(node_pointer node,			//list all nodes that belong to this mesh element
											  std::vector<node_pointer>& storage,
											  std::vector<double>& distances, 
											  double threshold_distance) = 0;	//list only the nodes whose current distance is larger than the threshold

	std::vector<Node> m_nodes;	//nodes of the graph 

	typedef std::set<node_pointer, Node> queue_type;
	queue_type m_queue;

	std::vector<SurfacePoint> m_sources;		//for simplicity, we keep sources as they are
};

template<class Node>
void GeodesicAlgorithmGraphBase<Node>::propagate(std::vector<SurfacePoint>& sources,
   												 double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
												 std::vector<SurfacePoint>* stop_points) //or after ensuring that all the stop_points are covered
{
	set_stop_conditions(stop_points, max_propagation_distance);
	set_sources(sources);
	
	m_queue.clear();
	m_propagation_distance_stopped = GEODESIC_INF;
	for(unsigned i=0; i<m_nodes.size(); ++i)
	{
		m_nodes[i].clear();
	}

	clock_t start = clock();

	std::vector<node_pointer> visible_nodes;		//initialize vertices directly visible from sources
	for(unsigned i=0; i<m_sources.size(); ++i)
	{
		SurfacePoint* source = &m_sources[i];
		list_nodes_visible_from_source(source->base_element(), 
									   visible_nodes);			

		for(unsigned j=0; j<visible_nodes.size(); ++j)
		{
			node_pointer node = visible_nodes[j];
			double distance = node->distance(source);
			if(distance < node->distance_from_source())
			{
				node->distance_from_source() = distance;
				node->source_index() = i;
				node->previous() = NULL;
			}
		}
		visible_nodes.clear();
	}

	for(unsigned i=0; i<m_nodes.size(); ++i)		//initialize the queue
	{
		if(m_nodes[i].distance_from_source() < GEODESIC_INF)
		{
			m_queue.insert(&m_nodes[i]);
		}
	}

	unsigned counter = 0;
	
	std::vector<double> distances_between_nodes;
	while(!m_queue.empty())					//main cycle
	{
		if(counter++ % 10 == 0)		//check if we covered all required vertices
		{
			if (check_stop_conditions())
			{
				break;
			}
		}

		node_pointer min_node = *m_queue.begin();
		m_queue.erase(m_queue.begin());
		assert(min_node->distance_from_source() < GEODESIC_INF);

		visible_nodes.clear();
		distances_between_nodes.clear();
		list_nodes_visible_from_node(min_node, 
									 visible_nodes, 
									 distances_between_nodes, 
									 min_node->distance_from_source());

		for(unsigned i=0; i<visible_nodes.size(); ++i)		//update all the adgecent vertices
		{
			node_pointer next_node = visible_nodes[i];

			if(next_node->distance_from_source() > min_node->distance_from_source() +
												   distances_between_nodes[i])
			{
				if(next_node->distance_from_source() < GEODESIC_INF)		//remove it from the queue
				{
					typename queue_type::iterator iter = m_queue.find(next_node);
					assert(iter != m_queue.end());
					m_queue.erase(iter);
				}
				next_node->distance_from_source() = min_node->distance_from_source() +
													distances_between_nodes[i];
				next_node->source_index() = min_node->source_index();
				next_node->previous() = min_node;
				m_queue.insert(next_node);
			}
		}
	} 

	m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->distance_from_source();
	clock_t finish = clock();
	m_time_consumed = (static_cast<double>(finish)-static_cast<double>(start))/CLOCKS_PER_SEC;
	//std::cout << std::endl;
}

template<class Node>
inline bool GeodesicAlgorithmGraphBase<Node>::check_stop_conditions()
{
	double queue_min_distance = (*m_queue.begin())->distance_from_source();

	if(queue_min_distance > m_max_propagation_distance)
	{
		return true;
	}

	for (unsigned index = 0; index < m_stop_vertices.size(); ++index)
	{
		vertex_pointer v = m_stop_vertices[index].first;
		Node& node = m_nodes[node_index(v)];
		if(queue_min_distance > node.distance_from_source() + m_stop_vertices[index].second)
		{
			return true;
		}
	}
	return false;
}

template<class Node>
inline void GeodesicAlgorithmGraphBase<Node>::trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
														 std::vector<SurfacePoint>& path) 
{
	path.clear();

	double total_path_length;
	node_pointer node = best_first_node(destination, total_path_length);

	if(total_path_length>GEODESIC_INF/2.0)		//unable to find the path
	{
		return;
	}

	path.push_back(destination);

	if(node->distance(&destination) > 1e-50)
	{
		path.push_back(node->surface_point());
	}

	while(node->previous())			//follow the path
	{
		node = node->previous();
		path.push_back(node->surface_point());
	}

	SurfacePoint& source = m_sources[node->source_index()];		//add source to the path if it is not already there
	if(node->distance(&source) > 1e-50)
	{
		path.push_back(source);
	}
}


template<class Node>
inline unsigned GeodesicAlgorithmGraphBase<Node>::best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
																 double& best_source_distance)
{
	node_pointer node = best_first_node(point, best_source_distance);
	return node ? node->source_index() : 0;
};

}		//geodesic

#endif //GEODESIC_ALGORITHM_GRAPH_BASE_010907

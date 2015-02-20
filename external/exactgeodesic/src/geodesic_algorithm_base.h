//Copyright (C) 2008 Danil Kirsanov, MIT License

#ifndef GEODESIC_ALGORITHM_BASE_122806
#define GEODESIC_ALGORITHM_BASE_122806

#include "geodesic_mesh.h"
#include "geodesic_constants_and_simple_functions.h"
#include <iostream>
#include <ctime>

namespace geodesic{

class GeodesicAlgorithmBase
{
public:
    enum AlgorithmType
    {
        EXACT,
		DIJKSTRA,
        SUBDIVISION,
		UNDEFINED_ALGORITHM
    };

	GeodesicAlgorithmBase(geodesic::Mesh* mesh = 0):
		m_type(UNDEFINED_ALGORITHM),
		m_max_propagation_distance(1e100),
		m_mesh(mesh)
	{};	

	virtual ~GeodesicAlgorithmBase(){};

	virtual void propagate(std::vector<SurfacePoint>& sources,
   						   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
						   std::vector<SurfacePoint>* stop_points = NULL) = 0; //or after ensuring that all the stop_points are covered

	virtual void trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
							std::vector<SurfacePoint>& path) = 0;

	void geodesic(SurfacePoint& source,
						  SurfacePoint& destination,
						  std::vector<SurfacePoint>& path); //lazy people can find geodesic path with one function call

	void geodesic(std::vector<SurfacePoint>& sources,
						  std::vector<SurfacePoint>& destinations,
						  std::vector<std::vector<SurfacePoint> >& paths); //lazy people can find geodesic paths with one function call

	virtual unsigned best_source(SurfacePoint& point,			//after propagation step is done, quickly find what source this point belongs to and what is the distance to this source
								 double& best_source_distance) = 0; 

	virtual void print_statistics()		//print info about timing and memory usage in the propagation step of the algorithm
	{
		std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;
	};	

	AlgorithmType type(){return m_type;};

	virtual std::string name();

	geodesic::Mesh* mesh(){return m_mesh;};
protected:

	void set_stop_conditions(std::vector<SurfacePoint>* stop_points, 
						     double stop_distance);
	double stop_distance()
	{
		return m_max_propagation_distance;
	}

	AlgorithmType m_type;					   // type of the algorithm

	typedef std::pair<vertex_pointer, double> stop_vertex_with_distace_type;
	std::vector<stop_vertex_with_distace_type> m_stop_vertices; // algorithm stops propagation after covering certain vertices
	double m_max_propagation_distance;			 // or reaching the certain distance

	geodesic::Mesh* m_mesh;

	double m_time_consumed;		//how much time does the propagation step takes
	double m_propagation_distance_stopped;		//at what distance (if any) the propagation algorithm stopped 
};

inline double length(std::vector<SurfacePoint>& path)
{
	double length = 0;
	if(!path.empty())
	{
		for(unsigned i=0; i<path.size()-1; ++i)
		{
			length += path[i].distance(&path[i+1]);
		}
	}
	return length;
}

inline void print_info_about_path(std::vector<SurfacePoint>& path)
{
	std::cout << "number of the points in the path = " << path.size()
			  << ", length of the path = " << length(path) 
			  << std::endl;
}

inline std::string GeodesicAlgorithmBase::name()
{
	switch(m_type)
	{
	case EXACT:
		return "exact";
	case DIJKSTRA:
		return "dijkstra";
	case SUBDIVISION:
		return "subdivision";
	default:
	case UNDEFINED_ALGORITHM:
		return "undefined";
	}
}

inline void GeodesicAlgorithmBase::geodesic(SurfacePoint& source,
											SurfacePoint& destination,
											std::vector<SurfacePoint>& path) //lazy people can find geodesic path with one function call
{
	std::vector<SurfacePoint> sources(1, source);
	std::vector<SurfacePoint> stop_points(1, destination);
	double const max_propagation_distance = GEODESIC_INF;

	propagate(sources, 
			  max_propagation_distance,
			  &stop_points);

	trace_back(destination, path);
}

inline void GeodesicAlgorithmBase::geodesic(std::vector<SurfacePoint>& sources,
											std::vector<SurfacePoint>& destinations,
											std::vector<std::vector<SurfacePoint> >& paths) //lazy people can find geodesic paths with one function call
{
	double const max_propagation_distance = GEODESIC_INF;

	propagate(sources, 
			  max_propagation_distance,
			  &destinations);		//we use desinations as stop points

	paths.resize(destinations.size());

	for(unsigned i=0; i<paths.size(); ++i)
	{
		trace_back(destinations[i], paths[i]);
	}
}

inline void GeodesicAlgorithmBase::set_stop_conditions(std::vector<SurfacePoint>* stop_points, 
														double stop_distance)
{
	m_max_propagation_distance = stop_distance;

	if(!stop_points)
	{
		m_stop_vertices.clear();
		return;
	}

	m_stop_vertices.resize(stop_points->size());

	std::vector<vertex_pointer> possible_vertices;
	for(unsigned i = 0; i < stop_points->size(); ++i)
	{
		SurfacePoint* point = &(*stop_points)[i];

		possible_vertices.clear();
		m_mesh->closest_vertices(point, &possible_vertices);
		
		vertex_pointer closest_vertex = NULL;
		double min_distance = 1e100;
		for(unsigned j = 0; j < possible_vertices.size(); ++j)
		{
			double distance = point->distance(possible_vertices[j]);
			if(distance < min_distance)
			{
				min_distance = distance;
				closest_vertex = possible_vertices[j];
			}
		}
		assert(closest_vertex);

		m_stop_vertices[i].first = closest_vertex;
		m_stop_vertices[i].second = min_distance;
	}
}

}//geodesic

#endif //GEODESIC_ALGORITHM_BASE_122806

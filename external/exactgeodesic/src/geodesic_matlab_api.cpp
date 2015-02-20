#define GEODESIC_DLL_IMPORT __declspec(dllexport)

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <fstream>
#include <string>

#include "geodesic_mesh.h"
#include "geodesic_algorithm_dijkstra_alternative.h"
#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"
#include "geodesic_matlab_api.h"

typedef boost::shared_ptr<geodesic::Mesh> mesh_shared_pointer;
std::vector<mesh_shared_pointer> meshes;

typedef boost::shared_ptr<geodesic::GeodesicAlgorithmBase> algorithm_shared_pointer;
std::vector<algorithm_shared_pointer> algorithms;

std::vector<geodesic::SurfacePoint> output_path;
geodesic::OutputBuffer output_buffer, output_buffer1;


std::size_t find_mesh_id(geodesic::Mesh* mesh)
{
	for(std::size_t i=0; i<meshes.size(); ++i)
	{
		if(meshes[i].get() == mesh)
		{
			return i;
		}
	}
}

long distance_and_source(long algorithm_id,		//quickly find what source this point belongs to and what is the distance to this source
											 double* destination,			
											 double* best_source_distance)
{
	geodesic::SurfacePoint point;
	geodesic::GeodesicAlgorithmBase* algorithm = algorithms[algorithm_id].get();
	std::size_t mesh_id = find_mesh_id(algorithm->mesh());
	geodesic::fill_surface_point_structure(&point, 
										   destination, 
										   algorithm->mesh());

	std::size_t best_source = algorithm->best_source(point,
													 *best_source_distance);
	return best_source;
}

long distance_and_source_for_all_vertices(long algorithm_id,	//same idea as in the previous function
															  double** distances,	//list distance/source info for all vertices of the mesh
															  long** sources)
{

	geodesic::GeodesicAlgorithmBase* algorithm = algorithms[algorithm_id].get();
	geodesic::Mesh* mesh = algorithm->mesh();

	output_buffer.allocate<double>(mesh->vertices().size());	//allocate space for distances
	*distances = output_buffer.get<double>();

	output_buffer1.allocate<long>(mesh->vertices().size());		//allocate space for sources
	*sources = output_buffer1.get<long>();


	for(std::size_t i = 0; i < mesh->vertices().size(); ++i)
	{
		geodesic::SurfacePoint point(&mesh->vertices()[i]);
		double* distace_location = *distances + i;
		(*sources)[i] = algorithm->best_source(point, *distace_location);
	}

	return mesh->vertices().size();
}


long trace_back(long algorithm_id,
									double* destination,
									double** path)
{
	geodesic::SurfacePoint point;
	geodesic::GeodesicAlgorithmBase* algorithm = algorithms[algorithm_id].get();
	geodesic::fill_surface_point_structure(&point, 
										   destination, 
										   algorithm->mesh());

	algorithm->trace_back(point,
						  output_path);

	std::size_t mesh_id = find_mesh_id(algorithm->mesh());
	output_buffer.allocate<double>(output_path.size()*5);
	for(std::size_t i=0; i<output_path.size(); ++i)
	{
		geodesic::fill_surface_point_double(&output_path[i], 
											output_buffer.get<double>() + 5*i, 
											mesh_id);
	}

	*path = output_buffer.get<double>();

	return output_path.size();
}

void propagate(long algorithm_id,
									double* source_points,	
									long num_sources,
									double* stop_points,	
									long num_stop_points,
									double max_propagation_distance)
{
	std::vector<geodesic::SurfacePoint> sources(num_sources);

	geodesic::Mesh* mesh = algorithms[algorithm_id]->mesh();
	for(std::size_t i=0; i<num_sources; ++i)
	{
		geodesic::fill_surface_point_structure(&sources[i], 
												source_points + 5*i, 
												mesh);
	}

	std::vector<geodesic::SurfacePoint> stop(num_stop_points);
	for(std::size_t i=0; i<num_stop_points; ++i)
	{
		geodesic::fill_surface_point_structure(&stop[i], 
												stop_points + 5*i, 
												mesh);
	}

	algorithms[algorithm_id]->propagate(sources, 
										max_propagation_distance,
										&stop);
}


long new_mesh(long num_points,
								  double* points,	
								  long num_triangles,
								  long* triangles, 
								  long* num_edges, 
								  double** edges)
{
	mesh_shared_pointer new_mesh = mesh_shared_pointer(new geodesic::Mesh);
	meshes.push_back(new_mesh);

	new_mesh->initialize_mesh_data(num_points,
								   points,	
								   num_triangles,
								   triangles);

	*num_edges = new_mesh->edges().size();

	output_buffer.allocate<double>(*num_edges * 4);
	*edges = output_buffer.get<double>();
	
	for(std::size_t i=0; i<*num_edges; ++i)
	{
		geodesic::Edge& edge = new_mesh->edges()[i];
		double* buffer = *edges + 4*i;

		buffer[0] = edge.adjacent_vertices()[0]->id();
		buffer[1] = edge.adjacent_vertices()[1]->id();
		buffer[2] = edge.adjacent_faces()[0]->id();
		buffer[3] = edge.adjacent_faces().size() == 2 ? 
					edge.adjacent_faces()[1]->id() : 
					-1; 
	}

	return meshes.size() - 1;
}

long new_algorithm(long mesh_id,
						               long type,
									   long subdivision)
{
	if(subdivision == 0 && type == 1)//SUBDIVISION)
	{
		type = 2;		//DIJKSTRA;
	}

	geodesic::Mesh* mesh = meshes[mesh_id].get();
	geodesic::GeodesicAlgorithmBase* algorithm;
	switch(type)
	{
		case 2: //DIJKSTRA
		{
			algorithm = new geodesic::GeodesicAlgorithmDijkstra(mesh);
			break;
		}
		case 1: //SUBDIVISION:
		{
			algorithm = new geodesic::GeodesicAlgorithmSubdivision(mesh, subdivision);
			break;
		}
		case 0://EXACT:
		default:
		{
			algorithm = new geodesic::GeodesicAlgorithmExact(mesh);
			break;
		}
	}

	algorithms.push_back(algorithm_shared_pointer(algorithm));

	return algorithms.size() - 1;
}

void delete_mesh(long id)		//delete mesh and all related algorithms
{
	assert(id < meshes.size());

	geodesic::Mesh* mesh = meshes[id].get();
	for(std::size_t i=0; i<algorithms.size(); ++i)
	{
		geodesic::GeodesicAlgorithmBase* algorithm = algorithms[i].get();
		if(algorithm && algorithm->mesh() == mesh)
		{
			algorithms[i] = algorithm_shared_pointer();
		}
	}

	meshes[id] = mesh_shared_pointer();
}

void delete_algorithm(long id)
{
	if(id < algorithms.size())
	{
		algorithms[id] = algorithm_shared_pointer();
	}
}

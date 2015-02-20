//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_EXACT_20071231
#define GEODESIC_ALGORITHM_EXACT_20071231

#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
#include <vector>
#include <cmath>
#include <assert.h>
#include <set>

namespace geodesic{

class GeodesicAlgorithmExact : public GeodesicAlgorithmBase
{
public:
    GeodesicAlgorithmExact() {}
	GeodesicAlgorithmExact(geodesic::Mesh* mesh):
	  	GeodesicAlgorithmBase(mesh),
		m_memory_allocator(mesh->edges().size(), mesh->edges().size()),
		m_edge_interval_lists(mesh->edges().size())
	{
		m_type = EXACT;

		for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
		{
			m_edge_interval_lists[i].initialize(&mesh->edges()[i]);
		}
	};	

	~GeodesicAlgorithmExact(){};

	void propagate(std::vector<SurfacePoint>& sources,
   				   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
				   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

	void trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
					std::vector<SurfacePoint>& path);

	unsigned best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
		double& best_source_distance); 

	void print_statistics();

private:
	typedef std::set<interval_pointer, Interval> IntervalQueue;

	void update_list_and_queue(list_pointer list,
							   IntervalWithStop* candidates,	//up to two candidates
							   unsigned num_candidates);

	unsigned compute_propagated_parameters(double pseudo_x, 
											double pseudo_y, 
											double d,		//parameters of the interval
											double start, 
											double end,		//start/end of the interval
											double alpha,	//corner angle
											double L,		//length of the new edge
											bool first_interval,		//if it is the first interval on the edge
											bool last_interval,
											bool turn_left,
											bool turn_right,
											IntervalWithStop* candidates);		//if it is the last interval on the edge

	void construct_propagated_intervals(bool invert, 
									  edge_pointer edge, 
									  face_pointer face,		//constructs iNew from the rest of the data
									  IntervalWithStop* candidates,
									  unsigned& num_candidates,
									  interval_pointer source_interval);

	double compute_positive_intersection(double start,
										 double pseudo_x,
										 double pseudo_y,
										 double sin_alpha,
										 double cos_alpha);		//used in construct_propagated_intervals

	unsigned intersect_intervals(interval_pointer zero, 
								    IntervalWithStop* one);			//intersecting two intervals with up to three intervals in the end

	interval_pointer best_first_interval(SurfacePoint& point, 
										double& best_total_distance, 
										double& best_interval_position,
										unsigned& best_source_index);

	bool check_stop_conditions();

	void clear()
	{
		m_memory_allocator.clear();
		m_queue.clear();
		for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
		{
			m_edge_interval_lists[i].clear();
		}
		m_propagation_distance_stopped = GEODESIC_INF;
	};

	list_pointer interval_list(edge_pointer e)
	{
		return &m_edge_interval_lists[e->id()];
	};

	void set_sources(std::vector<SurfacePoint>& sources)
	{
		m_sources.initialize(sources);
	}

	void initialize_propagation_data();		

	void list_edges_visible_from_source(MeshElementBase* p,
										std::vector<edge_pointer>& storage); //used in initialization

	long visible_from_source(SurfacePoint& point);	//used in backtracing

	void best_point_on_the_edge_set(SurfacePoint& point, 
									std::vector<edge_pointer> const& storage,
									interval_pointer& best_interval,
									double& best_total_distance,
									double& best_interval_position);

	void possible_traceback_edges(SurfacePoint& point, 
								  std::vector<edge_pointer>& storage);

	bool erase_from_queue(interval_pointer p);

	IntervalQueue m_queue;	//interval queue

	MemoryAllocator<Interval> m_memory_allocator;			//quickly allocate and deallocate intervals 
	std::vector<IntervalList> m_edge_interval_lists;		//every edge has its interval data 

	enum MapType {OLD, NEW};		//used for interval intersection
	MapType map[5];		
	double start[6];
	interval_pointer i_new[5];

	unsigned m_queue_max_size;			//used for statistics
	unsigned m_iterations;			//used for statistics

	SortedSources m_sources;
};

inline void GeodesicAlgorithmExact::best_point_on_the_edge_set(SurfacePoint& point, 
															   std::vector<edge_pointer> const& storage,
															   interval_pointer& best_interval,
															   double& best_total_distance,
															   double& best_interval_position)
{
	best_total_distance = 1e100;
	for(unsigned i=0; i<storage.size(); ++i)
	{
		edge_pointer e = storage[i];
		list_pointer list = interval_list(e);

		double offset;
		double distance;
		interval_pointer interval;

		list->find_closest_point(&point, 
								 offset, 
								 distance, 
								 interval);

		if(distance < best_total_distance)
		{
			best_interval = interval;
			best_total_distance = distance;
			best_interval_position = offset;
		}
	}
}

inline void GeodesicAlgorithmExact::possible_traceback_edges(SurfacePoint& point, 
															 std::vector<edge_pointer>& storage)
{
	storage.clear();

	if(point.type() == VERTEX)
	{
		vertex_pointer v = static_cast<vertex_pointer>(point.base_element());
		for(unsigned i=0; i<v->adjacent_faces().size(); ++i)
		{
			face_pointer f = v->adjacent_faces()[i];
			storage.push_back(f->opposite_edge(v));
		}
	}
	else if(point.type() == EDGE)
	{
		edge_pointer e = static_cast<edge_pointer>(point.base_element());
		for(unsigned i=0; i<e->adjacent_faces().size(); ++i)
		{
			face_pointer f = e->adjacent_faces()[i];

			storage.push_back(f->next_edge(e,e->v0()));
			storage.push_back(f->next_edge(e,e->v1()));
		}
	}
	else
	{
		face_pointer f = static_cast<face_pointer>(point.base_element());
		storage.push_back(f->adjacent_edges()[0]);
		storage.push_back(f->adjacent_edges()[1]);
		storage.push_back(f->adjacent_edges()[2]);
	}
}


inline long GeodesicAlgorithmExact::visible_from_source(SurfacePoint& point)	//negative if not visible
{
	assert(point.type() != UNDEFINED_POINT);

	if(point.type() == EDGE)		
	{
		edge_pointer e = static_cast<edge_pointer>(point.base_element());
		list_pointer list = interval_list(e);
		double position = std::min(point.distance(e->v0()), e->length());
		interval_pointer interval = list->covering_interval(position);
		//assert(interval);
		if(interval && interval->visible_from_source())
		{
			return (long)interval->source_index();
		}
		else
		{
			return -1;
		}
	}
	else if(point.type() == FACE)		
	{
		return -1;
	}
	else if(point.type() == VERTEX)		
	{
		vertex_pointer v = static_cast<vertex_pointer>(point.base_element());
		for(unsigned i=0; i<v->adjacent_edges().size(); ++i)
		{
			edge_pointer e = v->adjacent_edges()[i];
			list_pointer list = interval_list(e);

			double position = e->v0()->id() == v->id() ? 0.0 : e->length();
			interval_pointer interval = list->covering_interval(position);
			if(interval && interval->visible_from_source())
			{
				return (long)interval->source_index();
			}
		}

		return -1;
	}

	assert(0);
	return 0;
}

inline double GeodesicAlgorithmExact::compute_positive_intersection(double start,
																	double pseudo_x,
																	double pseudo_y,
																	double sin_alpha,
																	double cos_alpha)
{
	assert(pseudo_y < 0);

	double denominator = sin_alpha*(pseudo_x - start) - cos_alpha*pseudo_y;
	if(denominator<0.0)
	{
		return -1.0;
	}

	double numerator = -pseudo_y*start;

	if(numerator < 1e-30)
	{
		return 0.0;
	}

	if(denominator < 1e-30)
	{
		return -1.0;
	}

	return numerator/denominator;
}

inline void GeodesicAlgorithmExact::list_edges_visible_from_source(MeshElementBase* p,
																   std::vector<edge_pointer>& storage)
{
	assert(p->type() != UNDEFINED_POINT);

	if(p->type() == FACE)
	{
		face_pointer f = static_cast<face_pointer>(p);
		for(unsigned i=0; i<3; ++i)
		{
			storage.push_back(f->adjacent_edges()[i]);
		}
	}
	else if(p->type() == EDGE)
	{
		edge_pointer e = static_cast<edge_pointer>(p);
		storage.push_back(e);
	}
	else			//VERTEX
	{
		vertex_pointer v = static_cast<vertex_pointer>(p);
		for(unsigned i=0; i<v->adjacent_edges().size(); ++i)
		{
			storage.push_back(v->adjacent_edges()[i]);
		}

	}
}

inline bool GeodesicAlgorithmExact::erase_from_queue(interval_pointer p)
{
	if(p->min() < GEODESIC_INF/10.0)// && p->min >= queue->begin()->first)
	{
		assert(m_queue.count(p)<=1);			//the set is unique

		IntervalQueue::iterator it = m_queue.find(p);

		if(it != m_queue.end())
		{
			m_queue.erase(it);
			return true;
		}
	}

	return false;
}

inline unsigned GeodesicAlgorithmExact::intersect_intervals(interval_pointer zero, 
															   IntervalWithStop* one)			//intersecting two intervals with up to three intervals in the end
{
	assert(zero->edge()->id() == one->edge()->id());
	assert(zero->stop() > one->start() && zero->start() < one->stop());
	assert(one->min() < GEODESIC_INF/10.0);

	double const local_epsilon = SMALLEST_INTERVAL_RATIO*one->edge()->length(); 

	unsigned N=0;
	if(zero->min() > GEODESIC_INF/10.0)
	{
		start[0] = zero->start();
		if(zero->start() < one->start() - local_epsilon)
		{
			map[0] = OLD;
			start[1] = one->start();
			map[1] = NEW;
			N = 2;
		}
		else
		{
			map[0] = NEW;
			N = 1;
		}

		if(zero->stop() > one->stop() + local_epsilon)
		{
			map[N] = OLD;							//"zero" interval
			start[N++] = one->stop();
		}

		start[N+1] = zero->stop();
		return N;
	}

	double const local_small_epsilon = 1e-8*one->edge()->length(); 

	double D = zero->d() - one->d();
	double x0 = zero->pseudo_x();
	double x1 = one->pseudo_x();
	double R0 = x0*x0 + zero->pseudo_y()*zero->pseudo_y();
	double R1 = x1*x1 + one->pseudo_y()*one->pseudo_y();

	double inter[2];									//points of intersection
	char Ninter=0;										//number of the points of the intersection

	if(std::abs(D)<local_epsilon)					//if d1 == d0, equation is linear
	{
		double denom = x1 - x0;
		if(std::abs(denom)>local_small_epsilon)
		{
			inter[0] =  (R1 - R0)/(2.*denom);					//one solution
			Ninter = 1;
		}
	}
	else
	{
		double D2 = D*D;
		double Q = 0.5*(R1-R0-D2);
		double X = x0 - x1;

		double A = X*X - D2;
		double B = Q*X + D2*x0;
		double C = Q*Q - D2*R0;

		if (std::abs(A)<local_small_epsilon)							//if A == 0, linear equation
		{
			if(std::abs(B)>local_small_epsilon)
			{
				inter[0] =  -C/B;							//one solution
				Ninter = 1;
			}
		}
		else
		{
			double det = B*B-A*C;
			if(det>local_small_epsilon*local_small_epsilon)			//two roots
			{
				det = sqrt(det);
				if(A>0.0)								//make sure that the roots are ordered
				{
					inter[0] = (-B - det)/A;
					inter[1] = (-B + det)/A;
				}
				else
				{
					inter[0] = (-B + det)/A;
					inter[1] = (-B - det)/A;
				}
				
				if(inter[1] - inter[0] > local_small_epsilon)
				{
					Ninter = 2;
				}
				else
				{
					Ninter = 1;
				}
			}
			else if(det>=0.0)					//single root
			{
				inter[0] = -B/A;
				Ninter = 1;
			}
		}
	}
	//---------------------------find possible intervals---------------------------------------
	double left = std::max(zero->start(), one->start());		//define left and right boundaries of the intersection of the intervals
	double right = std::min(zero->stop(), one->stop());

	double good_start[4];										//points of intersection within the (left, right) limits +"left" + "right"
	good_start[0] = left;
	char Ngood_start=1;										//number of the points of the intersection	

	for(char i=0; i<Ninter; ++i)							//for all points of intersection
	{
		double x = inter[i];
		if(x > left + local_epsilon && x < right - local_epsilon)
		{
			good_start[Ngood_start++] = x;
		}
	}
	good_start[Ngood_start++] = right;

	MapType mid_map[3];
	for(char i=0; i<Ngood_start-1; ++i)
	{
		double mid = (good_start[i] + good_start[i+1])*0.5;
		mid_map[i] = zero->signal(mid) <= one->signal(mid) ? OLD : NEW;
	}

	//-----------------------------------output----------------------------------
	N = 0;
	if(zero->start() < left - local_epsilon)						//additional "zero" interval
	{
		if(mid_map[0] == OLD)				//first interval in the map is already the old one
		{
			good_start[0] = zero->start();
		}
		else
		{
			map[N] = OLD;					//"zero" interval
			start[N++] = zero->start();
		}
	}

	for(long i=0;i<Ngood_start-1;++i)							//for all intervals
	{
		MapType current_map = mid_map[i];
		if(N==0 || map[N-1] != current_map)
		{
			map[N] = current_map;				
			start[N++] = good_start[i];
		}
	}

	if(zero->stop() > one->stop() + local_epsilon)
	{
		if(N==0 || map[N-1] == NEW)
		{
			map[N] = OLD;							//"zero" interval
			start[N++] = one->stop();
		}
	}

	start[0] = zero->start();		// just to make sure that epsilons do not damage anything
	//start[N] = zero->stop();

	return N; 
}

inline void GeodesicAlgorithmExact::initialize_propagation_data()
{
	clear();

	IntervalWithStop candidate;
	std::vector<edge_pointer> edges_visible_from_source;
	for(unsigned i=0; i<m_sources.size(); ++i)		//for all edges adjacent to the starting vertex			
	{
		SurfacePoint* source = &m_sources[i];
		
		edges_visible_from_source.clear();
		list_edges_visible_from_source(source->base_element(), 
									   edges_visible_from_source);
		
		for(unsigned j=0; j<edges_visible_from_source.size(); ++j)
		{
			edge_pointer e = edges_visible_from_source[j];
			candidate.initialize(e, source, i);
            candidate.stop() = e->length();
			candidate.compute_min_distance(candidate.stop());
			candidate.direction() = Interval::FROM_SOURCE;

			update_list_and_queue(interval_list(e), &candidate, 1);
		}
	}
}

inline void GeodesicAlgorithmExact::propagate(std::vector<SurfacePoint>& sources,
   									   double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
									   std::vector<SurfacePoint>* stop_points)
{
	set_stop_conditions(stop_points, max_propagation_distance);
	set_sources(sources);
	initialize_propagation_data();

	clock_t start = clock();

	m_iterations = 0;		//for statistics
	m_queue_max_size = 0;

	IntervalWithStop candidates[2];

	while(!m_queue.empty())
	{
		m_queue_max_size = std::max((unsigned int)m_queue.size(), m_queue_max_size);

		unsigned const check_period = 10;
    	if(++m_iterations % check_period == 0)		//check if we covered all required vertices
		{
			if (check_stop_conditions())
			{
				break;
			}
		}

		interval_pointer min_interval = *m_queue.begin();
		m_queue.erase(m_queue.begin());
		edge_pointer edge = min_interval->edge();
		list_pointer list = interval_list(edge);

		assert(min_interval->d() < GEODESIC_INF);

		bool const first_interval = min_interval->start() == 0.0;
		//bool const last_interval = min_interval->stop() == edge->length();
		bool const last_interval = min_interval->next() == NULL;

		bool const turn_left = edge->v0()->saddle_or_boundary();
		bool const turn_right = edge->v1()->saddle_or_boundary();

		for(unsigned i=0; i<edge->adjacent_faces().size(); ++i)		//two possible faces to propagate
		{
			if(!edge->is_boundary())		//just in case, always propagate boundary edges
			{
				if((i == 0 && min_interval->direction() == Interval::FROM_FACE_0) ||
					(i == 1 && min_interval->direction() == Interval::FROM_FACE_1))
				{
					continue;
				}
			}

			face_pointer face = edge->adjacent_faces()[i];			//if we come from 1, go to 2
			edge_pointer next_edge = face->next_edge(edge,edge->v0());

			unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), 
																	 min_interval->pseudo_y(), 
																	 min_interval->d(),		//parameters of the interval
																	 min_interval->start(), 
																	 min_interval->stop(),		//start/end of the interval
																	 face->vertex_angle(edge->v0()),	//corner angle
																	 next_edge->length(),		//length of the new edge
																	 first_interval,		//if it is the first interval on the edge
																	 last_interval,
																	 turn_left,
																	 turn_right,
																	 candidates);		//if it is the last interval on the edge
			bool propagate_to_right = true;

			if(num_propagated)
			{
				if(candidates[num_propagated-1].stop() != next_edge->length()) 
				{
					propagate_to_right = false;
				}
				
				bool const invert = next_edge->v0()->id() != edge->v0()->id(); //if the origins coinside, do not invert intervals

				construct_propagated_intervals(invert,		//do not inverse 
											 next_edge, 
											 face,
											 candidates,
											 num_propagated,
											 min_interval);
				
				update_list_and_queue(interval_list(next_edge), 
									  candidates, 
									  num_propagated);
			}

			if(propagate_to_right)
			{
									//propogation to the right edge
				double length = edge->length();
				next_edge = face->next_edge(edge,edge->v1());

				num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), 
															 min_interval->pseudo_y(), 
															 min_interval->d(),		//parameters of the interval
															 length - min_interval->stop(), 
															 length - min_interval->start(),		//start/end of the interval
															 face->vertex_angle(edge->v1()),	//corner angle
															 next_edge->length(),		//length of the new edge
															 last_interval,		//if it is the first interval on the edge
															 first_interval,
															 turn_right,
															 turn_left,
															 candidates);		//if it is the last interval on the edge

				if(num_propagated)
				{
					bool const invert = next_edge->v0()->id() != edge->v1()->id();		//if the origins coinside, do not invert intervals

					construct_propagated_intervals(invert,		//do not inverse 
												 next_edge, 
												 face,
												 candidates,
												 num_propagated,
												 min_interval);

					update_list_and_queue(interval_list(next_edge), 
									      candidates, 
										  num_propagated);
				}
			}
		} 
	} 

	m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
	clock_t stop = clock();
	m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;

/*	for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
	{
		list_pointer list = &m_edge_interval_lists[i];
		interval_pointer p = list->first();
		assert(p->start() == 0.0);
		while(p->next())
		{
			assert(p->stop() == p->next()->start());
			assert(p->d() < GEODESIC_INF);
			p = p->next();
		}
	}*/
}


inline bool GeodesicAlgorithmExact::check_stop_conditions()
{
	double queue_distance = (*m_queue.begin())->min();
	if(queue_distance > stop_distance())
	{
		return true;
	}

	for (unsigned index = 0; index < m_stop_vertices.size(); ++index)
	{
		vertex_pointer v = m_stop_vertices[index].first;
		edge_pointer edge = v->adjacent_edges()[0];				//take any edge

		double distance = edge->v0()->id() == v->id() ? 
						  interval_list(edge)->signal(0.0) :
						  interval_list(edge)->signal(edge->length());

		if(queue_distance > distance + m_stop_vertices[index].second)
		{
			return true;
		}
	}
	return false;
}


inline void GeodesicAlgorithmExact::update_list_and_queue(list_pointer list,
												IntervalWithStop* candidates,	//up to two candidates
												unsigned num_candidates)
{
	assert(num_candidates <= 2);
	//assert(list->first() != NULL);
	edge_pointer edge = list->edge();
	double const local_epsilon = SMALLEST_INTERVAL_RATIO * edge->length(); 

	if(list->first() == NULL) 
	{
		interval_pointer* p = &list->first();
		IntervalWithStop* first;
		IntervalWithStop* second; 

		if(num_candidates == 1)
		{
			first = candidates;
			second = candidates;
			first->compute_min_distance(first->stop());
		}
		else 
		{	
			if(candidates->start() <= (candidates+1)->start())
			{
				first = candidates;
				second = candidates+1;
			}
			else 
			{
				first = candidates+1;
				second = candidates;
			}
			assert(first->stop() == second->start());

			first->compute_min_distance(first->stop());
			second->compute_min_distance(second->stop());
		}

		if(first->start() > 0.0)
		{
			*p = m_memory_allocator.allocate();
			(*p)->initialize(edge);
			p = &(*p)->next();
		}

		*p = m_memory_allocator.allocate();
		memcpy(*p,first,sizeof(Interval));
		m_queue.insert(*p);

		if(num_candidates == 2)
		{
			p = &(*p)->next();
			*p = m_memory_allocator.allocate();
			memcpy(*p,second,sizeof(Interval));
			m_queue.insert(*p);
		}

		if(second->stop() < edge->length())
		{
			p = &(*p)->next();
			*p = m_memory_allocator.allocate();
			(*p)->initialize(edge);
			(*p)->start() = second->stop();
		}
		else
		{
			(*p)->next() = NULL;
		}
		return;
	}

	bool propagate_flag;

	for(unsigned i=0; i<num_candidates; ++i)				//for all new intervals
	{
		IntervalWithStop* q = &candidates[i];
	
		interval_pointer previous = NULL;

		interval_pointer p = list->first();
		assert(p->start() == 0.0);

		while(p != NULL && p->stop() - local_epsilon < q->start())
		{
			p = p->next(); 
		}

		while(p != NULL && p->start() < q->stop() - local_epsilon)			//go through all old intervals
		{
			unsigned const N = intersect_intervals(p, q);								//interset two intervals

			if(N == 1)			
			{
				if(map[0]==OLD)	//if "p" is always better, we do not need to update anything)
				{
					if(previous)		//close previous interval and put in into the queue
					{
						previous->next() = p;
						previous->compute_min_distance(p->start());
						m_queue.insert(previous);
						previous = NULL;
					}

					p = p->next(); 
					
				}
				else if(previous)	//extend previous interval to cover everything; remove p
				{
					previous->next() = p->next(); 
					erase_from_queue(p);
					m_memory_allocator.deallocate(p);

					p = previous->next();
				}
				else				//p becomes "previous"
				{
					previous = p;
					interval_pointer next = p->next();
					erase_from_queue(p);

					memcpy(previous,q,sizeof(Interval));

					previous->start() = start[0];
					previous->next() = next;

					p = next; 
				}
				continue;
			}

			//update_flag = true;

			Interval swap(*p);							//used for swapping information
			propagate_flag = erase_from_queue(p);

			for(unsigned j=1; j<N; ++j)				//no memory is needed for the first one
			{
				i_new[j] = m_memory_allocator.allocate();	//create new intervals
			}

			if(map[0]==OLD)	//finish previous, if any
			{
				if(previous)
				{
					previous->next() = p;
					previous->compute_min_distance(previous->stop());
					m_queue.insert(previous);
					previous = NULL;
				}
				i_new[0] = p;
				p->next() = i_new[1];
				p->start() = start[0];
			}
			else if(previous)	//extend previous interval to cover everything; remove p
			{
				i_new[0] = previous;
				previous->next() = i_new[1]; 
				m_memory_allocator.deallocate(p);
				previous = NULL;
			}
			else				//p becomes "previous"
			{
				i_new[0] = p;
				memcpy(p,q,sizeof(Interval));

				p->next() = i_new[1];
				p->start() = start[0];
			}

			assert(!previous);

			for(unsigned j=1; j<N; ++j)					
			{
				interval_pointer current_interval = i_new[j];

				if(map[j] == OLD)	
				{
					memcpy(current_interval,&swap,sizeof(Interval));
				}
				else
				{
					memcpy(current_interval,q,sizeof(Interval));
				}
				
				if(j == N-1)	
				{
					current_interval->next() = swap.next();
				}
				else			
				{
					current_interval->next() = i_new[j+1];
				}

				current_interval->start() = start[j];
			}

			for(unsigned j=0; j<N; ++j)								//find "min" and add the intervals to the queue
			{
				if(j==N-1 && map[j]==NEW)
				{
					previous = i_new[j];
				}
				else
				{
					interval_pointer current_interval = i_new[j];

					current_interval->compute_min_distance(current_interval->stop());					//compute minimal distance

					if(map[j]==NEW || (map[j]==OLD && propagate_flag))
					{
						m_queue.insert(current_interval);
					}
				}
			}

			p = swap.next();
		}

		if(previous)		//close previous interval and put in into the queue
		{
			previous->compute_min_distance(previous->stop());
			m_queue.insert(previous);
			previous = NULL;
		}
	}
}

inline unsigned GeodesicAlgorithmExact::compute_propagated_parameters(double pseudo_x, 
																		double pseudo_y, 
																		double d,		//parameters of the interval
																		double begin, 
																		double end,		//start/end of the interval
																		double alpha,	//corner angle
																		double L,		//length of the new edge
																		bool first_interval,		//if it is the first interval on the edge
																		bool last_interval,
																		bool turn_left,
																		bool turn_right,
																		IntervalWithStop* candidates)		//if it is the last interval on the edge
{				
	assert(pseudo_y<=0.0);
	assert(d<GEODESIC_INF/10.0);
	assert(begin<=end);
	assert(first_interval ? (begin == 0.0) : true);

	IntervalWithStop* p = candidates;

	if(std::abs(pseudo_y) <= 1e-30)				//pseudo-source is on the edge
	{
		if(first_interval && pseudo_x <= 0.0)
		{
			p->start() = 0.0;
			p->stop() = L;
			p->d() = d - pseudo_x;
			p->pseudo_x() = 0.0;
			p->pseudo_y() = 0.0;
			return 1;
		}
		else if(last_interval && pseudo_x >= end)
		{
			p->start() = 0.0;
			p->stop() = L;
			p->d() = d + pseudo_x-end;
			p->pseudo_x() = end*cos(alpha);
			p->pseudo_y() = -end*sin(alpha);
			return 1;
		}
		else if(pseudo_x >= begin && pseudo_x <= end)
		{
			p->start() = 0.0;
			p->stop() = L;
			p->d() = d;
			p->pseudo_x() = pseudo_x*cos(alpha);
			p->pseudo_y() = -pseudo_x*sin(alpha);
			return 1;
		}
		else
		{
			return 0;
		}
	}

	double sin_alpha = sin(alpha);
	double cos_alpha = cos(alpha);

	//important: for the first_interval, this function returns zero only if the new edge is "visible" from the source
	//if the new edge can be covered only after turn_over, the value is negative (-1.0)
	double L1 = compute_positive_intersection(begin, 
											  pseudo_x, 
											  pseudo_y, 
										      sin_alpha, 
										      cos_alpha);

	if(L1 < 0 || L1 >= L)
	{
		if(first_interval && turn_left)
		{
			p->start() = 0.0;
			p->stop() = L;
			p->d() = d + sqrt(pseudo_x*pseudo_x + pseudo_y*pseudo_y);
			p->pseudo_y() = 0.0;
			p->pseudo_x() = 0.0;
			return 1;
		}
		else
		{
			return 0;
		}
	}

	double L2 = compute_positive_intersection(end, 
											  pseudo_x, 
											  pseudo_y, 
											  sin_alpha, 
											  cos_alpha);

	if(L2 < 0 || L2 >= L)
	{
		p->start() = L1;
		p->stop() = L;
		p->d() = d;
		p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
		p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;

		return 1;
	}

	p->start() = L1;
	p->stop() = L2;
	p->d() = d;
	p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
	p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;
	assert(p->pseudo_y() <= 0.0);

	if(!(last_interval && turn_right))
	{
		return 1;
	}
	else
	{
		p = candidates + 1;

		p->start() = L2;
		p->stop() = L;
		double dx = pseudo_x - end;
		p->d() = d + sqrt(dx*dx + pseudo_y*pseudo_y);
		p->pseudo_x() = end*cos_alpha;
		p->pseudo_y() = -end*sin_alpha;

		return 2;
	}
}

inline void GeodesicAlgorithmExact::construct_propagated_intervals(bool invert, 
																	edge_pointer edge, 
																	face_pointer face,		//constructs iNew from the rest of the data
																	IntervalWithStop* candidates,
																	unsigned& num_candidates,
																	interval_pointer source_interval)	//up to two candidates
{
	double edge_length = edge->length();
	double local_epsilon = SMALLEST_INTERVAL_RATIO * edge_length;
		
		//kill very small intervals in order to avoid precision problems
	if(num_candidates == 2)		
	{
		double start = std::min(candidates->start(), (candidates+1)->start());
		double stop = std::max(candidates->stop(), (candidates+1)->stop());
		if(candidates->stop()-candidates->start() < local_epsilon) // kill interval 0
		{
			*candidates = *(candidates+1);
			num_candidates = 1;
			candidates->start() = start;
			candidates->stop() = stop;
		}
		else if ((candidates+1)->stop() - (candidates+1)->start() < local_epsilon)
		{
			num_candidates = 1;
			candidates->start() = start;
			candidates->stop() = stop;
		}
	}

	IntervalWithStop* first; 
	IntervalWithStop* second; 
	if(num_candidates == 1)
	{
		first = candidates;
		second = candidates;
	}
	else
	{
		if(candidates->start() <= (candidates+1)->start())
		{
			first = candidates;
			second = candidates+1;
		}
		else 
		{
			first = candidates+1;
			second = candidates;
		}
		assert(first->stop() == second->start());
	}

	if(first->start() < local_epsilon)
	{
		first->start() = 0.0;
	}
	if(edge_length - second->stop() < local_epsilon)
	{
		second->stop() = edge_length;
	}

		//invert intervals if necessary; fill missing data and set pointers correctly
	Interval::DirectionType direction = edge->adjacent_faces()[0]->id() == face->id() ? 
										Interval::FROM_FACE_0 : 
										Interval::FROM_FACE_1;

	if(!invert)					//in this case everything is straighforward, we do not have to invert the intervals
	{
		for(unsigned i=0; i<num_candidates; ++i)
		{
			IntervalWithStop* p = candidates + i;

			p->next() = (i == num_candidates - 1) ? NULL : candidates + i + 1;
			p->edge() = edge;
			p->direction() = direction;
			p->source_index() = source_interval->source_index();

			p->min() = 0.0;					//it will be changed later on

			assert(p->start() < p->stop());
		}
	}
	else				//now we have to invert the intervals 
	{
		for(unsigned i=0; i<num_candidates; ++i)
		{
			IntervalWithStop* p = candidates + i;

			p->next() = (i == 0) ? NULL : candidates + i - 1;
			p->edge() = edge;
			p->direction() = direction;
			p->source_index() = source_interval->source_index();

			double length = edge_length;
			p->pseudo_x() = length - p->pseudo_x();

			double start = length - p->stop();
			p->stop() = length - p->start();
			p->start() = start;

			p->min() = 0;

			assert(p->start() < p->stop());
			assert(p->start() >= 0.0);
			assert(p->stop() <= edge->length());
		}
	}
}


inline unsigned GeodesicAlgorithmExact::best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
													   double& best_source_distance)
{
	double best_interval_position; 
	unsigned best_source_index;

	best_first_interval(point, 
						best_source_distance, 
						best_interval_position, 
						best_source_index);
	
	return best_source_index;
} 

inline interval_pointer GeodesicAlgorithmExact::best_first_interval(SurfacePoint& point, 
															 double& best_total_distance, 
															 double& best_interval_position, 
															 unsigned& best_source_index)
{
	assert(point.type() != UNDEFINED_POINT);

	interval_pointer best_interval = NULL;	
	best_total_distance = GEODESIC_INF;

	if(point.type() == EDGE)		
	{
		edge_pointer e = static_cast<edge_pointer>(point.base_element());
		list_pointer list = interval_list(e);

		best_interval_position = point.distance(e->v0());
		best_interval = list->covering_interval(best_interval_position);
		if(best_interval)
		{
			//assert(best_interval && best_interval->d() < GEODESIC_INF);
			best_total_distance = best_interval->signal(best_interval_position);
			best_source_index = best_interval->source_index();
		}
	}
	else if(point.type() == FACE)		
	{
		face_pointer f = static_cast<face_pointer>(point.base_element());
		for(unsigned i=0; i<3; ++i)
		{
			edge_pointer e = f->adjacent_edges()[i];
			list_pointer list = interval_list(e);

			double offset;
			double distance;
			interval_pointer interval;

			list->find_closest_point(&point, 
									 offset, 
									 distance, 
									 interval);

			if(interval && distance < best_total_distance)
			{
				best_interval = interval;
				best_total_distance = distance;
				best_interval_position = offset;
				best_source_index = interval->source_index();
			}
		}

			//check for all sources that might be located inside this face
		SortedSources::sorted_iterator_pair local_sources = m_sources.sources(f);
		for(SortedSources::sorted_iterator it=local_sources.first; it != local_sources.second; ++it)
		{
			SurfacePointWithIndex* source = *it;
			double distance = point.distance(source);
			if(distance < best_total_distance)
			{
				best_interval = NULL;
				best_total_distance = distance;
				best_interval_position = 0.0;
				best_source_index = source->index();
			}
		}
	}
	else if(point.type() == VERTEX)		
	{
		vertex_pointer v = static_cast<vertex_pointer>(point.base_element());
		for(unsigned i=0; i<v->adjacent_edges().size(); ++i)
		{
			edge_pointer e = v->adjacent_edges()[i];
			list_pointer list = interval_list(e);
			
			double position = e->v0()->id() == v->id() ? 0.0 : e->length();
			interval_pointer interval = list->covering_interval(position);
			if(interval)
			{
				double distance = interval->signal(position);

				if(distance < best_total_distance)
				{
					best_interval = interval;
					best_total_distance = distance;
					best_interval_position = position;
					best_source_index = interval->source_index();
				}
			}
		}
	}

	if(best_total_distance > m_propagation_distance_stopped)		//result is unreliable
	{
		best_total_distance = GEODESIC_INF;
		return NULL;
	}
	else
	{
		return best_interval;
	}
}

inline void GeodesicAlgorithmExact::trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
										std::vector<SurfacePoint>& path)
{					
	path.clear();
	double best_total_distance; 
	double best_interval_position;
	unsigned source_index = std::numeric_limits<unsigned>::max();
	interval_pointer best_interval = best_first_interval(destination, 
														 best_total_distance, 
														 best_interval_position, 
														 source_index);

	if(best_total_distance >= GEODESIC_INF/2.0)		//unable to find the right path
	{
		return;
	}

	path.push_back(destination);

	if(best_interval)	//if we did not hit the face source immediately  
	{
		std::vector<edge_pointer> possible_edges;
		possible_edges.reserve(10);
	
		while(visible_from_source(path.back()) < 0)		//while this point is not in the direct visibility of some source (if we are inside the FACE, we obviously hit the source)
		{
			SurfacePoint& q = path.back();

			possible_traceback_edges(q, possible_edges);

			interval_pointer interval;
			double total_distance;
			double position;

			best_point_on_the_edge_set(q, 
									   possible_edges,
									   interval,
									   total_distance,
									   position);

			//std::cout << total_distance + length(path) << std::endl;
			assert(total_distance<GEODESIC_INF);
			source_index = interval->source_index();

			edge_pointer e = interval->edge();
			double local_epsilon = SMALLEST_INTERVAL_RATIO*e->length();
			if(position < local_epsilon)
			{
				path.push_back(SurfacePoint(e->v0()));
			}
			else if(position > e->length()-local_epsilon)
			{
				path.push_back(SurfacePoint(e->v1()));
			}
			else
			{
				double normalized_position = position/e->length();
				path.push_back(SurfacePoint(e, normalized_position));
			}
		}
	}

	SurfacePoint& source = static_cast<SurfacePoint&>(m_sources[source_index]);
	if(path.back().distance(&source) > 0)
	{
		path.push_back(source);
	}
}

inline void GeodesicAlgorithmExact::print_statistics()
{
	GeodesicAlgorithmBase::print_statistics();

	unsigned interval_counter = 0;
	for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
	{
		interval_counter += m_edge_interval_lists[i].number_of_intervals();
	}
	double intervals_per_edge = (double)interval_counter/(double)m_edge_interval_lists.size();

	double memory = m_edge_interval_lists.size()*sizeof(IntervalList) + 
					interval_counter*sizeof(Interval);

	std::cout << "uses about " << memory/1e6 << "Mb of memory" <<std::endl;
	std::cout << interval_counter << " total intervals, or " 
			  << intervals_per_edge << " intervals per edge"
			  << std::endl;
	std::cout << "maximum interval queue size is " << m_queue_max_size << std::endl;
	std::cout << "number of interval propagations is " << m_iterations << std::endl;
}

}		//geodesic

#endif //GEODESIC_ALGORITHM_EXACT_20071231

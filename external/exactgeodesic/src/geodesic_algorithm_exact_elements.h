//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_EXACT_ELEMENTS_20071231
#define GEODESIC_ALGORITHM_EXACT_ELEMENTS_20071231

#include "geodesic_memory.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <cmath>
#include <assert.h>
#include <algorithm>

namespace geodesic{

class Interval;
class IntervalList;
typedef Interval* interval_pointer;
typedef IntervalList* list_pointer;

class Interval						//interval of the edge
{
public:
	
	Interval(){};
	~Interval(){};

	enum DirectionType
    {
        FROM_FACE_0,
		FROM_FACE_1,
		FROM_SOURCE,
		UNDEFINED_DIRECTION
    };

	double signal(double x)		//geodesic distance function at point x
	{
		assert(x>=0.0 && x <= m_edge->length());

		if(m_d == GEODESIC_INF)
		{
			return GEODESIC_INF;
		}
		else
		{
			double dx = x - m_pseudo_x;
			if(m_pseudo_y == 0.0)
			{
				return m_d + std::abs(dx);
			}
			else
			{
				return m_d + sqrt(dx*dx + m_pseudo_y*m_pseudo_y);
			}
		}
	}

	double max_distance(double end)	
	{
		if(m_d == GEODESIC_INF)
		{
			return GEODESIC_INF;
		}
		else
		{
			double a = std::abs(m_start - m_pseudo_x);
			double b = std::abs(end - m_pseudo_x);

			return a > b ? m_d + sqrt(a*a + m_pseudo_y*m_pseudo_y): 
						   m_d + sqrt(b*b + m_pseudo_y*m_pseudo_y);
		}
	}

	void compute_min_distance(double stop)			//compute min, given c,d theta, start, end.
	{
		assert(stop > m_start);

		if(m_d == GEODESIC_INF)
		{
			m_min = GEODESIC_INF;
		}
		else if(m_start > m_pseudo_x)
		{
			m_min = signal(m_start);
		}
		else if(stop < m_pseudo_x)
		{
			m_min = signal(stop);
		}
		else
		{
			assert(m_pseudo_y<=0);
			m_min = m_d - m_pseudo_y;
		} 
	}
			//compare two intervals in the queue
	bool operator()(interval_pointer const x, interval_pointer const y) const
	{
		if(x->min() != y->min())
		{
			return x->min() < y->min();
		}
		else if(x->start() != y->start())
		{
			return x->start() < y->start();
		}
		else
		{
			return x->edge()->id() < y->edge()->id();
		}
	}

	double stop()		//return the endpoint of the interval
	{
		return m_next ? m_next->start() : m_edge->length();
	}		

	double hypotenuse(double a, double b)
	{
		return sqrt(a*a + b*b);
	}

	void find_closest_point(double const x, 
						    double const y,
						    double& offset, 
						    double& distance);			//find the point on the interval that is closest to the point (alpha, s)

	double& start(){return m_start;}; 
	double& d(){return m_d;};
	double& pseudo_x(){return m_pseudo_x;};
	double& pseudo_y(){return m_pseudo_y;};
	double& min(){return m_min;};
	interval_pointer& next(){return m_next;};
	edge_pointer& edge(){return m_edge;};
	DirectionType& direction(){return m_direction;};
	bool visible_from_source(){return m_direction == FROM_SOURCE;};
	unsigned& source_index(){return m_source_index;};

	void initialize(edge_pointer edge, 
					SurfacePoint* point = NULL, 
					unsigned source_index = 0);

protected:
	double m_start;						//initial point of the interval on the edge
	double m_d;							//distance from the source to the pseudo-source
	double m_pseudo_x;					//coordinates of the pseudo-source in the local coordinate system
	double m_pseudo_y;					//y-coordinate should be always negative
	double m_min;						//minimum distance on the interval

	interval_pointer m_next;			//pointer to the next interval in the list	
	edge_pointer m_edge;				//edge that the interval belongs to
	unsigned m_source_index;			//the source it belongs to
	DirectionType m_direction;			//where the interval is coming from
};

struct IntervalWithStop : public Interval
{
public:
	double& stop(){return m_stop;};
protected:
	double m_stop;
};

class IntervalList						//list of the of intervals of the given edge
{
public:
	IntervalList(){m_first = NULL;};	
	~IntervalList(){};

	void clear()
	{
		m_first = NULL;
	};

	void initialize(edge_pointer e)
	{
		m_edge = e;
		m_first = NULL;
	};

	interval_pointer covering_interval(double offset)			//returns the interval that covers the offset
	{
		assert(offset >= 0.0 && offset <= m_edge->length());

		interval_pointer p = m_first; 
		while(p && p->stop() < offset)
		{
			p = p->next();
		}

		return p;// && p->start() <= offset ? p : NULL;
	};

	void find_closest_point(SurfacePoint* point, 
							double& offset, 
							double& distance, 
							interval_pointer& interval)
	{
		interval_pointer p = m_first; 
		distance = GEODESIC_INF;
		interval = NULL;

		double x,y;
		m_edge->local_coordinates(point, x, y);

		while(p)
		{
			if(p->min()<GEODESIC_INF)
			{
				double o, d;
				p->find_closest_point(x, y, o, d);
				if(d < distance)
				{
					distance = d;
					offset = o;
					interval = p;
				}
			}
			p = p->next();
		}
	};

	unsigned number_of_intervals()
	{
		interval_pointer p = m_first; 
		unsigned count = 0;
		while(p)
		{
			++count;
			p = p->next();
		}
		return count;
	}

	interval_pointer last()
	{
		interval_pointer p = m_first; 
		if(p)
		{
			while(p->next())
			{
				p = p->next();
			}
		}
		return p;
	}

	double signal(double x)
	{
		interval_pointer interval = covering_interval(x);

		return interval ? interval->signal(x) : GEODESIC_INF;
	}

	interval_pointer& first(){return m_first;};
	edge_pointer& edge(){return m_edge;};
private:
	interval_pointer m_first;			//pointer to the first member of the list
	edge_pointer m_edge;				//edge that owns this list
};

class SurfacePointWithIndex : public SurfacePoint
{
public:
	unsigned index(){return m_index;};

	void initialize(SurfacePoint& p, unsigned index)
	{
		SurfacePoint::initialize(p);
		m_index = index;
	} 

	bool operator()(SurfacePointWithIndex* x, SurfacePointWithIndex* y) const //used for sorting
	{
		assert(x->type() != UNDEFINED_POINT && y->type() !=UNDEFINED_POINT);

		if(x->type() != y->type())
		{
			return x->type() < y->type();
		}
		else 
		{
			return x->base_element()->id() < y->base_element()->id();
		}
	}

private:
	unsigned m_index;
}; 

class SortedSources : public std::vector<SurfacePointWithIndex>
{
private: 
	typedef std::vector<SurfacePointWithIndex*> sorted_vector_type;	
public:
	typedef sorted_vector_type::iterator sorted_iterator;
	typedef std::pair<sorted_iterator, sorted_iterator> sorted_iterator_pair;

	sorted_iterator_pair sources(base_pointer mesh_element)
	{
		m_search_dummy.base_element() = mesh_element;

		return equal_range(m_sorted.begin(), 
						   m_sorted.end(), 
						   &m_search_dummy, 
						   m_compare_less);
	}

	void initialize(std::vector<SurfacePoint>& sources)	//we initialize the sources by copie
	{
		resize(sources.size());
		m_sorted.resize(sources.size());
		for(unsigned i=0; i<sources.size(); ++i)
		{
			SurfacePointWithIndex& p = *(begin() + i);

			p.initialize(sources[i],i);
			m_sorted[i] = &p;
		}

		std::sort(m_sorted.begin(), m_sorted.end(), m_compare_less);
	};

	SurfacePointWithIndex& operator[](unsigned i)
	{
		assert(i < size());
		return *(begin() + i);
	}

private:
	sorted_vector_type m_sorted;
	SurfacePointWithIndex m_search_dummy;		//used as a search template
	SurfacePointWithIndex m_compare_less;			//used as a compare functor
};


inline void Interval::find_closest_point(double const rs, 
										 double const hs,
										 double& r, 
										 double& d_out)			//find the point on the interval that is closest to the point (alpha, s)
	{
		if(m_d == GEODESIC_INF)
		{
			r = GEODESIC_INF;
			d_out = GEODESIC_INF;
			return;
		}

		double hc = -m_pseudo_y;
		double rc = m_pseudo_x;
		double end = stop();

		double local_epsilon = SMALLEST_INTERVAL_RATIO*m_edge->length();
		if(std::abs(hs+hc) < local_epsilon)
		{
			if(rs<=m_start)
			{
				r = m_start;
				d_out = signal(m_start) + std::abs(rs - m_start);
			}
			else if(rs>=end)
			{
				r = end;
				d_out = signal(end) + fabs(end - rs);
			}
			else
			{
				r = rs;
				d_out = signal(rs);
			}
		}
		else
		{
			double ri = (rs*hc + hs*rc)/(hs+hc);

			if(ri<m_start)
			{
				r = m_start;
				d_out = signal(m_start) + hypotenuse(m_start - rs, hs);
			}
			else if(ri>end)
			{
				r = end;
				d_out = signal(end) + hypotenuse(end - rs, hs);
			}
			else
			{
				r = ri;
				d_out = m_d + hypotenuse(rc - rs, hc + hs);
			}
		}
	}


inline void Interval::initialize(edge_pointer edge, 
								 SurfacePoint* source,		
								 unsigned source_index)
{
	m_next = NULL;
	//m_geodesic_previous = NULL;	
	m_direction = UNDEFINED_DIRECTION;
	m_edge = edge;
	m_source_index = source_index;

	m_start = 0.0;
	//m_stop = edge->length();
	if(!source)
	{
		m_d = GEODESIC_INF;
		m_min = GEODESIC_INF;
		return;
	}
	m_d = 0;

	if(source->base_element()->type() == VERTEX)
	{
		if(source->base_element()->id() == edge->v0()->id())
		{
			m_pseudo_x = 0.0;
			m_pseudo_y = 0.0;
			m_min = 0.0;
			return;
		}
		else if(source->base_element()->id() == edge->v1()->id())
		{
			m_pseudo_x = stop();
			m_pseudo_y = 0.0;
			m_min = 0.0;
			return;
		}
	}

	edge->local_coordinates(source, m_pseudo_x, m_pseudo_y);
	m_pseudo_y = -m_pseudo_y;

	compute_min_distance(stop());
} 

}		//geodesic

#endif //GEODESIC_ALGORITHM_EXACT_ELEMENTS_20071231

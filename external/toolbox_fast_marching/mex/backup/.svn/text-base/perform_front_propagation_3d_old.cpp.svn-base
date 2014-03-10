/*=================================================================
% perform_front_propagation_3d - perform a Fast Marching front propagation.
%
%   [D,S] = perform_front_propagation_2d(W,start_points,end_points,nb_iter_max,H);
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 3 x num_start_points matrix where k is the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

#define kDead -1
#define kOpen 0
#define kFar 1

#include <math.h>
#include "mex.h"
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>

#define D_(i,j,k) D[(i)+n*(j)+n*p*(k)]
#define S_(i,j,k) S[(i)+n*(j)+n*p*(k)]
#define W_(i,j,k) W[(i)+n*(j)+n*p*(k)]
#define H_(i,j,k) H[(i)+n*(j)+n*p*(k)]
#define start_points_(i,s) start_points[(i)+3*(s)]
#define end_points_(i,s) end_points[(i)+3*(s)]

/* Global variables */
int n;			// size on X
int p;			// size on Y
int q;			// size on Z
double* D = NULL;
double* S = NULL;
double* W = NULL;
double* start_points = NULL;
double* end_points = NULL;
double* H = NULL;
int nb_iter_max = 100000;
int nb_start_points = 0;
int nb_end_points = 0;

struct point
{
	point( int ii, int jj, int kk )
	{ i = ii; j = jj; k = kk; }
	int i;
	int j;
	int k;
};
typedef std::vector<point> point_list;

inline bool end_points_reached(const int i, const int j, const int k )
{
	for( int s=0; s<nb_end_points; ++s )
	{
		if( i==((int)end_points_(0,s)) && j==((int)end_points_(1,s)) && k==((int)end_points_(2,s)) )
			return true;
	}
	return false;
}

inline bool compare_points(const point& a, const point& b)
{
	if( H==NULL )
		return D_(a.i,a.j,a.k) > D_(b.i,b.j,b.k);
	else
		return (D_(a.i,a.j,a.k)+H_(a.i,a.j,a.k)) > (D_(b.i,b.j,b.k)+H_(b.i,b.j,b.k));
}

// test the heap validity
void check_heap( int i, int j, int k, point_list& open_list )
{
	// test the heap
	for( int s=0; s<(int)open_list.size(); ++s  )
	{
		point pt = open_list[s];
		if( H==NULL )
		{
			if( D_(i,j,k)>D_(pt.i,pt.j,pt.k) )
				printf("Problem with heap.\n");
		}
		else
		{
			if( D_(i,j,k)+H_(i,j,k)>D_(pt.i,pt.j,pt.k)+H_(pt.i,pt.j,pt.k) )
				printf("Problem with heap.\n");
		}
	}
}

// select to test or not to test (debug purpose)
// #define CHECK_HEAP check_heap(i,j,k,open_list);
#define CHECK_HEAP

void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 or 5 input arguments are required."); 
	if( nlhs<1 ) 
		mexErrMsgTxt("1 or 2 output arguments are required."); 

	// first argument : weight list
	if( mxGetNumberOfDimensions(prhs[0])!= 3 )
		mexErrMsgTxt("W must be a 3D array.");
	n = mxGetDimensions(prhs[0])[0];
	p = mxGetDimensions(prhs[0])[1];
	q = mxGetDimensions(prhs[0])[2];
	W = mxGetPr(prhs[0]);
	// second argument : start_points
	start_points = mxGetPr(prhs[1]);
	int tmp = mxGetM(prhs[1]); 
	nb_start_points = mxGetN(prhs[1]);
	if( nb_start_points==0 || tmp!=3 )
		mexErrMsgTxt("start_points must be of size 3 x nb_start_poins."); 
	// third argument : end_points
	end_points = mxGetPr(prhs[2]);
	tmp = mxGetM(prhs[2]); 
	nb_end_points = mxGetN(prhs[2]);
	if( nb_end_points!=0 && tmp!=3 )
		mexErrMsgTxt("end_points must be of size 3 x nb_end_poins."); 
	// third argument : nb_iter_max
	nb_iter_max = (int) *mxGetPr(prhs[3]);
	// second argument : heuristic
	if( nrhs==5 )
	{
		H = mxGetPr(prhs[4]);
		if( mxGetNumberOfDimensions(prhs[4])!= 3 )
			mexErrMsgTxt("H must be a 3D array.");
		if( mxGetDimensions(prhs[4])[0]!=n || mxGetDimensions(prhs[4])[1]!=p || mxGetDimensions(prhs[4])[2]!=q )
			mexErrMsgTxt("H must be of size n x p x q."); 
	}
	else
		H = NULL;
	// first ouput : distance
	int dims[3] = {n,p,q};
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );
	D = mxGetPr(plhs[0]);
	// second output : state
	if( nlhs>=2 )
	{
		plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );
		S = mxGetPr(plhs[1]);
	}
	else
	{
		S = new double[n*p*q];
	}

	double h = 1.0/n;

	// initialize points
	for( int i=0; i<n; ++i )
	for( int j=0; j<p; ++j )
	for( int k=0; k<q; ++k )
	{
		D_(i,j,k) = GW_INFINITE;
		S_(i,j,k) = kFar;
	}
	point_list open_list;
	for( int s=0; s<nb_start_points; ++s )
	{
		int i = (int) start_points_(0,s);
		int j = (int) start_points_(1,s);
		int k = (int) start_points_(2,s);
		open_list.push_back( point( i,j,k ) );
		if( D_( i,j,k )==0 )
			mexErrMsgTxt("start_points should not contain duplicates."); 
		D_( i,j,k ) = 0;
		S_(i,j,k) = kOpen;
	}

	// set up the heap
	std::make_heap( open_list.begin(), open_list.end(), compare_points );

	// perform the front propagation
	int num_iter = 0;
	bool stop_iteration = GW_False;
	while( !open_list.empty() && num_iter<nb_iter_max && !stop_iteration )
	{
		num_iter++;

		// current point
		point cur_point = open_list.front();
		int i = cur_point.i;
		int j = cur_point.j;
		int k = cur_point.k;

		CHECK_HEAP;

		// remove from open list and set up state to dead
		std::pop_heap( open_list.begin(), open_list.end(), compare_points );
		open_list.pop_back();
		S_(i,j,k) = kDead;

		stop_iteration = end_points_reached(i,j,k);

		// recurse on each neighbor
		int nei_i[6] = {i+1,i,i-1,i,i,i};
		int nei_j[6] = {j,j+1,j,j-1,j,j};
		int nei_k[6] = {k,k,k,k,k-1,k+1};
		for( int s=0; s<6; ++s )
		{
			int ii = nei_i[s];
			int jj = nei_j[s];
			int kk = nei_k[s];
			if( ii>=0 && jj>=0 && ii<n && jj<p && kk>=0 && kk<q  )
			{
				double P = h/W_(ii,jj,kk);
				// compute its neighboring values
				double a1 = GW_INFINITE;
				if( ii<n-1 )
					a1 = D_(ii+1,jj,kk);
				if( ii>0 )
					a1 = GW_MIN( a1, D_(ii-1,jj,kk) );
				double a2 = GW_INFINITE;
				if( jj<p-1 )
					a2 = D_(ii,jj+1,kk);
				if( jj>0 )
					a2 = GW_MIN( a2, D_(ii,jj-1,kk) );
				double a3 = GW_INFINITE;
				if( kk<q-1 )
					a3 = D_(ii,jj,kk+1);
				if( kk>0 )
					a3 = GW_MIN( a3, D_(ii,jj,kk-1) );
				// order so that a1<a2<a3
				double tmp = 0;
				#define SWAP(a,b) tmp = a; a = b; b = tmp
				#define SWAPIF(a,b) if(a>b) { SWAP(a,b); }
				SWAPIF(a2,a3)
				SWAPIF(a1,a2)
				SWAPIF(a2,a3)
				// update its distance
				// now the equation is   (a-a1)^2+(a-a2)^2+(a-a3)^2 - P^2 = 0, with a >= a3 >= a2 >= a1.
				// =>    3*a^2 - 2*(a2+a1+a3)*a - P^2 + a1^2 + a3^2 + a2^2
				// => delta = (a2+a1+a3)^2 - 3*(a1^2 + a3^2 + a2^2 - P^2)
				double delta = (a2+a1+a3)*(a2+a1+a3) - 3*(a1*a1 + a2*a2 + a3*a3 - P*P);
				double A1 = 0;
				if( delta>=0 )
					A1 = ( a2+a1+a3 + sqrt(delta) )/3.0;
				if( A1<=a3 )
				{
					// at least a3 is too large, so we have
					// a >= a2 >= a1  and  a<a3 so the equation is 
					//		(a-a1)^2+(a-a2)^2 - P^2 = 0
					//=> 2*a^2 - 2*(a1+a2)*a + a1^2+a2^2-P^2
					// delta = (a2+a1)^2 - 2*(a1^2 + a2^2 - P^2)
					delta = (a2+a1)*(a2+a1) - 2*(a1*a1 + a2*a2 - P*P);
					A1 = 0;
					if( delta>=0 )
						A1 = 0.5 * ( a2+a1 +sqrt(delta) );
					if( A1<=a2 )
						A1 = a1 + P;
				}
				// update the value
				if( ((int) S_(ii,jj,kk)) == kDead )
				{
					// check if action has change. Should not appen for FM
					// if( A1<D_(ii,jj) )
					//	mexWarnMsgTxt("The update is not monotone");
				}
				else if( ((int) S_(ii,jj,kk)) == kOpen )
				{
					// check if action has change.
					if( A1<D_(ii,jj,kk) )
					{
						D_(ii,jj,kk) = A1;
						// TODO : see if we can skip it
						// std::make_heap( open_list.begin(), open_list.end(), compare_points );
					}
				}
				else if( ((int) S_(ii,jj,kk)) == kFar )
				{
					if( D_(ii,jj,kk)!=GW_INFINITE )
						mexWarnMsgTxt("Distance must be initialized to Inf");  
					S_(ii,jj,kk) = kOpen;
					// distance must have change.
					D_(ii,jj,kk) = A1;
					// add to open list
					open_list.push_back( point(ii,jj,kk) );
					std::push_heap( open_list.begin(), open_list.end(), compare_points );
				}
				else 
					mexErrMsgTxt("Unkwnown state."); 
			}	// end swich
		}		// end for
	}			// end while

	if( nlhs<2 )
		GW_DELETEARRAY(S);
	return;
}
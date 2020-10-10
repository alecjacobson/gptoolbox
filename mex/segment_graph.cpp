
/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef DISJOINT_SET
#define DISJOINT_SET

// disjoint-set forests using union-by-rank and path compression (sort of).

typedef struct {
  int rank;
  int p;
  int size;
} uni_elt;

class universe {
public:
  universe(int elements);
  ~universe();
  int find(int x);  
  void join(int x, int y);
  int size(int x) const { return elts[x].size; }
  int num_sets() const { return num; }

private:
  uni_elt *elts;
  int num;
};

universe::universe(int elements) {
  elts = new uni_elt[elements];
  num = elements;
  for (int i = 0; i < elements; i++) {
    elts[i].rank = 0;
    elts[i].size = 1;
    elts[i].p = i;
  }
}
  
universe::~universe() {
  delete [] elts;
}

int universe::find(int x) {
  int y = x;
  while (y != elts[y].p)
    y = elts[y].p;
  elts[x].p = y;
  return y;
}

void universe::join(int x, int y) {
  if (elts[x].rank > elts[y].rank) {
    elts[y].p = x;
    elts[x].size += elts[y].size;
  } else {
    elts[x].p = y;
    elts[y].size += elts[x].size;
    if (elts[x].rank == elts[y].rank)
      elts[y].rank++;
  }
  num--;
}

#endif

/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_GRAPH
#define SEGMENT_GRAPH

#include <algorithm>
#include <cmath>
//#include "disjoint-set.h"

// threshold function
#define THRESHOLD(size, c) (c/size)

typedef struct {
  float w;
  int a, b;
} edge;

bool operator<(const edge &a, const edge &b) {
  return a.w < b.w;
}

/*
 * Segment a graph
 *
 * Returns a disjoint-set forest representing the segmentation.
 *
 * num_vertices: number of vertices in graph.
 * num_edges: number of edges in graph
 * edges: array of edges.
 * c: constant for treshold function.
 */
universe *segment_graph(int num_vertices, int num_edges, edge *edges, 
			float c) { 
  // sort edges by weight
  std::sort(edges, edges + num_edges);

  // make a disjoint-set forest
  universe *u = new universe(num_vertices);

  // init thresholds
  float *threshold = new float[num_vertices];
  for (int i = 0; i < num_vertices; i++)
    threshold[i] = THRESHOLD(1,c);

  // for each edge, in non-decreasing weight order...
  for (int i = 0; i < num_edges; i++) {
    edge *pedge = &edges[i];
    
    // components conected by this edge
    int a = u->find(pedge->a);
    int b = u->find(pedge->b);
    if (a != b) {
      if ((pedge->w <= threshold[a]) &&
	  (pedge->w <= threshold[b])) {
	u->join(a, b);
	a = u->find(a);
	threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
      }
    }
  }

  // free up
  delete[] threshold;
  return u;
}

#endif

// mexopts = gptoolbox_mexopts('Static',false,'Debug',true);
// mex('segment_graph.cpp',mexopts{:});
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/unique.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/requires_arg.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <Eigen/Sparse>
void mexFunction(
 int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  using namespace igl::matlab;
  using namespace Eigen;
  using namespace std;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);


  mexErrMsgTxt(nrhs>0,"Too few inputs");
  mexErrMsgTxt(mxIsSparse(prhs[0]),"Matrix should be sparse");
  const mxArray * mx_data = prhs[0];
  const int m = mxGetM(mx_data);
  const int n = mxGetN(mx_data);
  mexErrMsgTxt(n == mxGetM(prhs[0]), "Matrix should be square");
  assert(mxIsSparse(mx_data));
  assert(mxGetNumberOfDimensions(mx_data) == 2);
  // TODO: It should be possible to directly load the data into the sparse
  // matrix without going through the triplets
  // Copy data immediately
  double * pr = mxGetPr(mx_data);
  mwIndex * ir = mxGetIr(mx_data);
  mwIndex * jc = mxGetJc(mx_data);
  const int num_edges = mxGetNzmax(mx_data);
  edge * edges = new edge[num_edges];
  int k = 0;
  for(int j=0; j<n;j++)
  {
    // Iterate over inside
    while(k<(int)jc[j+1])
    {
      //cout<<ir[k]<<" "<<j<<" "<<pr[k]<<endl;
      assert((int)ir[k]<m);
      assert((int)j<n);
      edges[k].a = ir[k];
      edges[k].b = j;
      edges[k].w = pr[k];
      k++;
    }
  }

  // defaults 
  int min_size = 0;
  // Threshold 
  int c = sqrt((double)n);
  {
    int i = 1;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Threshold",name) == 0)
      {
        requires_arg(i,nrhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_double(i,nrhs,prhs,name);
        c = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("MinSize",name) == 0)
      {
        requires_arg(i,nrhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_double(i,nrhs,prhs,name);
        min_size = (int)((double)*mxGetPr(prhs[++i]));
      }
      i++;
    }
  }

  universe *u = segment_graph(n, num_edges, edges, c);

  // post process small components
  for (int i = 0; i < num_edges; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }

  switch(nlhs)
  {
    case 1:
    {
      plhs[0] = mxCreateDoubleMatrix(m,1, mxREAL);
      Eigen::VectorXi C(m);
      for(int i = 0;i<m;i++)
      {
        C(i) = u->find(i);
      }
      Eigen::VectorXi uC,I,J;
      igl::unique(C,uC,I,J);
      prepare_lhs_index(J,plhs);
    }
    default: break;
  }

  delete[] edges;
  delete u;
  std::cout.rdbuf(outbuf);
}

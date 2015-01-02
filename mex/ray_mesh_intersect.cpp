#include "ray_mesh_intersect.h"

/*
 * RAY_MESH_INTERSECT  Ray/mesh intersection using the algorithm proposed by
 * Möller and Trumbore (1997).
 *
 * [flag, t, lambda_1, lambda_2, lambda_3] = ray_mesh_intersect(o, d, V, F);
 *
 * Input:
 *    o  3D vector ray origin.
 *    d  3D vector ray direction.
 *    V  #V by 3 list of vertex positions
 *    F  #F by 3 list of triangle indices
 * Output:
 *    flag  #F list of bools: (false) Reject, (true) Intersect.
 *    t  #F list of distances from the ray origin.
 *    lambda: #F by 3 list of barycentric coordinate of hits on each triangle
 * Modified from code by: 
 *    Jesus Mena
 *
 * mex ray_mesh_intersect.cpp
 * 
 */

void parse_rhs(
  int nrhs, 
  const mxArray *prhs[], 
  double *& o, 
  double *& d,
  double *& V,
  int & n,
  double *& F,
  int & m)
{
  mexErrMsgTxt(
     (mxGetM(prhs[0]) * mxGetN(prhs[0])) == 3, 
    "Ray origin must be 3d row or column vector");
  // set origin pointer
  o = mxGetPr(prhs[0]);

  mexErrMsgTxt(
     (mxGetM(prhs[1]) * mxGetN(prhs[1])) == 3, 
    "Ray direction must be 3d row or column vector");
  // set direction pointer
  d = mxGetPr(prhs[1]);

  mexErrMsgTxt(mxGetN(prhs[2]) == 3, 
    "Mesh vertex list must be #V by 3 list of 3D vertex positions");
  // set number of mesh vertices
  n = mxGetM(prhs[2]);
  // set vertex position pointers
  V = mxGetPr(prhs[2]);

  mexErrMsgTxt(mxGetN(prhs[3]) == 3, 
    "Mesh face list must be #F by 3 list of triangle indices");
  // set number of faces
  m = mxGetM(prhs[3]);
  // set face index list pointer
  F = mxGetPr(prhs[3]);
}

void prepare_lhs(
  int nlhs,
  mxArray *plhs[],
  int m,
  mxLogical *& flag,
  double *& t,
  double *& lambda)
{
  // Create pointers for output arrays
  if(nlhs > 0)
  {
    // hit or not
    plhs[0] = mxCreateLogicalMatrix(m,1);
    flag = mxGetLogicals(plhs[0]);
  }
  if(nlhs > 1)
  {
    // distances
    plhs[1] = mxCreateDoubleMatrix(m,1, mxREAL);
    t = mxGetPr(plhs[1]);
  }
  if(nlhs > 2)
  {
    // barycentric coordinates
    plhs[2] = mxCreateDoubleMatrix(m,3, mxREAL);
    lambda = mxGetPr(plhs[2]);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  mexErrMsgTxt(nrhs == 4, "The number of input arguments must be 4.");

  
  // 3D ray origin
  double * o;
  // 3D ray direction
  double * d;
  // vertex position list
  double * V;
  // face list
  double * F;
  // number of mesh vertices
  int n;
  // number of mesh faces
  int m;
  parse_rhs(nrhs,prhs,o,d,V,n,F,m);

  // #F list of bools whether face intersects ray
  mxLogical * flag = NULL;
  // #F list of distances from the ray origin
  double * t = NULL;
  // #F list of barycentric coordinate cooresponding to hits on each triangle
  double * lambda = NULL;
  prepare_lhs(nlhs,plhs,m,flag,t,lambda);
  // Just do nothing if no outputs
  if( flag == NULL)
  {
    return;
  }

  ray_mesh_intersect(o,d,V,F,n,m,flag,t,lambda);
}


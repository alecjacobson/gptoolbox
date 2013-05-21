#include "ray_mesh_intersect.h"
/*
 * BONE_VISIBLE  test whether vertices of mesh are visible to a set of 
 * points
 *
 * [flag] = bone_visible(V,F,s,d);
 *
 * Input:
 *    V  #V by 3 list of vertex positions
 *    F  #F by 3 list of triangle indices
 *    s  row vector of position of start end point of bone
 *    d  row vector of position of dest end point of bone
 * Output:
 *    flag  #V by 1 list of bools (true) visible, (false) obstructed
 *
 * mex -o bone_visible_mex bone_visible.cpp
 */

void parse_rhs(
  int nrhs, 
  const mxArray *prhs[], 
  double *& V,
  int & n,
  double *& F,
  int & m,
  double *& s,
  double *& d)
{
  mexErrMsgTxt(mxGetN(prhs[0]) == 3, 
    "Mesh vertex list must be #V by 3 list of 3D vertex positions");
  // set number of mesh vertices
  n = mxGetM(prhs[0]);
  // set vertex position pointers
  V = mxGetPr(prhs[0]);
  mexErrMsgTxt(mxGetN(prhs[1]) == 3, 
    "Mesh face list must be #F by 3 list of triangle indices");
  // set number of faces
  m = mxGetM(prhs[1]);
  // set face index list pointer
  F = mxGetPr(prhs[1]);
  mexErrMsgTxt(
     (mxGetM(prhs[2]) * mxGetN(prhs[2])) == 3, 
    "Bone start point must be 3d row or column vector");
  // set bone start point
  s = mxGetPr(prhs[2]);
  mexErrMsgTxt(
     (mxGetM(prhs[3]) * mxGetN(prhs[3])) == 3, 
    "Bone dest point must be 3d row or column vector");
  // set bone dest point
  d = mxGetPr(prhs[3]);
}

void prepare_lhs(
  int nlhs,
  mxArray *plhs[],
  int n,
  mxLogical *& flag)
{
  // Create pointers for output arrays
  if(nlhs > 0)
  {
    // hit or not
    plhs[0] = mxCreateLogicalMatrix(n,1);
    flag = mxGetLogicals(plhs[0]);
  }
}

// c = a+b
void plus(
  const double * a,
  const double * b,
  double * c)
{
  c[0] = b[0] + a[0]; c[1] = b[1] + a[1]; c[2] = b[2] + a[2]; 
}

// c = a-b
void minus(
  const double * a,
  const double * b,
  double * c)
{
  c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2]; 
}

// b =<<copy= a
void copy(
  const double * a,
  double * b)
{
  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
}

void scalar_times(
  const double s,
  const double * v,
  double * u)
{
  u[0] = s*v[0];
  u[1] = s*v[1];
  u[2] = s*v[2];
}

// returns q the closest point to p on the line segment from a to b 
void project_to_line_segment(
  const double * p,
  const double * a,
  const double * b,
  double * q)
{
  // vector from a to b
  double ab[3];
  minus(b,a,ab);
  double ab_sqr = dot(ab,ab);
  if( ab_sqr == 0 )
  {
    // a and b are the same point
    copy(a,q);
  }else
  {
    // vector from a to p
    double pa[3];
    minus(a,p,pa);
    // vector from b to p
    double t = - dot(pa,ab)/ab_sqr;
    if(t < 0)
    {
      copy(a,q);
    }else if (t>1)
    {
      copy(b,q);
    }else
    {
      scalar_times(t,ab,q);
      plus(a,q,q);
    }

  }
}



void bone_visible(
  const double * s,
  const double * d,
  const double * V,
  const double * F,
  int n,
  int m,
  mxLogical * flag)
{
  // loop over mesh vertices
  for(int j = 0; j < n; j++)
  {
    double q[3];
    q[0] = V[n*0+j];
    q[1] = V[n*1+j];
    q[2] = V[n*2+j];
    // query: can [s,d] see q?
    // we approximate this by asking if projection of q onto [s,d] can see q
    double o[3];
    project_to_line_segment(q,s,d,o);
    // ray direction
    double d[3];
    d[0] = q[0] - o[0]; d[1] = q[1] - o[1]; d[2] = q[2] - o[2];
    // normalize d
    double ud[3];
    double d_len = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    ud[0] = d[0]/d_len;
    ud[1] = d[1]/d_len;
    ud[2] = d[2]/d_len;
    // reusable corner positions
    double p1[3], p2[3], p3[3];
    // assume visible
    flag[j] = true;
    // loop over triangles
    for(int i = 0; i < m; i++)
    {
      // grab corner indices into V
      int i1 = (int)F[m*0+i]-1;
      int i2 = (int)F[m*1+i]-1;
      int i3 = (int)F[m*2+i]-1;
      if( j == i1 || j == i2 || j == i3)
      {
        // skip incident triangles
        continue;
      }
      // grab corner positions
      p1[0] = V[n*0 + i1]; p1[1] = V[n*1 + i1]; p1[2] = V[n*2 + i1];
      p2[0] = V[n*0 + i2]; p2[1] = V[n*1 + i2]; p2[2] = V[n*2 + i2];
      p3[0] = V[n*0 + i3]; p3[1] = V[n*1 + i3]; p3[2] = V[n*2 + i3];
      // determine if this triangle intersects line segment o --> q
      double t,lambda0,lambda1,lambda2;
      bool hit = 
        ray_triangle_intersect(o,ud,p1,p2,p3,t,lambda0,lambda1,lambda2);
      // if there was a hit, check that it wasn't past q
      if(hit)
      {
        // was hit triangle closer than q
        if(t<d_len)
        {
          // then q is not visible 
          flag[j] = false;
          break;
        }
      }
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  mexErrMsgTxt(nrhs == 4, "The number of input arguments must be 4.");
  // 3D query bone start and dest point
  double * s;
  double * d;
  // vertex position list
  double * V;
  // face list
  double * F;
  // number of mesh vertices
  int n;
  // number of mesh faces
  int m;
  parse_rhs(nrhs,prhs,V,n,F,m,s,d);
  // #F list of bools whether face intersects ray
  mxLogical * flag = NULL;
  prepare_lhs(nlhs,plhs,n,flag);
  // Just do nothing if no outputs
  if( flag == NULL)
  {
    return;
  }
  bone_visible(s,d,V,F,n,m,flag);
}

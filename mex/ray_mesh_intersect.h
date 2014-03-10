#include "mex.h"
#include <cmath>

// Overload mexErrMsgTxt to check an assertion then print text only if
// assertion fails
void mexErrMsgTxt(bool assertion, const char * text)
{
  if(!assertion)
  {
    mexErrMsgTxt(text);
  }
}
void mexPrintVector(const double * a, bool newline = true)
{
  mexPrintf("%g %g %g",a[0],a[1],a[2]);
  if(newline)
  {
    mexPrintf("\n");
  }
}

// return dot(a,b)
double dot(const double * a, const double * b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// c = cross(a,b)
void cross(const double * a, const double * b, double *c)
{
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}


#define EPSILON 2e-16
mxLogical ray_triangle_intersect(
  const double o[3],
  const double ud[3],
  const double p1[3],
  const double p2[3],
  const double p3[3],
  double & ti,
  double & lambda0,
  double & lambda1,
  double & lambda2)
{

    // compute edge vectors
    double e31[3],e32[3];
    e31[0] = p1[0] - p3[0]; e31[1] = p1[1] - p3[1]; e31[2] = p1[2] - p3[2];
    e32[0] = p2[0] - p3[0]; e32[1] = p2[1] - p3[1]; e32[2] = p2[2] - p3[2];
    // compute d cross e32
    double q[3];
    cross(ud,e32,q);
    // compute e31 dot q
    double a = dot(e31,q);
    // direction parallel to plane
    if(a > -EPSILON && a < EPSILON)
    {
      return false;
    }
    double f = 1/a;
    double s[3];
    s[0] = o[0]-p3[0];
    s[1] = o[1]-p3[1];
    s[2] = o[2]-p3[2];
    lambda0 = f*dot(s,q);
    // first BC says the intersection is outside the triangle
    if( lambda0 < 0)
    {
      return false;
    }
    double r[3];
    cross(s,e31,r);
    lambda1 = f*dot(r,ud);
    // third BC is implicitly defined by first two since they always sum to one
    lambda2 = 1-(lambda0 + lambda1);
    // second or implicitly third BC says the intersection is outside the
    // triangle
    if( lambda1 < 0 || lambda2 < 0) 
    {
      return false;
    }
    // Now we know that the line of the ray intersects the triangle, we just
    // need to check the direction
    ti = f*dot(e32,r);
    // wrong side of ray line
    if(ti < 0)
    {
      return false;
    }
    // Passed all tests: we know it's a hit
    return true;
}


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
 */

void ray_mesh_intersect(
  const double * o,
  const double * d,
  const double * V,
  const double * F,
  int n,
  int m,
  mxLogical * flag,
  double * t,
  double * lambda)
{
  // normalize d
  double ud[3];
  double d_len = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
  ud[0] = d[0]/d_len;
  ud[1] = d[1]/d_len;
  ud[2] = d[2]/d_len;
  // reusable corner positions
  double p1[3], p2[3], p3[3];
  // loop over triangles
  for(size_t i = 0; i < m; i++)
  {
    // grab corner indices into V
    int i1 = (int)F[m*0+i]-1;
    int i2 = (int)F[m*1+i]-1;
    int i3 = (int)F[m*2+i]-1;
    // grab corner positions
    p1[0] = V[n*0 + i1]; p1[1] = V[n*1 + i1]; p1[2] = V[n*2 + i1];
    p2[0] = V[n*0 + i2]; p2[1] = V[n*1 + i2]; p2[2] = V[n*2 + i2];
    p3[0] = V[n*0 + i3]; p3[1] = V[n*1 + i3]; p3[2] = V[n*2 + i3];
    double ti, lambda0, lambda1, lambda2;
    flag[i] = 
      ray_triangle_intersect(o,ud,p1,p2,p3,ti,lambda0,lambda1,lambda2);

    // only set other outputs if hit
    if(NULL != lambda)
    {
      lambda[m*0 + i] = lambda0;
      lambda[m*1 + i] = lambda1;
      lambda[m*2 + i] = lambda2;
    }
    if(NULL != t)
    {
      t[i] = ti;
    }
  }
}

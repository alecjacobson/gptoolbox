#include "mex.h"
#include <cmath>
#include <iostream>

// SOLID_ANGLE Compute the solid angle of a triangle (tetrahedron) described
// by points (vectors) V
//
// Inputs:
//  V  n by 3 list of vertex positions
//  n  number of mesh vertices
//  F  #F by 3 list of triangle indices
//  m  number of faces
//  O  3 by no list of origin positions
//  no  number of origins
// Outputs:
//  S  no by f list of solid angles
//
void solid_angle_mex_2(
  const double * V,
  const int n,
  const double * F,
  const int m,
  const double * O,
  const int no,
  double * S);

// Overload mexErrMsgTxt to check an assertion then print text only if
// assertion fails
static void mexErrMsgTxt(bool assertion, const char * text)
{
  if(!assertion)
  {
    mexErrMsgTxt(text);
  }
}
static void mexPrintVector(const double * a, bool newline = true)
{
  mexPrintf("%g %g %g",a[0],a[1],a[2]);
  if(newline)
  {
    mexPrintf("\n");
  }
}

// Parse right hand side arguments
//
// Inputs:
//   nrhs  number of right hand side arguments
//   prhs  pointer to right hand side arguments
// Outputs:
//   V  n by 3 list of mesh vertex positions
//   n  number of vertices
//   F  m by 3 list of mesh face indices
//   m  number of faces
//   O  o by 3 list of origins
//   no  number of origins
void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  double *& V,
  int & n,
  double *& F,
  int & m,
  double *& O,
  int & no)
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

  mexErrMsgTxt(mxGetN(prhs[2]) == 3, 
    "Origin list must be #O by 3 list of 3D positions");
  // set number of faces
  no = mxGetM(prhs[2]);
  // set face index list pointer
  O = mxGetPr(prhs[2]);
}

// Prepare left hand side
//
// Inputs:
//   nlhs  number of left hand side arguments
//   plhs  pointer to left hand side arguments
//   m  number of faces
//   no  number of origins
// Outputs:
//   S  o by m  list of solid angles
void prepare_lhs(
  const int nlhs,
  mxArray *plhs[],
  const int m,
  const int no,
  double *& S)
{
  // Create pointers for output arrays
  if(nlhs > 0)
  {
    // solid angles
    plhs[0] = mxCreateDoubleMatrix(no,m, mxREAL);
    S = mxGetPr(plhs[0]);
  }
}

void solid_angle_mex2(
  const double * V,
  const int n,
  const double * F,
  const int m,
  const double * O,
  const int no,
  double * S)
{
  // loop over faces
#pragma omp parallel for
  for(int f = 0;f<m;f++)
  {
    //if(f==0)
    //{
    //  mexPrintf("Num threads: %d\n",omp_get_num_threads());
    //}
    // Gather corners 
    double C[3][3];
    // loop around triangle
    for(int t=0;t<3;t++)
    {
      // loop over dimensions
      for(int d = 0;d<3;d++)
      {
        // Indices are offset by 1
        int Ff = F[m*t + f]-1;
        C[t][d] = V[d*n + Ff];
      }
    }
    // loop over origins
    for(int o = 0;o<no;o++)
    {
      // Gather vectors to corners
      double v[3][3];
      double vl[3];
      // loop around triangle
      for(int t=0;t<3;t++)
      {
        vl[t] = 0;
        // loop over dimensions
        for(int d = 0;d<3;d++)
        {
          v[t][d] = C[t][d] - O[d*no + o];
          // compute edge length contribution
          vl[t] += v[t][d]*v[t][d];
        }
        //mexPrintVector(v[t],true);
        // finish edge length computation
        // Matlab crashes on NaN
        if(vl[t]!=0)
        {
          vl[t] /= sqrt(vl[t]);
        }
      }
      //printf("\n");
      
      // Compute determinant
      double detf = 
        v[0][0]*v[1][1]*v[2][2]+
        v[1][0]*v[2][1]*v[0][2]+
        v[2][0]*v[0][1]*v[1][2]-
        v[2][0]*v[1][1]*v[0][2]-
        v[1][0]*v[0][1]*v[2][2]-
        v[0][0]*v[2][1]*v[1][2];

      // Compute pairwise dotproducts
      double dp[3];

      dp[0] = v[1][0]*v[2][0];
      dp[0] += v[1][1]*v[2][1];
      dp[0] += v[1][2]*v[2][2];

      dp[1] = v[2][0]*v[0][0];
      dp[1] += v[2][1]*v[0][1];
      dp[1] += v[2][2]*v[0][2];

      dp[2] = v[0][0]*v[1][0];
      dp[2] += v[0][1]*v[1][1];
      dp[2] += v[0][2]*v[1][2];

      S[f*no+o] = 2*atan2(detf,
        vl[0]*vl[1]*vl[2] + 
        dp[0]*vl[0] +
        dp[1]*vl[1] +
        dp[2]*vl[2]);
      //S[o*m +f] = detf;
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  using namespace std;
  mexErrMsgTxt(nrhs == 3, "The number of input arguments must be 3.");
  
  // vertex position list
  double * V;
  // face list
  double * F;
  // origin list
  double * O;
  // number of mesh vertices
  int n;
  // number of mesh faces
  int m;
  // number of origins
  int no;
  parse_rhs(nrhs,prhs,V,n,F,m,O,no);

  //// Set up openmp
  //int nProcessors=omp_get_max_threads();
  ////std::cout<<nProcessors<<std::endl;
  //omp_set_num_threads(min(nProcessors,2));
  //mexPrintf("Using %d out of %d threads\n",min(nProcessors,2),nProcessors);

  // List of solid angles
  double * S = NULL;
  if(nlhs == 0)
  {
    return;
  }
  prepare_lhs(nlhs,plhs,m,no,S);
  solid_angle_mex2(V,n,F,m,O,no,S);
}

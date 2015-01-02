#ifdef MEX
#include "parse_rhs.h"
#include <igl/C_STR.h>
#include <cstring>
#include <iostream>

// Overload mexErrMsgTxt to check an assertion then print text only if
// assertion fails
static void mexErrMsgTxt(bool assertion, const char * text)
{
  if(!assertion)
  {
    mexErrMsgTxt(text);
  }
}

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  int & dim,
  double *& V,
  int & n,
  double *& F,
  int & m,
  double *& O,
  int & no,
  winding_number_params & params)
{
  using namespace std;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3 || dim == 2,
    "Mesh vertex list must be #V by 2 or 3 list of vertex positions");

  mexErrMsgTxt(dim == mxGetN(prhs[1]),
    "Mesh \"face\" simplex size must equal dimension");

  // set number of mesh vertices
  n = mxGetM(prhs[0]);
  // set vertex position pointers
  V = mxGetPr(prhs[0]);

  // set number of faces
  m = mxGetM(prhs[1]);
  // set face index list pointer
  F = mxGetPr(prhs[1]);

  mexErrMsgTxt(mxGetN(prhs[2]) == dim, 
    "Origin list dimension must match vertex list dimension");
  // set number of faces
  no = mxGetM(prhs[2]);
  // set face index list pointer
  O = mxGetPr(prhs[2]);

  // Default values
  params.hierarchical = true;
  params.ray_cast = false;
  params.ray_parity = false;
  params.num_rays = 32;
  params.twod_rays = false;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),
        "Parameter names should be char strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Hierarchical",name) == 0)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
        i++;
        mexErrMsgTxt(mxIsLogical(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Logical argument"));
        mxLogical * v = (mxLogical *)mxGetData(prhs[i]);
        params.hierarchical = *v;
        mexErrMsgTxt(!params.hierarchical || dim==3,
          "Hierarchical only supported for dim = 3");
      }else if(strcmp("RayCast",name) == 0)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
        i++;
        mexErrMsgTxt(mxIsLogical(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Logical argument"));
        mxLogical * v = (mxLogical *)mxGetData(prhs[i]);
        params.ray_cast = *v;
        mexErrMsgTxt(!params.ray_cast || dim==3,
          "RayCast only supported for dim = 3");
      }else if(strcmp("RayParity",name) == 0)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
        i++;
        mexErrMsgTxt(mxIsLogical(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Logical argument"));
        mxLogical * v = (mxLogical *)mxGetData(prhs[i]);
        params.ray_parity = *v;
        mexErrMsgTxt(!params.ray_parity || dim==3,
          "RayParity only supported for dim = 3");
      }else if(strcmp("TwoDRays",name) == 0)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
        i++;
        mexErrMsgTxt(mxIsLogical(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Logical argument"));
        mxLogical * v = (mxLogical *)mxGetData(prhs[i]);
        params.twod_rays= *v;
      }else if(strcmp("NumRays",name) == 0)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
        i++;
        mexErrMsgTxt(mxIsDouble(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Double argument"));
        mexErrMsgTxt(mxGetN(prhs[i])==1, "NumRays should be size 1");
        mexErrMsgTxt(mxGetM(prhs[i])==1, "NumRays should be size 1");
        params.num_rays = (int) *mxGetPr(prhs[i]);

      }else
      {
        mexErrMsgTxt(false,
          C_STR("Unsupported parameter: "<<name));
      }
      i++;
    }
  }
}

#endif

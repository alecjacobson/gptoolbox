#ifdef MEX
// mex -v -largeArrayDims -DMEX ...
//   -I/usr/local/igl/libigl/include ...
//   -I/opt/local/include/eigen3 ...
//   -I/opt/local/include ...
//   -L/opt/local/lib ...
//   -lCGAL -lCGAL_Core -lgmp -lmpfr ...
//   -lboost_thread-mt -lboost_system-mt ...
//   signed_distance_isosurface.cpp

#include <igl/cgal/signed_distance_isosurface.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>

#include <cstring>

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & IV,
  Eigen::MatrixXi & IF,
  double & level,
  double & angle_bound,
  double & radius_bound,
  double & distance_bound,
  igl::SignedDistanceType & type)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 2, "The number of input arguments must be >=2.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3,
    "Mesh vertex list must be #V by 3 list of vertex positions");

  parse_rhs_double(prhs,IV);
  parse_rhs_index(prhs+1,IF);

  // defaults
  level = 0.0;
  angle_bound = 28.0;
  double bbd = 1.0;
  // Radius and distance in terms of fraction of bbd
  if(IV.size() > 0)
  {
    bbd = (IV.colwise().maxCoeff()-IV.colwise().minCoeff()).norm();
  }
  radius_bound = 0.02*bbd;
  distance_bound = 0.02*bbd;
  type = SIGNED_DISTANCE_TYPE_DEFAULT;

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      const auto requires_arg = 
        [](const int i, const int nrhs, const char * name)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
      };
      const auto validate_double = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsDouble(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Double argument"));
      };
      const auto validate_scalar = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxGetN(prhs[i])==1 && mxGetM(prhs[i])==1,
          C_STR("Parameter '"<<name<<"' requires scalar argument"));
      };
      const auto validate_char = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsChar(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires char argument"));
      };
      if(strcmp("Level",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        validate_scalar(i,prhs,name);
        level = (double)*mxGetPr(prhs[i]);
      }else if(strcmp("AngleBound",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        validate_scalar(i,prhs,name);
        angle_bound = (double)*mxGetPr(prhs[i]);
      }else if(strcmp("RadiusBound",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        validate_scalar(i,prhs,name);
        radius_bound = ((double)*mxGetPr(prhs[i])) * bbd;
      }else if(strcmp("DistanceBound",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        validate_scalar(i,prhs,name);
        distance_bound = ((double)*mxGetPr(prhs[i])) * bbd;
      }else if(strcmp("SignedDistanceType",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_char(i,prhs,name);
        const char * type_name = mxArrayToString(prhs[i]);
        if(strcmp("pseudonormal",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
        }else if(strcmp("winding_number",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
        }else if(strcmp("default",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_DEFAULT;
        }else if(strcmp("unsigned",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_UNSIGNED;
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown SignedDistanceType: "<<type_name));
        }
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  if(type != igl::SIGNED_DISTANCE_TYPE_UNSIGNED)
  {
    mexErrMsgTxt(dim == mxGetN(prhs[1]),
      "Mesh \"face\" simplex size must equal dimension");
  }
}

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::cgal;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd IV,V;
  MatrixXi IF,F;
  double level,angle_bound,radius_bound,distance_bound;
  SignedDistanceType type;
  parse_rhs(nrhs,prhs,IV,IF,level,angle_bound,radius_bound,distance_bound,type);
  signed_distance_isosurface(
    IV,IF,level,angle_bound,radius_bound,distance_bound,type,V,F);
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 2:
    {
      prepare_lhs_index(F,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(V,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#endif

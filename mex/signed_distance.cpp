#ifdef MEX
// mex -v -largeArrayDims -DMEX ...
//   -I/usr/local/igl/libigl/include ...
//   -I/opt/local/include/eigen3 ...
//   -I/opt/local/include ...
//   -L/opt/local/lib ...
//   -lCGAL -lCGAL_Core -lgmp -lmpfr ...
//   -lboost_thread-mt -lboost_system-mt ...
//   -o signed_distance signed_distance.cpp

#include <igl/cgal/signed_distance.h>
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
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  igl::SignedDistanceType & type)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3,
    "Mesh vertex list must be #P by 3 list of vertex positions");
  mexErrMsgTxt(mxGetN(prhs[1]) == 3,
    "Mesh vertex list must be #V by 3 list of vertex positions");

  mexErrMsgTxt(dim == mxGetN(prhs[2]),
    "Mesh \"face\" simplex size must equal dimension");

  parse_rhs_double(prhs,P);
  parse_rhs_double(prhs+1,V);
  parse_rhs_index(prhs+2,F);

  type = SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      const auto requires_arg = 
        [](const int i, const int nrhs, const char * name)
      {
        // Windows doesn't find igl::mexErrMsgTxt overload
        igl::mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
      };
      const auto validate_char = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        // Windows doesn't find igl::mexErrMsgTxt overload
        igl::mexErrMsgTxt(mxIsChar(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires char argument"));
      };
      if(strcmp("SignedDistanceType",name) == 0)
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
}

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;

  igl::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd P,V,C,N;
  MatrixXi F;
  VectorXi I;
  VectorXd S;
  SignedDistanceType type;
  parse_rhs(nrhs,prhs,P,V,F,type);
  signed_distance<CGAL::Simple_cartesian<double> >(P,V,F,type,S,I,C,N);
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 4:
    {
      prepare_lhs_double(N,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_double(C,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(I,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(S,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#endif


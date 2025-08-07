#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab_format.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/triangle/refine.h>
#include <Eigen/Core>

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd V;
  Eigen::MatrixXi E,F;
  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,E);
  parse_rhs_index(prhs+2,F);

  std::string flags = "";
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Flags",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        flags = mxArrayToString(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }


  Eigen::MatrixXd TV;
  Eigen::MatrixXi TF;
  try
  {
    igl::triangle::refine(V,E,F,flags,TV,TF);
  }catch(std::runtime_error e)
  {
    ::mexErrMsgTxt((std::string("triangle failed: ")+e.what()).c_str());
  }
  switch(nlhs)
  {
    case 2:
      prepare_lhs_index(TF,plhs+1);
    case 1:
      prepare_lhs_double(TV,plhs+0);
    default:break;
  }
}

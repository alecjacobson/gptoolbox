#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/triangle/triangulate.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

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
  Eigen::MatrixXi E;
  Eigen::MatrixXd H;
  Eigen::VectorXi VM,EM;
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,E);

  std::string flags = "";
  bool verbose = false;

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Flags",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        flags = mxArrayToString(prhs[++i]);
      }else if(strcmp("EdgeMarkers",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_index(prhs+(++i),EM);
      }else if(strcmp("Holes",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_double(prhs+(++i),H);
      }else if(strcmp("VertexMarkers",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_index(prhs+(++i),VM);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }


  Eigen::MatrixXd TV;
  Eigen::MatrixXi TF;
  Eigen::VectorXi TVM,TEM;
  try
  {
    igl::triangle::triangulate(V,E,H,VM,EM,flags,TV,TF,TVM,TEM);
  }catch(std::runtime_error e)
  {
    ::mexErrMsgTxt((std::string("triangle failed: ")+e.what()).c_str());
  }

  switch(nlhs)
  {
    case 4:
      prepare_lhs_index(TEM,plhs+3);
    case 3:
      prepare_lhs_index(TVM,plhs+2);
    case 2:
      prepare_lhs_index(TF,plhs+1);
    case 1:
      prepare_lhs_double(TV,plhs+0);
    default:break;
  }
}



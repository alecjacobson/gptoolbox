#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
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
  using namespace igl::copyleft::tetgen;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs+0,SV);
  parse_rhs_index(prhs+1,SF);

  Eigen::MatrixXd SH(0,3);
  std::string flags = "-q2";
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
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  Eigen::MatrixXd SR(0,5);
  size_t nr = -1;
  Eigen::MatrixXd TV,TR;
  Eigen::MatrixXi TT,TF,TN,PT,FT;
  int ret = tetrahedralize(SV,SF,SH,SR,flags,TV,TT,TF,TR,TN,PT,FT,nr);
  //int ret = tetrahedralize(SV,SF,flags,TV,TT,TF);
  mexErrMsgTxt((ret==0) || (ret==2),"Tetgen failed.");

  switch(nlhs)
  {
    case 7:
      prepare_lhs_index(FT,plhs+7);
    case 6:
      prepare_lhs_index(PT,plhs+5);
    case 5:
      prepare_lhs_index(TN,plhs+4);
    case 4:
      prepare_lhs_double(TR,plhs+3);
    case 3:
      prepare_lhs_index(TF,plhs+2);
    case 2:
      prepare_lhs_index(TT,plhs+1);
    case 1:
      prepare_lhs_double(TV,plhs+0);
    default:break;
  }
}


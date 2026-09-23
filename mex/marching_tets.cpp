#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/marching_tets.h>
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

  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT;
  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  parse_rhs_double(prhs+0,TV);
  parse_rhs_index(prhs+1,TT);
  Eigen::VectorXd S;
  parse_rhs_double(prhs+2,S);
  mexErrMsgTxt(TV.cols()==3,"TV should be #TV by 3");
  mexErrMsgTxt(TT.cols()==4,"TT should be #TT by 4");
  mexErrMsgTxt(S.rows()==TV.rows(),"S should be #TV by 1");

  double isovalue = 0;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("IsoValue",name) == 0)
      {
        validate_arg_scalar(i,nrhs,prhs,name);
        isovalue = *mxGetPr(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi J;
  Eigen::SparseMatrix<double> BC;

  igl::marching_tets(TV,TT,S,isovalue,SV,SF,J,BC);

  switch(nlhs)
  {
    case 4:
      prepare_lhs_double(BC,plhs+3);
    case 3:
      prepare_lhs_index(J,plhs+2);
    case 2:
      prepare_lhs_index(SF,plhs+1);
    case 1:
      prepare_lhs_double(SV,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

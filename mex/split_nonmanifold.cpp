#include <igl/split_nonmanifold.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace igl;
  using namespace igl::matlab;
  Eigen::MatrixXi F,SF;
  Eigen::VectorXi I;
    
  
  mexErrMsgTxt(nrhs>=1,"nrhs should be >= 1");
  parse_rhs_index(prhs,F);
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  
  igl::split_nonmanifold(F,SF,I);

  switch(nlhs)
  {
      case 2:
          prepare_lhs_index(I,plhs+1);
      case 1:
          prepare_lhs_index(SF,plhs);
      default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}


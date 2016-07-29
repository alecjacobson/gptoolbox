#include <igl/decimate.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#ifdef MEX
#  include <mex.h>
#  include <igl/C_STR.h>
#  include <igl/matlab/mexErrMsgTxt.h>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd V,W;
  MatrixXi F,G;
  VectorXi J;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,F);
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  mexErrMsgTxt(
    mxIsDouble(prhs[2]) && mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1,
    "fraction to decimate should be scalar");
  double ratio = * mxGetPr(prhs[2]);
  mexErrMsgTxt((ratio>0 && ratio<1) || (ratio>0 && ratio<F.rows()) ,"Ratio should be in (0,1) or [1,#F)");
  const size_t max_m = ratio<1 ? ratio*F.rows() : ratio;
  decimate(V,F,max_m,W,G,J);

  switch(nlhs)
  {
    case 3:
      prepare_lhs_index(J,plhs+2);
    case 2:
      prepare_lhs_index(G,plhs+1);
    case 1:
      prepare_lhs_double(W,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/blue_noise.h>
#include <Eigen/Core>

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
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexErrMsgTxt(nrhs==3,"nrhs should be == 3");
  Eigen::MatrixXd V;
  parse_rhs_double(prhs+0,V);
  Eigen::MatrixXi F;
  parse_rhs_index(prhs+1,F);
  const double r = (double)*mxGetPr(prhs[2]);
  Eigen::VectorXi FI;
  Eigen::MatrixXd B,P;
  igl::blue_noise(V,F,r,B,FI,P);

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_double(B,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(FI,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(P,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}


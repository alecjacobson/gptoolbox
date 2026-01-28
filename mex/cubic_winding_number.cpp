#include <Eigen/Core>
#include <iostream>
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/predicates/cubic_winding_number.h>


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

  mexErrMsgTxt(nrhs==2,"nrhs should be == 2");
  Eigen::MatrixXd C;
  parse_rhs_double(prhs+0,C);
  mexErrMsgTxt(C.rows()==4 && C.cols()==2,"C should be 4 x 2.");
  Eigen::RowVectorXd q;
  parse_rhs_double(prhs+1,q);
  // Q should be n x 2
  mexErrMsgTxt(q.cols()==2,"Q should have 2 columns.");

  const double wn = igl::predicates::cubic_winding_number(C,q);

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 1:
    {
      prepare_lhs_double(wn,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}



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
#include <igl/cycodebase/point_cubic_squared_distance.h>


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
  Eigen::MatrixXd Q;
  parse_rhs_double(prhs+0,Q);
  // Q should be n x 2
  Eigen::MatrixXd C;
  parse_rhs_double(prhs+1,C);
  // C should be 4 x 2
  mexErrMsgTxt(C.rows()==4,"C should have 4 rows");
  mexErrMsgTxt(C.cols() == Q.cols(), "C and Q should have the same number of columns");

  Eigen::VectorXd sqrD;
  Eigen::VectorXd S;
  Eigen::MatrixXd K;

  igl::cycodebase::point_cubic_squared_distance(Q, C, sqrD, S, K);

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_double(K,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_double(S,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(sqrD,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}


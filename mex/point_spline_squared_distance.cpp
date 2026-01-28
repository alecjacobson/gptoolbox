#include <Eigen/Core>
#include <iostream>
#include <mex.h>
#include "parse_spline_eytzinger_aabb.h"
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/cycodebase/point_spline_squared_distance.h>


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

  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  Eigen::MatrixXd Q;
  parse_rhs_double(prhs+0,Q);

  Eigen::MatrixXd P;
  parse_rhs_double(prhs+1,P);
  mexErrMsgTxt(P.cols() == Q.cols(), "P and Q should have the same number of columns");

  Eigen::MatrixXi C;
  parse_rhs_index(prhs+2,C);
  mexErrMsgTxt(C.cols()==4,"C should be #C x 4.");
  
  Eigen::MatrixXd B1,B2;
  Eigen::VectorXi leaf;
  bool acceleration_provided = false;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Acceleration",name) == 0)
      {
        // i will be incremented on success
        acceleration_provided = parse_spline_eytzinger_aabb(i,nrhs,prhs,name,B1,B2,leaf);
      }else
      {
        mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
      }
      i++;
    }
  }

  if(!acceleration_provided)
  {
    igl::cycodebase::spline_eytzinger_aabb(P, C, B1, B2,leaf);
  }

  Eigen::VectorXd sqrD;
  Eigen::VectorXd S;
  Eigen::MatrixXd K;
  Eigen::VectorXi I;

  igl::cycodebase::point_spline_squared_distance(Q, P, C, B1, B2, leaf, sqrD, I, S, K);

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 7:
    {
      prepare_lhs_index(leaf,plhs+6);
      // Fall through
    }
    case 6:
    {
      prepare_lhs_double(B2,plhs+5);
      // Fall through
    }
    case 5:
    {
      prepare_lhs_double(B1,plhs+4);
      // Fall through
    }
    case 4:
    {
      prepare_lhs_double(K,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_double(S,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_double(I,plhs+1);
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



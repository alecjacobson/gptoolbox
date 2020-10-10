
#include <Eigen/Core>
#include <iostream>
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/PI.h>
#include <igl/iterative_closest_point.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>


void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[]
)
{
  using namespace igl::matlab;
  Eigen::MatrixXd VX,VY;
  Eigen::MatrixXi FX,FY;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  mexErrMsgTxt(nrhs>=4,"nrhs should be >= 4");
  parse_rhs_double(prhs+0,VX);
  parse_rhs_index(prhs+1,FX);
  parse_rhs_double(prhs+2,VY);
  parse_rhs_index(prhs+3,FY);
  int num_iters = 100;
  int num_samples = 1000;

  {
    int i = 4;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("MaxIter",name)==0)
      {
        igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
        igl::matlab::validate_arg_double(i,nrhs,prhs,name);
        num_iters = (double)*mxGetPr(prhs[++i]);
      }
      else if(strcmp("NumSamples",name)==0)
      {
        igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
        igl::matlab::validate_arg_double(i,nrhs,prhs,name);
        num_samples = (double)*mxGetPr(prhs[++i]);
      }
      else
      {
        mexErrMsgTxt(false,C_STR("Unrecognized Parameter: "<<name));
      }
      i++;
    }
  }

  Eigen::Matrix3d R;
  Eigen::RowVector3d t;
  igl::iterative_closest_point(VX,FX,VY,FY,num_samples,num_iters,R,t);

  switch(nlhs)
  {
    case 3:
    {
      Eigen::MatrixXd VXRT = (VX*R).rowwise()+t;
      prepare_lhs_double(VXRT,plhs+2);
    }
    case 2:
      prepare_lhs_double(t,plhs+1);
    case 1:
      prepare_lhs_double(R,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

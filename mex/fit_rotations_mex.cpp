#include <mex.h>
#include <igl/C_STR.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/fit_rotations.h>

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  // If no args then assume we're checking whether function exists, return true
  if(nrhs == 0)
  {
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    double * lhs = mxGetPr(plhs[0]);
    *lhs = true;
    return;
  }

  MatrixXd S,R;

  parse_rhs_double(prhs,S);
  const size_t dim = S.cols();
  mexErrMsgTxt(dim==2 || dim==3,"dim must be 2 or 3");
  mexErrMsgTxt(nrhs>=1,"must provide S arg");
  mexErrMsgTxt((S.rows()%dim) == 0,"S must be #R*dim by dim");
  const size_t nr = S.rows()/dim;
  bool allow_flips = false;
  bool single_precision = false;

  if(nrhs>2)
  {
    int i = 1;
    while(i<nrhs)
    {
      if(!mxIsChar(prhs[i]))
      {
        mexErrMsgTxt("Parameter names should be char strings");
      }
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("AllowFlips",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        allow_flips = *v;
      }
      if(strcmp("SinglePrecision",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        single_precision = *v;
      }else
      {
        mexErrMsgTxt(C_STR("Unsupported parameter: "<<name));
      }
      i++;
    }
  }

  mexErrMsgTxt(!allow_flips,"allow_flips not supported");

  if(dim == 2)
  {
    igl::fit_rotations(S,single_precision,R);
  }else
  {
    if(single_precision)
    {
#ifdef __SSE__
      MatrixXf Sf = S.cast<float>();
      MatrixXf Rf;
      fit_rotations_SSE(Sf,Rf);
      R = Rf.cast<double>();
#else
      igl::fit_rotations(S,single_precision,R);
#endif
    }else
    {
      igl::fit_rotations(S,false,R);
    }
  }

  // Prepare left-hand side
  switch(nlhs)
  {
    case 1:
    {
      plhs[0] = mxCreateDoubleMatrix(R.rows(),R.cols(), mxREAL);
      double * Rp = mxGetPr(plhs[0]);
      copy(&R.data()[0],&R.data()[0]+R.size(),Rp);
      break;
    }
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);

}


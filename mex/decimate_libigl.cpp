#include <igl/decimate.h>
#include <igl/qslim.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

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
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd V,W;
  MatrixXi F,G;
  VectorXi J,I;

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
  mexErrMsgTxt((ratio>0 && ratio<1) || (ratio>0 && ratio<F.rows()) ,
    "Ratio should be in (0,1) or [1,#F)");
  const size_t max_m = ratio<1 ? ratio*F.rows() : ratio;


  enum DecimateMethod
  {
    DECIMATE_METHOD_NAIVE = 0,
    DECIMATE_METHOD_QSLIM = 1,
    NUM_DECIMATE_METHODS = 2
  } method = DECIMATE_METHOD_NAIVE;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Method",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * type_name = mxArrayToString(prhs[++i]);
        if(strcmp("naive",type_name)==0)
        {
          method = DECIMATE_METHOD_NAIVE;
        }else if(strcmp("qslim",type_name)==0)
        {
          method = DECIMATE_METHOD_QSLIM;
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown method: "<<method));
        }
      }else
      {
        mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
      }
      i++;
    }
  }


  switch(method)
  {
    case DECIMATE_METHOD_NAIVE:
      decimate(V,F,max_m,W,G,J,I);
      break;
    case DECIMATE_METHOD_QSLIM:
      qslim(V,F,max_m,W,G,J,I);
      break;
    default:
      mexErrMsgTxt(false,"Unkown method.");
      break;
  }

  switch(nlhs)
  {
    case 4:
      prepare_lhs_index(I,plhs+3);
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

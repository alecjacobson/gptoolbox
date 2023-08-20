#include <igl/exact_geodesic.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd V;
  MatrixXi F;
  MatrixXi VS,FS;
  MatrixXi VT,FT;
    
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  
  mexErrMsgTxt(nrhs>=6,"nrhs should be >= 6");
  parse_rhs_double(prhs,V); 
  parse_rhs_index(prhs+1,F);
  parse_rhs_index(prhs+2,VS);
  parse_rhs_index(prhs+3,FS);
  parse_rhs_index(prhs+4,VT);
  parse_rhs_index(prhs+5,FT);
  mexErrMsgTxt(V.cols()==3 || V.cols()==2,"V must be #V by 2|3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  Eigen::VectorXd D;
  igl::exact_geodesic(V,F,VS,FS,VT,FT,D);

  switch(nlhs)
  {
      case 1:
          prepare_lhs_double(D,plhs+0);
      default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

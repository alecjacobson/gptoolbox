#include <igl/copyleft/cgal/trim_with_solid.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/C_STR.h>

#include <mex.h>
#include <Eigen/Dense>

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  mexErrMsgTxt(nrhs>=4,"Must have four inputs");
  Eigen::MatrixXd VA,VB;
  Eigen::MatrixXi FA,FB;
  const auto & parse_mesh = [](
      const mxArray *prhs[],
      MatrixXd & V,
      MatrixXi & F)
  {
    parse_rhs_double(prhs,V);
    parse_rhs_index(prhs+1,F);
    mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
    mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  };
  parse_mesh(prhs+0,VA,FA);
  parse_mesh(prhs+2,VB,FB);
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi D,J;
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,V,F,D,J);
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 4:
    {
      prepare_lhs_index(J,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_logical(D,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(F,plhs+1);
      // Fall through
    }

    case 1:
    {
      prepare_lhs_double(V,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

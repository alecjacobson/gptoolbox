// 
//   V2 = [0 0 0;1 0 0;0 1 0;0 1.004 0];
//   F2 = [1 2 3;3 2 4];
//   [CF] = collapse_small_triangles(V,F,1e-3);
// 

#include <igl/collapse_small_triangles.h>
#include <igl/read_triangle_mesh.h>
#include <igl/EPS.h>
#include <igl/pathinfo.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#define IGL_REDRUM_NOOP
#include <igl/REDRUM.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/mexErrMsgTxt.h>

#include "mex.h"
#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  if(nrhs != 2 && nrhs != 3)
  {
    mexErrMsgTxt("The number of input arguments must be 2 or 3.");
  }
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  MatrixXd V;
  MatrixXi F;
  VectorXi IM;

  string prefix;
  bool use_obj_format = false;
  double eps = FLOAT_EPS;
  // This parses first two arguements
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,F);
  igl::matlab::mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  igl::matlab::mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  if(nrhs==3)
  {
    if(mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1)
    {
      mexErrMsgTxt("3rd argument should be scalar");
    }
    eps = *mxGetPr(prhs[2]);
  }


  // let's first try to merge small triangles
  {
    MatrixXd tempV;
    MatrixXi tempF;
    collapse_small_triangles(V,F,eps,tempF);
    //if(F.rows() == tempF.rows())
    //{
    //  cout<<GREENGIN("No small triangles detected.")<<endl;
    //}else
    //{
    //  cout<<BLUEGIN("Collapsed all small triangles, reducing "<<F.rows()<<
    //    " input triangles to "<<tempF.rows()<<" triangles.")<<endl;
    //}
    F=tempF;
  }

  // Prepare left-hand side
  nlhs = 1;

  // Treat indices as reals
  plhs[0] = mxCreateDoubleMatrix(F.rows(),F.cols(), mxREAL);
  double * Fp = mxGetPr(plhs[0]);
  MatrixXd Fd = (F.cast<double>().array()+1).matrix();
  copy(&Fd.data()[0],&Fd.data()[0]+Fd.size(),Fp);

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);

}



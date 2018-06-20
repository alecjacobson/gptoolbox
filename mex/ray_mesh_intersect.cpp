#include <igl/embree/EmbreeIntersector.h>
#include <igl/Hit.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <Eigen/Core>
#include <limits>
#include <mex.h>

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
  using namespace igl::embree;
  // Embree uses singles
  MatrixXf V,source,dir;
  MatrixXi F;
  mexErrMsgTxt(nrhs==4,"nrhs should == 4");

  parse_rhs_double(prhs+0,source);
  parse_rhs_double(prhs+1,dir);
  parse_rhs_double(prhs+2,V);
  parse_rhs_index( prhs+3,F);
  mexErrMsgTxt(source.cols()==3,"source must be #source by 3");
  mexErrMsgTxt(dir.cols()==3,"dir must be #dir by 3");
  mexErrMsgTxt(source.rows() == dir.rows(), "#source must equal #dir");
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  const int n = source.rows();
  Eigen::Matrix<int, Eigen::Dynamic, 1> id(n);
  Eigen::Matrix<float, Eigen::Dynamic, 1> t(n);
  Eigen::Matrix<float, Eigen::Dynamic, 3> lambda(n,3);

  EmbreeIntersector ei;
  ei.init(V,F,true);
  for(int si = 0;si<n;si++)
  {
    Eigen::Vector3f s = source.row(si);
    Eigen::Vector3f d = dir.row(si);
    igl::Hit hit;
    const float tnear = 1e-4f;
    if(ei.intersectRay(s,d,hit,tnear))
    {
      id(si) = hit.id;
      t(si) = hit.t;
      lambda(si,0) = 1.0-hit.u-hit.v;
      lambda(si,1) = hit.u;
      lambda(si,2) = hit.v;
    }else
    {
      id(si) = -1;
      t(si) = std::numeric_limits<float>::infinity();
      lambda.row(si).setZero();
    }
  }


  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_double(lambda,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_double(t,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_index(id,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

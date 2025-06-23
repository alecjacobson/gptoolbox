#include <igl/embree/EmbreeIntersector.h>
#include <igl/Hit.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/parallel_for.h>

#include <Eigen/Core>
#include <limits>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
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

  EmbreeIntersector ei;
  ei.init(V,F,true);
  std::vector<std::vector<igl::Hit<float>>> hits(n);
  //for(int si = 0;si<n;si++)
  igl::parallel_for(n,[&](const int si)
  {
    Eigen::Vector3f s = source.row(si);
    Eigen::Vector3f d = dir.row(si);
    igl::Hit<float> hit;
    const float tnear = 1e-4f;
    int num_rays_shot;
    ei.intersectRay(s,d,hits[si],num_rays_shot,tnear);
  }
  );
  // total number of hits
  int total_hits = 0;
  for (int i = 0; i < n; ++i)
  {
    total_hits += hits[i].size();
  }
  Eigen::VectorXi I(total_hits);
  Eigen::VectorXi J(total_hits);
  Eigen::VectorXd T(total_hits);
  Eigen::MatrixXd lambda(total_hits,3);


  int si = 0;
  for(int i = 0; i < n; ++i)
  {
    for(const auto& hit : hits[i])
    {
      I(si) = i;
      J(si) = hit.id;
      T(si) = hit.t;
      lambda(si,0) = 1.0-hit.u-hit.v;
      lambda(si,1) = hit.u;
      lambda(si,2) = hit.v;
      si++;
    }
  }

  


  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 4:
    {
      prepare_lhs_double(lambda,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_double(T,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_double(J,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_index(I,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

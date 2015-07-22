// 
// mex -v -o ambient_occlusion -DMEX -largeArrayDims ...
//   -I/opt/local/include/eigen3 -I/usr/local/igl/libigl/include ...
//   -L/usr/local/igl/libigl/lib -ligl -liglmatlab -liglembree ...
//   CXXFLAGS="\$CXXFLAGS -m64 -msse4.2 -fopenmp" ...
//   -I/usr/local/igl/libigl/external/embree/include/ ...
//   -I/usr/local/igl/libigl/external/embree/ ...
//   -L/usr/local/igl/libigl/external/embree/build -lembree -lsys ... 
//   ambient_occlusion.cpp
//

#include <igl/matlab/MexStream.h>
#include <igl/embree/ambient_occlusion.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/Core>
#include <mex.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>

// Parse right hand side arguments for a matlab mex function.
//
// Inputs:
//   nrhs  number of right hand side arguments
//   prhs  pointer to right hand side arguments
// Outputs:
//   V  n by dim list of mesh vertex positions
//   F  m by dim list of mesh face indices
//   P  #P by 3 list of origin points
//   N  #P by 3 list of origin normals
//   num_samples  number of samples
// Throws matlab errors if dimensions are not sane.
void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N,
  int & num_samples)
{
  using namespace std;
  if(nrhs < 5)
  {
    mexErrMsgTxt("nrhs < 5");
  }

  const int dim = mxGetN(prhs[0]);
  if(dim != 3)
  {
    mexErrMsgTxt("Mesh vertex list must be #V by 3 list of vertex positions");
  }
  if(dim != (int)mxGetN(prhs[1]))
  {
   mexErrMsgTxt("Mesh facet size must be 3");
  }
  if(mxGetN(prhs[2]) != dim)
  {
    mexErrMsgTxt("Point list must be #P by 3 list of origin locations");
  }
  if(mxGetN(prhs[3]) != dim)
  {
    mexErrMsgTxt("Normal list must be #P by 3 list of origin normals");
  }
  if(mxGetN(prhs[4]) != 1 || mxGetM(prhs[4]) != 1)
  {
    mexErrMsgTxt("Number of samples must be scalar.");
  }


  V.resize(mxGetM(prhs[0]),mxGetN(prhs[0]));
  copy(mxGetPr(prhs[0]),mxGetPr(prhs[0])+V.size(),V.data());
  F.resize(mxGetM(prhs[1]),mxGetN(prhs[1]));
  copy(mxGetPr(prhs[1]),mxGetPr(prhs[1])+F.size(),F.data());
  F.array() -= 1;
  P.resize(mxGetM(prhs[2]),mxGetN(prhs[2]));
  copy(mxGetPr(prhs[2]),mxGetPr(prhs[2])+P.size(),P.data());
  N.resize(mxGetM(prhs[3]),mxGetN(prhs[3]));
  copy(mxGetPr(prhs[3]),mxGetPr(prhs[3])+N.size(),N.data());
  if(*mxGetPr(prhs[4]) != (int)*mxGetPr(prhs[4]))
  {
    mexErrMsgTxt("Number of samples should be non negative integer.");
  }
  num_samples = (int) *mxGetPr(prhs[4]);
}

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
  MatrixXd V,P,N;
  VectorXd S;
  MatrixXi F;
  int num_samples;
  parse_rhs(nrhs,prhs,V,F,P,N,num_samples);
  // Prepare left-hand side
  nlhs = 1;
  igl::embree::ambient_occlusion(V,F,P,N,num_samples,S);
  plhs[0] = mxCreateDoubleMatrix(S.rows(),S.cols(), mxREAL);
  copy(S.data(),S.data()+S.size(),mxGetPr(plhs[0]));
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

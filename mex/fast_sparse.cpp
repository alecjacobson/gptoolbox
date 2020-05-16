#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>

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

  mexErrMsgTxt(nrhs==5,"nrhs should be == 5");
  // unnecessary copy
  int m,n;
  m = (double)*mxGetPr(prhs[3]);
  n = (double)*mxGetPr(prhs[4]);
  mexErrMsgTxt(mxGetM(prhs[0]) == mxGetM(prhs[1]) && mxGetM(prhs[0]) == mxGetM(prhs[2]),"I J V must be same #rows");
  mexErrMsgTxt(mxGetN(prhs[0]) == mxGetN(prhs[1]) && mxGetN(prhs[0]) == mxGetN(prhs[2]),"I J V must be same #cols");
  const int nr = mxGetM(prhs[0]);
  const int nc = mxGetN(prhs[0]);
  std::vector<Eigen::Triplet<double> > IJV;
  IJV.reserve(nr*nc);
  for(int i = 0;i<nr;i++)
  {
    for(int j = 0;j<nc;j++)
    {
      IJV.emplace_back(
        mxGetPr(prhs[0])[i+nr*j]-1,
        mxGetPr(prhs[1])[i+nr*j]-1,
        mxGetPr(prhs[2])[i+nr*j]);
    }
  }
  Eigen::SparseMatrix<double> S(m,n);
  // S.reserve doesn't seem to help much
  S.setFromTriplets(IJV.begin(),IJV.end());
  // This final copy is not that expensive (e.g., 0.5s gather, 1.1s scatter,
  // 0.05s prepare lhs)
  switch(nlhs)
  {
    case 1:
      prepare_lhs_double(S,plhs+0);
    default:break;
  }
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

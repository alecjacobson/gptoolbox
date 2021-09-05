#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <Eigen/Core>

template <int _N>
void psd_project_rows(Eigen::MatrixXd & H)
{
  const Eigen::Index m = H.rows();
  const int N = _N==Eigen::Dynamic ? round(sqrt(H.cols())) : _N;
  static_assert( _N == Eigen::Dynamic || _N == N,"");
  for(Eigen::Index r = 0;r<m;r++)
  {
    typedef Eigen::Matrix<double,_N,_N> Matrix;
    Matrix Hr(N,N);
    for(int i=0;i<N;i++)
    {
      for(int j=i;j<N;j++)
      {
        Hr(i,j) = H(r,i+N*j);
        Hr(j,i) = H(r,i+N*j);
      }
    }
    Eigen::SelfAdjointEigenSolver<Matrix> es(Hr);
    Hr = 
      (es.eigenvectors()*
      es.eigenvalues().cwiseMax(0).asDiagonal() *
      es.eigenvectors().transpose()).eval();
    for(int i=0;i<N;i++)
    {
      for(int j=i;j<N;j++)
      {
        H(r,i+N*j) = Hr(i,j);
        H(r,j+N*i) = Hr(i,j);
      }
    }

  }
}

// https://stackoverflow.com/a/1549960/148668
bool is_perfect_square(int n) {
    if (n < 0)
        return false;
    int root(round(sqrt(n)));
    return n == root * root;
}


void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[]
  )
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=1,"nrhs should be == 1");
  Eigen::MatrixXd H;
  parse_rhs_double(prhs+0,H);

  mexErrMsgTxt(is_perfect_square(H.cols()),"size(H,2) must be perfect square");

  const int N(round(sqrt(H.cols())));
  switch(N)
  {
    case  2: psd_project_rows< 2>(H); break;
    case  3: psd_project_rows< 3>(H); break;
    case  4: psd_project_rows< 4>(H); break;
    case  5: psd_project_rows< 5>(H); break;
    case  6: psd_project_rows< 6>(H); break;
    case  7: psd_project_rows< 7>(H); break;
    case  8: psd_project_rows< 8>(H); break;
    case  9: psd_project_rows< 9>(H); break;
    case 10: psd_project_rows<10>(H); break;
    case 11: psd_project_rows<11>(H); break;
    case 12: psd_project_rows<12>(H); break;
    default:
      psd_project_rows<Eigen::Dynamic>(H);
      break;
  }


  switch(nlhs)
  {
    case 1:
      prepare_lhs_double(H,plhs+0);
    default:break;
  }
  std::cout.rdbuf(outbuf);
}


#include <Eigen/Core>
#include <iostream>
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <cyPolynomial.h>

template<int N>
bool call_PolynomialRoots(
  int deg,
  Eigen::RowVectorXd &R,
  const Eigen::RowVectorXd &coef,
  const double xmin,
  const double xmax)
{
  using cy::PolynomialRoots;
  if constexpr (N == 0)
  {
    return false;
  }
  else
  {
    if (deg == N)
    {
      double r[N];
      int nr = PolynomialRoots<N>(r, coef.data(),xmin,xmax);
      R.setConstant(N,std::numeric_limits<double>::quiet_NaN());
      for (int i = 0; i < nr; ++i)
        R[i] = r[i];
      return true;
    }
    return call_PolynomialRoots<N - 1>(deg, R, coef,xmin,xmax);
  }
}

Eigen::RowVectorXd roots(const Eigen::RowVectorXd &coef, const double xmin, const double xmax)
{
  constexpr int MAX_DEG = 16;
  const int deg = coef.size() - 1;

  if (deg < 1 || deg > MAX_DEG)
    throw std::runtime_error("Polynomial degree out of range");

  Eigen::RowVectorXd R;
  if (!call_PolynomialRoots<MAX_DEG>(deg, R, coef, xmin, xmax))
    throw std::runtime_error("Root computation failed");

  return R;
}


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

  mexErrMsgTxt(nrhs==1 || nrhs==3,"nrhs should be == 1 | 3");
  Eigen::MatrixXd P;
  parse_rhs_double(prhs+0,P);
  Eigen::VectorXd Xmin(1);Xmin(0) = std::numeric_limits<double>::lowest();
  Eigen::VectorXd Xmax(1);Xmax(0) = std::numeric_limits<double>::max();
  const int n = P.rows();
  if(nrhs>1)
  {
    parse_rhs_double(prhs+1,Xmin);
    parse_rhs_double(prhs+2,Xmax);
    mexErrMsgTxt(Xmin.size()==1 || Xmin.size()==n,"Xmin size invalid");
    mexErrMsgTxt(Xmax.size()==1 || Xmax.size()==n,"Xmax size invalid");
  }
  const int deg = P.cols()-1;
  Eigen::MatrixXd X(n,deg-1);
  for(int i = 0;i<n;i++)
  {
    const Eigen::RowVectorXd coef = P.row(i).reverse();
    const double xmin = Xmin(Xmin.size()==1?0:i);
    const double xmax = Xmax(Xmax.size()==1?0:i);
    X.row(i) = roots(coef,xmin,xmax);
  }

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 1:
    {
      prepare_lhs_double(X,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

// https://igl.ethz.ch/projects/LB3D/dualLaplace.cpp
#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <Eigen/Core>
#include <Eigen/Sparse>


void circumcenter(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, Eigen::Vector3d& cc)
{
    const double l[3]{
        (b - c).squaredNorm(),
        (a - c).squaredNorm(),
        (a - b).squaredNorm()
    };
    
    const double ba[3]{l[0] * (l[1] + l[2] - l[0]), l[1] * (l[2] + l[0] - l[1]), l[2] * (l[0] + l[1] - l[2])};
    const double sum = ba[0] + ba[1] + ba[2];
    
    cc = (ba[0] / sum) * a + (ba[1] / sum) * b + (ba[2] / sum) * c;
}

void circumcenter(const Eigen::Matrix<double, 4, 3>& t, Eigen::Vector3d& c)
{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    
    const double n0 = t.row(0).squaredNorm();
    
    for(int k = 0; k < 3; ++k)
    {
        A.row(k) = t.row(k + 1) - t.row(0);
        b(k) = t.row(k + 1).squaredNorm() - n0;
    }
    
    c = 0.5 * A.fullPivHouseholderQr().solve(b);
}

double volume(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d)
{
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;
    
    return A.determinant() / 6.;
}

void dualLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>& M)
{
    const size_t nt = T.rows();
    const size_t nv = V.rows();
    
    const int turn[4][4]
    {
        {-1, 2, 3, 1},
        {3, -1, 0, 2},
        {1, 3, -1, 0},
        {2, 0, 1, -1}
    };

    auto getTet = [&](const int i, Eigen::Matrix<double, 4, 3>& t)
    {
        for(int k = 0; k < 4; ++k)
        {
            t.row(k) = V.row(T(i, k));
        }
    };
    
    std::vector<Eigen::Triplet<double>> tripL, tripM;
    
    Eigen::Vector3d cc;
    Eigen::Matrix<double, 4, 3> t;
    
    for(int k = 0; k < nt; ++k)
    {
        getTet(k, t);
        circumcenter(t, cc);
        
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(i != j)
                {
                    Eigen::Vector3d cf;
                    circumcenter(t.row(i), t.row(j), t.row(turn[i][j]), cf);
                    
                    const Eigen::Vector3d ce = 0.5 * (t.row(i) + t.row(j));
                    
                    const double vol = volume(t.row(i), ce, cf, cc);
                    const double wij = 6. * vol / (t.row(i) - t.row(j)).squaredNorm();
                    
                    tripL.emplace_back(T(k, i), T(k, j), wij);
                    tripL.emplace_back(T(k, j), T(k, i), wij);
                    
                    tripL.emplace_back(T(k, i), T(k, i), -wij);
                    tripL.emplace_back(T(k, j), T(k, j), -wij);
                    
                    tripM.emplace_back(T(k, i), T(k, i), vol);
                    tripM.emplace_back(T(k, j), T(k, j), vol);
                }
            }
        }
    }
    
    L.resize(nv, nv);
    M.resize(nv, nv);
    
    L.setFromTriplets(tripL.begin(), tripL.end());
    M.setFromTriplets(tripM.begin(), tripM.end());
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

  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  mexErrMsgTxt(nrhs>=2,"nrhs should be == 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,T);
  Eigen::SparseMatrix<double> L,M;
  mexErrMsgTxt(T.cols() == 4, "T should be indices into a tet-mesh");
  dualLaplace(V,T,L,M);

  switch(nlhs)
  {
    case 2:
      prepare_lhs_double(M,plhs+1);
    case 1:
      prepare_lhs_double(L,plhs+0);
    default:break;
  }
  std::cout.rdbuf(outbuf);
}

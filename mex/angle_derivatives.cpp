
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Geometry>
// From libshell by Etienne Vouga
double angle(const Eigen::Vector3d &v, const Eigen::Vector3d &w, const Eigen::Vector3d &axis,
    Eigen::Matrix<double, 1, 6> *derivative, // v, w
    Eigen::Matrix<double, 6, 6> *hessian
)
{
  const auto crossMatrix = [](Eigen::Vector3d v)->Eigen::Matrix3d 
  {
      Eigen::Matrix3d ret;
      ret << 0, -v(2), v(1),
          v(2), 0, -v(0),
          -v(1), v(0), 0;
      return ret;
  };

  // This is a bad idea, if v and w are colinear, this derivatives below still
  // make sense but this will make axis 0 and nuke everything.
  //// This is unnecessary if the caller would promise that axis is unit and
  //// orthogonal to v and w
  const Eigen::Vector3d vcrossw = v.cross(w);
  //const Eigen::Vector3d axis = ((vcrossw.dot(in_axis)>0?1:-1) * vcrossw).eval().normalized();
  //std::cout<<"in: "<<in_axis.transpose()<<std::endl;
  //std::cout<<"vcrossw: "<<vcrossw.transpose()<<std::endl;
  //std::cout<<"ef: "<<axis.transpose()<<std::endl;

    double theta = 2.0 * atan2((vcrossw.dot(axis)), v.dot(w) + v.norm() * w.norm());

    if (derivative)
    {
        derivative->segment(0, 3) = -axis.cross(v) / v.squaredNorm();
        derivative->segment(3, 3) = axis.cross(w) / w.squaredNorm();
    }
    if (hessian)
    {
        hessian->setZero();
        hessian->block(0, 0, 3, 3) += 2.0 * (axis.cross(v))*v.transpose() / v.squaredNorm() / v.squaredNorm();
        hessian->block(3, 3, 3, 3) += -2.0 * (axis.cross(w))*w.transpose() / w.squaredNorm() / w.squaredNorm();
        hessian->block(0, 0, 3, 3) += -crossMatrix(axis) / v.squaredNorm();
        hessian->block(3, 3, 3, 3) += crossMatrix(axis) / w.squaredNorm();

        double sigma = 1.0;
        if (v.cross(w).dot(axis) < 0)
            sigma = -1.0;

        double vwnorm = v.cross(w).norm();
        if (vwnorm > 1e-8)
        {
            Eigen::Matrix3d da = sigma * (1.0 / vwnorm * Eigen::Matrix3d::Identity() - 1.0 / vwnorm / vwnorm / vwnorm * (v.cross(w)) * (v.cross(w)).transpose());
            hessian->block(0, 0, 3, 3) += crossMatrix(v) / v.squaredNorm() * da * -crossMatrix(w);
            hessian->block(3, 0, 3, 3) += crossMatrix(v) / v.squaredNorm() * da * crossMatrix(v);
            hessian->block(0, 3, 3, 3) += -crossMatrix(w) / w.squaredNorm() * da * -crossMatrix(w);
            hessian->block(3, 3, 3, 3) += -crossMatrix(w) / w.squaredNorm() * da * crossMatrix(v);
        }
    }

    return theta;
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

  Eigen::MatrixXd V,W,A;
  mexErrMsgTxt(nrhs>=2,"nrhs should be == 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_double(prhs+1,W);
  mexErrMsgTxt(V.cols()==W.cols(),"dims should be the same");
  mexErrMsgTxt(V.rows()==W.rows(),"dims should be the same");
  const int n = V.rows();
  const int dim = V.cols();
  if(dim == 3)
  {
    mexErrMsgTxt(nrhs>=3,"nrhs should be == 3");
    parse_rhs_double(prhs+2,A);
  }else if(dim == 2)
  {
    V.conservativeResize(V.rows(),3);
    W.conservativeResize(W.rows(),3);
    V.col(2).setZero();
    W.col(2).setZero();
    A.resize(n,3);
    A.col(0).setConstant(0);
    A.col(1).setConstant(0);
    A.col(2).setConstant(1);
  }
  Eigen::MatrixXd theta(n,1);
  Eigen::MatrixXd dthetadV(n,dim);
  Eigen::MatrixXd dthetadW(n,dim);
  Eigen::MatrixXd d2thetadV2(n,dim==2?3:6);
  Eigen::MatrixXd d2thetadW2(n,dim==2?3:6);
  Eigen::MatrixXd d2thetadVW(n,dim==2?3:6);
  const Eigen::MatrixXi voigt2 = 
    (Eigen::MatrixXi(3,2)<<0,0,1,1,0,1).finished();
  const Eigen::MatrixXi voigt = (dim==2? voigt2 :
    (Eigen::MatrixXi(6,2)<<0,0,1,1,2,2,1,2,0,2,0,1).finished());
  for(int i = 0;i<n;i++)
  {
    Eigen::Matrix<double, 1, 6> derivative;
    Eigen::Matrix<double, 6, 6> hessian;
    theta(i) = angle(
        V.row(i).transpose(),W.row(i).transpose(),A.row(i).transpose(),
        &derivative,&hessian);
    dthetadV.row(i) = derivative.segment(0,dim);
    dthetadW.row(i) = derivative.segment(3,dim);
    for(int b = 0;b<voigt2.rows();b++)
    {
      const int bi = voigt2(b,0);
      const int bj = voigt2(b,1);
      const auto & block = hessian.block(bi*3,bj*3,dim,dim);
      for(int v = 0;v<voigt.rows();v++)
      {
        const int vi = voigt(v,0);
        const int vj = voigt(v,1);
        const double h = block(vi,vj);
        switch(b)
        {
          case 0: d2thetadV2(i,v) = h; break;
          case 1: d2thetadW2(i,v) = h; break;
          case 2: d2thetadVW(i,v) = h; break;
        }
      }
    }
  }

  switch(nlhs)
  {
    case 6:
      prepare_lhs_double(d2thetadVW,plhs+5);
    case 5:
      prepare_lhs_double(d2thetadW2,plhs+4);
    case 4:
      prepare_lhs_double(d2thetadV2,plhs+3);
    case 3:
      prepare_lhs_double(dthetadW,plhs+2);
    case 2:
      prepare_lhs_double(dthetadV,plhs+1);
    case 1:
      prepare_lhs_double(theta,plhs+0);
    default:break;
  }


  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

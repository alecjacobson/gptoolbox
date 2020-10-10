#include <igl/AABB.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>

#include <cstring>
#include <iostream>
void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & Ele,
  Eigen::MatrixXd & bb_mins,
  Eigen::MatrixXd & bb_maxs,
  Eigen::VectorXi & elements)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3 || dim == 2,
    "Mesh vertex list must be #V by 2 or 3 list of vertex positions");

  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,Ele);
  if(nrhs > 2)
  {
    mexErrMsgTxt(nrhs == 2, "The number of input arguments must be 2.");
  }
}

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

  MatrixXd V,bb_mins,bb_maxs;
  MatrixXi Ele;
  VectorXi elements;
  VectorXi I;
  parse_rhs(nrhs,prhs,V,Ele,bb_mins,bb_maxs,elements);
  bool was_serialized = bb_mins.size()>0;

  switch(V.cols())
  {
    default:
      mexErrMsgTxt(false,"Un-supported dimension.");
      break;
    case 3:
    {
      AABB<MatrixXd,3> aabb;
      aabb.init(V,Ele,bb_mins,bb_maxs,elements);
      aabb.serialize(bb_mins,bb_maxs,elements);
      break;
    }
    case 2:
    {
      AABB<MatrixXd,2> aabb;
      aabb.init(V,Ele,bb_mins,bb_maxs,elements);
      aabb.serialize(bb_mins,bb_maxs,elements);
      break;
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
      prepare_lhs_index(elements,plhs+2);
      //fallthrough
    }
    case 2:
    {
      prepare_lhs_double(bb_maxs,plhs+1);
      //fallthrough
    }
    case 1:
    {
      prepare_lhs_double(bb_mins,plhs+0);
      //fallthrough
    }
    case 0: break;
  }

  std::cout.rdbuf(outbuf);
}

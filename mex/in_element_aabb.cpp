#ifdef MEX

#include <igl/AABB.h>
#include <igl/in_element.h>
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
  Eigen::MatrixXd & Q,
  Eigen::MatrixXd & bb_mins,
  Eigen::MatrixXd & bb_maxs,
  Eigen::VectorXi & elements)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3 || dim == 2,
    "Mesh vertex list must be #V by 2 or 3 list of vertex positions");

  mexErrMsgTxt(dim+1 == mxGetN(prhs[1]),
    "Mesh \"face\" simplex size must equal dimension+1");

  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,Ele);
  parse_rhs_double(prhs+2,Q);
  mexErrMsgTxt(Q.cols() == dim,"Dimension of Q should match V");
  if(nrhs > 3)
  {
    mexErrMsgTxt(nrhs >= 6, "The number of input arguments must be 3 or >=6.");
    parse_rhs_double(prhs+3,bb_mins);
    if(bb_mins.size()>0)
    {
      mexErrMsgTxt(bb_mins.cols() == dim,"Dimension of bb_mins should match V");
      mexErrMsgTxt(bb_mins.rows() >= Ele.rows(),"|bb_mins| should be > |Ele|");
    }
    parse_rhs_double(prhs+4,bb_maxs);
    mexErrMsgTxt(bb_maxs.cols() == bb_mins.cols(),
      "|bb_maxs| should match |bb_mins|");
    mexErrMsgTxt(bb_mins.rows() == bb_maxs.rows(),
      "|bb_mins| should match |bb_maxs|");
    parse_rhs_index(prhs+5,elements);
    mexErrMsgTxt(elements.cols() == 1,"Elements should be column vector");
    mexErrMsgTxt(bb_mins.rows() == elements.rows(),
      "|bb_mins| should match |elements|");
  }else
  {
    // Defaults
    bb_mins.resize(0,dim);
    bb_maxs.resize(0,dim);
    elements.resize(0,1);
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

  MatrixXd V,Q,bb_mins,bb_maxs;
  MatrixXi Ele;
  VectorXi elements;
  VectorXi I;
  parse_rhs(nrhs,prhs,V,Ele,Q,bb_mins,bb_maxs,elements);
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
      in_element(V,Ele,Q,aabb,I);
      if(nlhs>1 && !was_serialized)
      {
        aabb.serialize(bb_mins,bb_maxs,elements);
      }
      break;
    }
    case 2:
    {
      AABB<MatrixXd,2> aabb;
      aabb.init(V,Ele,bb_mins,bb_maxs,elements);
      in_element(V,Ele,Q,aabb,I);
      if(nlhs>1 && !was_serialized)
      {
        aabb.serialize(bb_mins,bb_maxs,elements);
      }
      break;
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
      prepare_lhs_index(elements,plhs+3);
      //fallthrough
    }
    case 3:
    {
      prepare_lhs_double(bb_maxs,plhs+2);
      //fallthrough
    }
    case 2:
    {
      prepare_lhs_double(bb_mins,plhs+1);
      //fallthrough
    }
    case 1:
    {
      prepare_lhs_index(I,plhs+0);
      //fallthrough
    }
    case 0: break;
  }

  std::cout.rdbuf(outbuf);
}

#endif

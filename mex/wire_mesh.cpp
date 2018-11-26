#include <igl/avg_edge_length.h>
#include <igl/edges.h>
#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>

#include <cstring>

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & WV,
  Eigen::MatrixXi & WE,
  double & th,
  int & poly_size,
  bool & solid)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 2, "The number of input arguments must be >=2.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3,
    "Mesh vertex list must be #V by 3 list of vertex positions");

  parse_rhs_double(prhs,WV);
  parse_rhs_index(prhs+1,WE);
  if(WE.cols()>2)
  {
    igl::edges(MatrixXi(WE),WE);
  }

  // defaults
  // Thickness
  th = 0.1*igl::avg_edge_length(WV,WE);
  // Size of extrusion polygon
  poly_size = 4;
  solid = true;

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Thickness",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        th = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("PolySize",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        poly_size = (int)*mxGetPr(prhs[++i]);
      }else if(strcmp("Solid",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        solid = (bool)*mxGetPr(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
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
  using namespace igl::copyleft::cgal;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd WV,V;
  MatrixXi WE,F;
  VectorXi J;
  double th;
  int poly_size;
  bool solid;
  parse_rhs(
    nrhs,prhs,
    WV,WE,
    th,poly_size,solid);
  wire_mesh(WV,WE,th,poly_size,solid,V,F,J);
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_index(J,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(F,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(V,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

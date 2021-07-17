#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab_format.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/triangle/triangulate.h>
#include <igl/copyleft/cgal/triangulate.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXd H;
  Eigen::VectorXi VM,EM;
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,E);

  std::string flags = "";
  bool verbose = false;
  enum MethodType
  {
    TRIANGLE,
    CGAL_EPECK,
    CGAL_EPICK,
  } method = TRIANGLE;

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Flags",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        flags = mxArrayToString(prhs[++i]);
      }else if(strcmp("EdgeMarkers",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_index(prhs+(++i),EM);
      }else if(strcmp("Holes",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_double(prhs+(++i),H);
      }else if(strcmp("VertexMarkers",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_index(prhs+(++i),VM);
      }else if(strcmp("Method",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * method_name = mxArrayToString(prhs[++i]);
        if(strcmp("triangle",method_name)==0)
        {
          method = TRIANGLE;
        }else if(strcmp("epeck",method_name)==0)
        {
          method = CGAL_EPECK;
#ifndef WITH_CGAL
          mexErrMsgTxt(false,"Recompile with cgal to enabel.");
#endif
        }else if(strcmp("epick",method_name)==0)
        {
          method = CGAL_EPICK;
#ifndef WITH_CGAL
          mexErrMsgTxt(false,"Recompile with cgal to enabel.");
#endif
        }else if(strcmp("cgal",method_name)==0)
        {
          method = CGAL_EPECK;
#ifndef WITH_CGAL
          mexErrMsgTxt(false,"Recompile with cgal to enabel.");
#endif
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown Method: "<<method_name));
        }
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }


  Eigen::MatrixXd TV;
  Eigen::MatrixXi TF;
  Eigen::MatrixXi TE;
  Eigen::VectorXi TVM,TEM;
  // triangle will annoyingly mark all boundary edges (e.g., "-c" flag) with 1.
  // Shift all the input labels so that min is at 2. Then we'll shift back and
  // mark boundaries with EM.maxCoeff()+1
  int shift,max_EM;
  if(EM.size())
  {
    shift = EM.minCoeff()+2;
    max_EM = EM.maxCoeff();
    EM.array() += shift;
  }
  switch(method)
  {
    default:
    case TRIANGLE:
    try
    {
      igl::triangle::triangulate(V,E,H,VM,EM,flags,TV,TF,TVM,TE,TEM);
    }catch(std::runtime_error e)
    {
      ::mexErrMsgTxt((std::string("triangle failed: ")+e.what()).c_str());
    }
    break;
    case CGAL_EPECK:
      igl::copyleft::cgal::triangulate<CGAL::Epeck>(V,E,H,flags.find('c')!=std::string::npos,TV,TF);
    case CGAL_EPICK:
      igl::copyleft::cgal::triangulate<CGAL::Epick>(V,E,H,flags.find('c')!=std::string::npos,TV,TF);
    break;
  }
  for(int e = 0;e<TEM.size();e++)
  {
    if(TEM(e) == 1){ TEM(e) = max_EM+1;}
    else{ TEM(e) -= shift;}
  }

  switch(nlhs)
  {
    case 5:
      prepare_lhs_index(TEM,plhs+4);
    case 4:
      prepare_lhs_index(TE,plhs+3);
    case 3:
      prepare_lhs_index(TVM,plhs+2);
    case 2:
      prepare_lhs_index(TF,plhs+1);
    case 1:
      prepare_lhs_double(TV,plhs+0);
    default:break;
  }
}



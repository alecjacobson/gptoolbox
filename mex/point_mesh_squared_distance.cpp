#include <mex.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/C_STR.h>
#include <igl/copyleft/cgal/point_mesh_squared_distance.h>

#include <iostream>
#include <string>
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
  using namespace igl::matlab;

  MatrixXd P,V,C;
  VectorXi I;
  VectorXd sqrD;
  MatrixXi F;
  enum PointMeshSquaredDistanceMethod
  {
    POINT_MESH_SQUARED_DISTANCE_METHOD_LIBIGL = 0,
    POINT_MESH_SQUARED_DISTANCE_METHOD_CGAL = 1,
    NUM_POINT_MESH_SQUARED_DISTANCE_METHODS = 2
  } method = POINT_MESH_SQUARED_DISTANCE_METHOD_LIBIGL;
  if(nrhs < 3)
  {
    mexErrMsgTxt("nrhs < 3");
  }
  parse_rhs_double(prhs,P);
  parse_rhs_double(prhs+1,V);
  parse_rhs_index(prhs+2,F);
  mexErrMsgTxt(P.cols()==3 || P.cols()==2,"P must be #P by (3|2)");
  mexErrMsgTxt(V.cols()==3 || V.cols()==2,"V must be #V by (3|2)");
  mexErrMsgTxt(V.cols()==P.cols(),"dim(V) must be dim(P)");
  mexErrMsgTxt(F.cols()==3 || F.cols()==2 || F.cols()==1,"F must be #F by (3|2|1)");

  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Method",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * type_name = mxArrayToString(prhs[++i]);
        if(strcmp("libigl",type_name)==0)
        {
          method = POINT_MESH_SQUARED_DISTANCE_METHOD_LIBIGL;
        }else if(strcmp("cgal",type_name)==0)
        {
          method = POINT_MESH_SQUARED_DISTANCE_METHOD_CGAL;
          mexErrMsgTxt(F.cols() == 3, "'cgal' method only works for triangles");
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown method: "<<method));
        }
      }else
      {
        mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
      }
      i++;
    }
  }

  switch(method)
  {
    case POINT_MESH_SQUARED_DISTANCE_METHOD_LIBIGL:
      igl::point_mesh_squared_distance(P,V,F,sqrD,I,C);
      break;
    case POINT_MESH_SQUARED_DISTANCE_METHOD_CGAL:
      igl::copyleft::cgal::point_mesh_squared_distance<CGAL::Epeck>(P,V,F,sqrD,I,C);
      break;
    default:
      mexErrMsgTxt(false,"Unkown method.");
      break;
  }

  // Prepare left-hand side
  switch(nlhs)
  {
    case 3:
    {
      // Treat indices as reals
      prepare_lhs_double(C,plhs+2);
      // Fallthrough
    }
    case 2:
    {
      prepare_lhs_index(I,plhs+1);
      // Fallthrough
    }
    case 1:
    {
      prepare_lhs_double(sqrD,plhs+0);
      break;
    }
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

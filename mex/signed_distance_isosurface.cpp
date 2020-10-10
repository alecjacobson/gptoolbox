#include <igl/copyleft/cgal/signed_distance_isosurface.h>
#include <igl/copyleft/offset_surface.h>
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

enum ContouringMethod
{
  CONTOURING_METHOD_MARCHING_CUBES = 0,
  CONTOURING_METHOD_CGAL = 1,
  NUM_CONTOURING_METHOD = 2,
};

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & IV,
  Eigen::MatrixXi & IF,
  double & level,
  double & angle_bound,
  double & radius_bound,
  double & distance_bound,
  igl::SignedDistanceType & type,
  ContouringMethod & contouring_method,
  int & grid_size)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 2, "The number of input arguments must be >=2.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3,
    "Mesh vertex list must be #V by 3 list of vertex positions");

  parse_rhs_double(prhs,IV);
  parse_rhs_index(prhs+1,IF);

  // defaults
  level = 0.0;
  angle_bound = 28.0;
  double bbd = 1.0;
  // Radius and distance in terms of fraction of bbd
  if(IV.size() > 0)
  {
    bbd = (IV.colwise().maxCoeff()-IV.colwise().minCoeff()).norm();
  }
  radius_bound = 0.02*bbd;
  distance_bound = 0.02*bbd;
  type = SIGNED_DISTANCE_TYPE_DEFAULT;
  contouring_method = CONTOURING_METHOD_MARCHING_CUBES;
  grid_size = 40;

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Level",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        level = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("GridSize",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        grid_size = (int)*mxGetPr(prhs[++i]);
      }else if(strcmp("AngleBound",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        angle_bound = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("RadiusBound",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        radius_bound = ((double)*mxGetPr(prhs[++i])) * bbd;
      }else if(strcmp("DistanceBound",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        distance_bound = ((double)*mxGetPr(prhs[++i])) * bbd;
      }else if(strcmp("ContouringMethod",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * contouring_name = mxArrayToString(prhs[++i]);
        if(strcmp("cgal",contouring_name)==0)
        {
          contouring_method = CONTOURING_METHOD_CGAL;
        }else if(strcmp("marching_cubes",contouring_name)==0)
        {
          contouring_method = CONTOURING_METHOD_MARCHING_CUBES;
        }
      }else if(strcmp("SignedDistanceType",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * type_name = mxArrayToString(prhs[++i]);
        if(strcmp("pseudonormal",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
        }else if(strcmp("winding_number",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
        }else if(strcmp("default",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_DEFAULT;
        }else if(strcmp("unsigned",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_UNSIGNED;
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown SignedDistanceType: "<<type_name));
        }
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  if(type != igl::SIGNED_DISTANCE_TYPE_UNSIGNED)
  {
    mexErrMsgTxt(dim == mxGetN(prhs[1]),
      "Mesh \"face\" simplex size must equal dimension");
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

  MatrixXd IV,V;
  MatrixXi IF,F;
  double level,angle_bound,radius_bound,distance_bound;
  SignedDistanceType type;
  ContouringMethod contouring_method;
  int grid_size;
  parse_rhs(
    nrhs,prhs,
    IV,IF,
    level,angle_bound,radius_bound,distance_bound,type,
    contouring_method, grid_size);
  Eigen::MatrixXd GV;
  Eigen::RowVector3i side;
  Eigen::VectorXd S;
  switch(contouring_method)
  {
    default:
    case CONTOURING_METHOD_MARCHING_CUBES:
      {
        igl::copyleft::offset_surface(IV,IF,level,grid_size,type,V,F,GV,side,S);
      }
      break;
    case CONTOURING_METHOD_CGAL:
      signed_distance_isosurface(
        IV,IF,level,angle_bound,radius_bound,distance_bound,type,V,F);
      break;
  }
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 5:
    {
      prepare_lhs_double(S,plhs+4);
      // Fall through
    }
    case 4:
    {
      prepare_lhs_double(side,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_double(GV,plhs+2);
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

#ifdef MEX
#ifdef WITH_EMBREE
#  include "winding_number_ray.h"
#  include "winding_number_ray.cpp"
#endif

#include <igl/parallel_for.h>
#include <igl/winding_number.h>
#include <igl/WindingNumberAABB.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/writeOBJ.h>
#include <igl/get_seconds.h>
#include <igl/per_face_normals.h>
#ifdef WITH_EMBREE
#  include <igl/embree/EmbreeIntersector.h>
#endif
#include <igl/random_dir.h>

#include <mex.h>
#include <Eigen/Dense>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/C_STR.h>
#include <iostream>

// Inputs:
//   nrhs  number of right hand side arguments
//   prhs  pointer to right hand side arguments
// Outputs:
//   V  #V by dim list of vertex positiobns
//   F  #F by dim+1 list of facet indices
//   O  #O by dim list of origin positions
//   hierarchical  Whether to use hierarchical evaluation
//   ray_cast  Whether to use Ray tracing
//   ray_parity  Use parity version
//   num_rays  Number of rays
//   twod_rays  2d rays only
//
void parse_rhs(
  const int nrhs,
  const mxArray *prhs[],
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & O,
  bool & hierarchical,
  bool & ray_cast,
  bool & ray_parity,
  int & num_rays,
  bool & twod_rays)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3 || dim == 2,
    "Mesh vertex list must be #V by 2 or 3 list of vertex positions");

  mexErrMsgTxt(dim == mxGetN(prhs[1]),
    "Mesh \"face\" simplex size must equal dimension");

  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);
  parse_rhs_double(prhs+2,O);
  mexErrMsgTxt(mxGetN(prhs[2]) == dim, 
    "Origin list dimension must match vertex list dimension");
  // set number of faces

  // Default values
  hierarchical = true;
  ray_cast = false;
  ray_parity = false;
  num_rays = 32;
  twod_rays = false;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),
        "Parameter names should be char strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Hierarchical",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        hierarchical = (mxLogical)*mxGetPr(prhs[++i]);
        mexErrMsgTxt(!hierarchical || dim==3,
          "Hierarchical only supported for dim = 3");
      }else if(strcmp("RayCast",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        ray_cast = (mxLogical)*mxGetPr(prhs[++i]);
        mexErrMsgTxt(!ray_cast || dim==3,
          "RayCast only supported for dim = 3");
      }else if(strcmp("RayParity",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        ray_parity = (mxLogical)*mxGetPr(prhs[++i]);
        mexErrMsgTxt(!ray_parity|| dim==3,
          "RayParity  only supported for dim = 3");
      }else if(strcmp("TwoDRays",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        twod_rays = (mxLogical)*mxGetPr(prhs[++i]);
      }else if(strcmp("NumRays",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        num_rays = (int)*mxGetPr(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,
          C_STR("Unsupported parameter: "<<name));
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

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);

  Eigen::MatrixXd V,O;
  Eigen::MatrixXi F;
  bool hierarchical, ray_cast, ray_parity, twod_rays;
  int num_rays;
  parse_rhs(
    nrhs,prhs,V,F,O, hierarchical, ray_cast, ray_parity, num_rays, twod_rays);
  const int dim = V.cols();

  //// Set up openmp
  //int nProcessors=omp_get_max_threads();
  ////std::cout<<nProcessors<<std::endl;
  //omp_set_num_threads(min(nProcessors,2));

  Eigen::VectorXd W(O.rows(),1);
  double start_sec = igl::get_seconds();
  if(ray_cast)
  {
#ifdef WITH_EMBREE
    // Initialize Embree
    //cout<<"Flipping faces..."<<endl;
    MatrixXi FF;
    FF.resize(F.rows()*2,F.cols());
    FF << F, F.rowwise().reverse().eval();
    //cout<<"Computing normals..."<<endl;
    Eigen::MatrixXd N;
    per_face_normals(V,F,N);
    //cout<<"Initializing Embree..."<<endl;
    // Initialize intersector
    igl::embree::EmbreeIntersector ei;
    ei.init(V.cast<float>(),FF);

    // loop over origins
#   pragma omp parallel for if (no>IGL_WINDING_NUMBER_OMP_MIN_VALUE)
    for(int o = 0;o<no;o++)
    {
      RowVector3d p = O.row(o);
      for(int r = 0;r<num_rays;r++)
      {
        int num_rays_shot = -1;
        Vector3d dir = random_dir();
        if(twod_rays)
        {
          dir[2] = 0;
          dir.normalize();
        }
        double w = 0;
        if(ray_parity)
        {
          winding_number_ray_parity(ei,N,p,dir,w,num_rays_shot);
        }else
        {
          winding_number_ray(ei,N,p,dir,w,num_rays_shot);
        }
        W(o) += w;
      }
      W(o) /= (double)num_rays;
    }


#else
    igl::matlab::mexErrMsgTxt(false,"Recompile with WITH_EMBREE defined");
#endif
  }else
  {
    if(hierarchical && dim == 3)
    {
      // Initialize hierarchy
      WindingNumberAABB< Eigen::RowVector3d, Eigen::MatrixXd, Eigen::MatrixXi>
        hier(V,F);
      // Build hierarchy
      hier.grow();
      // loop over origins
      igl::parallel_for(
        O.rows(),
        [&hier,&W,&O](const int o)
        {
          W(o) = hier.winding_number(O.row(o));
        },
        10000);
    }else
    {
      winding_number(V,F,O,W);
    }
  }
  //mexPrintf("Elapsed time is %g seconds\n",igl::get_seconds()-start_sec);

  switch(nlhs)
  {
    case 1:
    {
      igl::matlab::prepare_lhs_double(W,plhs+0);
      // Fall through
    }
    default:break;
  }


  std::cout.rdbuf(outbuf);
}

#endif


#include <igl/doublearea.h>
#include <igl/flipped_triangles.h>
#include <igl/EPS.h>
#include <igl/slice.h>
#include <igl/euler_characteristic.h>
#include <igl/C_STR.h>
#include <igl/slim.h>
#include <igl/is_edge_manifold.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/components.h>
#ifdef MEX
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>
#else
#include <igl/matlab/MatlabWorkspace.h>
#endif

#ifdef MEX
void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  const auto mexErrMsgTxt = [](const bool v ,const char * msg)
  {
    igl::matlab::mexErrMsgTxt(v,msg);
  };
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexErrMsgTxt(nrhs >= 4,"Four arguments expected");
#else
int main(int argc, char * argv[])
{
  const auto mexErrMsgTxt = [](const bool v ,const char * msg)
  {
    if(!v)
    {
      std::cerr<<msg<<std::endl;
      exit(EXIT_FAILURE);
    }
  };
#endif
  Eigen::MatrixXd V,U0,U,bc;
  Eigen::MatrixXi F;
  Eigen::VectorXi b;
  int iters = 100;
  bool align_guess = true;
  double p = 1e5;

#ifdef MEX
  igl::matlab::parse_rhs_double(prhs+0,V);
  igl::matlab::parse_rhs_index(prhs+1,F);
  igl::matlab::parse_rhs_index(prhs+2,b);
  igl::matlab::parse_rhs_double(prhs+3,bc);

  {
    int i = 4;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("P",name) == 0)
      {
        igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
        igl::matlab::validate_arg_double(i,nrhs,prhs,name);
        p = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("AlignGuess",name) == 0)
      {
        igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
        igl::matlab::validate_arg_logical(i,nrhs,prhs,name);
        align_guess = (bool)*mxGetLogicals(prhs[++i]);
      }else if(strcmp("Iters",name) == 0)
      {
        igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
        igl::matlab::validate_arg_double(i,nrhs,prhs,name);
        iters = (double)*mxGetPr(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
      }
      i++;
    }
  }

#else
  igl::matlab::MatlabWorkspace mw;
  mw.read(argv[1]);
  mw.find("V",V);
  mw.find_index("F",F);
  mw.find_index("b",b);
  mw.find("bc",bc);
#endif

  {
    Eigen::MatrixXi C;
    igl::components(F, C);
    mexErrMsgTxt(
      C.maxCoeff() == 0,
      "(V,F) should have exactly 1 connected component");
    const int ec = igl::euler_characteristic(F);
    mexErrMsgTxt(
      ec == 1,
      C_STR("(V,F) should have disk topology (euler characteristic = "<<ec));
    mexErrMsgTxt(
      igl::is_edge_manifold(F),
      "(V,F) should be edge-manifold");
    Eigen::VectorXd A;
    igl::doublearea(V,F,A);
    mexErrMsgTxt(
      (A.array().abs() > igl::EPS<double>()).all(),
      "(V,F) should have non-zero face areas");
  }


  // Initialize with Tutte embedding to disk
  Eigen::VectorXi bnd; 
  Eigen::MatrixXd bnd_uv;
  igl::boundary_loop(F,bnd);
  igl::map_vertices_to_circle(V,bnd,bnd_uv);
  igl::harmonic(V,F,bnd,bnd_uv,1,U0);
  if (igl::flipped_triangles(U0,F).size() != 0) 
  {
    // use uniform laplacian
    igl::harmonic(F,bnd,bnd_uv,1,U0);
  }

  if(align_guess)
  {
    Eigen::MatrixXd X,X0,Y,Y0;
    igl::slice(U0,b,1,X);
    Y = bc;
    Eigen::RowVectorXd Xmean = X.colwise().mean();
    Eigen::RowVectorXd Ymean = Y.colwise().mean();
    Y0 = Y.rowwise()-Ymean;
    X0 = X.rowwise()-Xmean;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X0.cols(),X0.cols())*1e-5;
    Eigen::MatrixXd T = (X0.transpose()*X0+I).inverse()*(X0.transpose()*Y0);
    U0.rowwise() -= Xmean;
    U0 = (U0*T).eval();
    U0.rowwise() += Ymean;
  }

  if(igl::flipped_triangles(U0,F).size() == F.rows())
  {
    F = F.array().rowwise().reverse().eval();
  }

  mexErrMsgTxt(
    igl::flipped_triangles(U0,F).size() == 0,
    "Failed to initialize to feasible guess");

  igl::SLIMData slim;
  slim.energy = igl::SLIMData::SYMMETRIC_DIRICHLET;
  igl::slim_precompute(V,F,U0,slim,igl::SLIMData::SYMMETRIC_DIRICHLET,b,bc,p);
  igl::slim_solve(slim,iters);
  U = slim.V_o;

#ifdef MEX
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 2:
    {
      igl::matlab::prepare_lhs_double(U0,plhs+1);
      // Fall through
    }
    case 1:
    {
      igl::matlab::prepare_lhs_double(U,plhs+0);
      // Fall through
    }
    case 0: break;
  }
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
#else
  mw.save(U,"U");
  mw.write(argv[1]);
#endif
}

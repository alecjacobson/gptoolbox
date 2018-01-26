// Based on `collide_eltopo_mex.cpp` by Leonardo Sacht, 2014.

#include <igl/massmatrix.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <eltopo.h>
#include <iostream>

typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> MatrixX3d;
typedef Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> MatrixX3i;

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  MatrixX3d & V0,
  MatrixX3i & F,
  MatrixX3d & V1)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  using namespace igl::matlab;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const int dim = mxGetN(prhs[0]);
  mexErrMsgTxt(dim == 3,
    "Mesh vertex list must be #V0 by 3 list of vertex positions");
  mexErrMsgTxt(mxGetN(prhs[1]) == 3,
    "Mesh vertex list must be #V1 by 3 list of vertex positions");

  parse_rhs_double(prhs,V0);
  parse_rhs_index(prhs+1,F);
  parse_rhs_double(prhs+2,V1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  MatrixX3d V0,V1,U;
  MatrixX3i F,G;
  Eigen::VectorXd M;
  double out_dt = 0;
  double eps_prox = 1e-4;
  double tol_dt = 1e-1;


  // Parse rhs
  parse_rhs(nrhs,prhs,V0,F,V1);
  {
    Eigen::SparseMatrix<double> MM;
    massmatrix(V0,F,MASSMATRIX_TYPE_DEFAULT,MM);
    M = MM.diagonal();
  }

  ElTopoMesh eltopo;
  eltopo.num_vertices = V0.rows();
  eltopo.vertex_locations = V0.data();
  eltopo.num_triangles = F.rows();
  eltopo.triangles = F.data();
  eltopo.vertex_masses = M.data();
  // Set general parameters
  ElTopoGeneralOptions gen_options;
  gen_options.m_verbose = 0;
  gen_options.m_collision_safety = 1;
  gen_options.m_proximity_epsilon = eps_prox;
  // Set Simulation parameters
  ElTopoIntegrationOptions int_options;
  int_options.m_friction_coefficient = 0.0;
  int_options.m_dt = 1.0;

  // For now topology optimization is not supported
  G = F;

  // run simulation (cut time step if it does not satisfy constraints)
  {
    double* Up;
    out_dt = 0.0;
    double rest_dt = 1.0;
    int attempts = 0;
    while (rest_dt > tol_dt)
    {
    el_topo_integrate(&eltopo, V1.data(), &gen_options, &int_options, &Up, &out_dt);
      eltopo.vertex_locations = Up;
      rest_dt = (1-out_dt)*rest_dt;
      //mexPrintf("current out_dt = %.4f\n", out_dt);
      //mexPrintf("running rest_dt = %.4f\n", rest_dt);
      attempts = attempts+1;
      if (out_dt<tol_dt)
      {
        Up = V0.data();
        break;
      }
    }
    U.resize(eltopo.num_vertices,3);
    copy(Up,Up+U.size(),U.data());
  }

  switch(nlhs)
  {
    default:
      mexErrMsgTxt(false,"Too many output parameters.");
      // fall through
    case 3:
      plhs[2] = mxCreateDoubleScalar(out_dt);
      // fall through
    case 2:
      prepare_lhs_index(G,plhs+1);
      // fall through
    case 1:
      prepare_lhs_double(U,plhs+0);
      // fall through
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

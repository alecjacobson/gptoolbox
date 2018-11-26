#include <igl/embree/bone_visible.h>

#include <igl/matlab/MexStream.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <functional>

// Parse right hand side arguments for a matlab mex function.
//
// Inputs:
//   nrhs  number of right hand side arguments
//   prhs  pointer to right hand side arguments
// Outputs:
//   V  n by dim list of mesh vertex positions
//   F  m by dim list of mesh face indices
//   s  1 by dim bone source vertex position
//   d  1 by dim bone dest vertex position
// "Throws" matlab errors if dimensions are not sane.
void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::VectorXd & s,
  Eigen::VectorXd & d)
{
  using namespace std;
  if(nrhs < 4)
  {
    mexErrMsgTxt("nrhs < 4");
  }
  const int dim = mxGetN(prhs[0]);
  if(dim != 3)
  {
    mexErrMsgTxt("Mesh vertex list must be #V by 3 list of vertex positions");
  }
  if(dim != (int)mxGetN(prhs[1]))
  {
   mexErrMsgTxt("Mesh facet size must equal dimension");
  }
  if(dim != (int)mxGetN(prhs[2]))
  {
   mexErrMsgTxt("Source dim must equal vertex dimension");
  }
  if(dim != (int)mxGetN(prhs[3]))
  {
   mexErrMsgTxt("Dest dim must equal vertex dimension");
  }
  // set number of mesh vertices
  const int n = mxGetM(prhs[0]);
  // set vertex position pointers
  double * Vp = mxGetPr(prhs[0]);
  // set number of faces
  const int m = mxGetM(prhs[1]);
  // set face index list pointer
  double * Fp = mxGetPr(prhs[1]);
  // set source and dest pointers
  double * sp = mxGetPr(prhs[2]);
  double * dp = mxGetPr(prhs[3]);
  // resize output to transpose
  V.resize(n,dim);
  copy(Vp,Vp+n*dim,V.data());
  // resize output to transpose
  F.resize(m,dim);
  // Q: Is this doing a cast?
  // A: Yes.
  copy(Fp,Fp+m*dim,F.data());
  // http://stackoverflow.com/a/4461466/148668
  transform(F.data(),F.data()+m*dim,F.data(),
    bind2nd(std::plus<double>(),-1.0));
  // resize output to transpose
  s.resize(dim);
  copy(sp,sp+dim,s.data());
  d.resize(dim);
  copy(dp,dp+dim,d.data());
}

// Prepare left hand side
//
// Inputs:
//   nlhs  number of left hand side arguments
//   plhs  pointer to left hand side arguments
//   no  number of origins
// Outputs:
//   S  no by 1  list of solid angles
void prepare_lhs(
  const int nlhs,
  mxArray *plhs[],
  const int no,
  double *& S)
{
  // Create pointers for output arrays
  if(nlhs > 0)
  {
    // solid angles
    plhs[0] = mxCreateDoubleMatrix(no,1, mxREAL);
    S = mxGetPr(plhs[0]);
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
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd V;
  MatrixXi F;
  VectorXd s;
  VectorXd d;
  parse_rhs(nrhs,prhs,V,F,s,d);

  //// Set up openmp
  //int nProcessors=omp_get_max_threads();
  ////std::cout<<nProcessors<<std::endl;
  //omp_set_num_threads(min(nProcessors,2));

  // List of flag
  Matrix<bool,Dynamic,1> flag;
  igl::embree::bone_visible(V,F,s,d,flag);

  nlhs = 1;
  plhs[0] = mxCreateLogicalMatrix(V.rows(),1);
  mxLogical * flagp = mxGetLogicals(plhs[0]);
  //copy(&flag.data()[0],&flag.data()[0]+flag.size(),flagp);
  // safe copy
  for(int f = 0;f<flag.size();f++)
  {
    flagp[f] = flag(f);
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

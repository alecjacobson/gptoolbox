#include <mex.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/simplex_simplex_squared_distance.h>
#include <igl/C_STR.h>
#include <Eigen/Core>

// Call igl::simplex_simplex_squared_distance with the simplex sizes and
// dimension baked in at compile time. That lets libigl's recursion run on
// stack-allocated fixed-size blocks, which is where most of the speedup over
// the pure MATLAB implementation comes from.
template <int N1, int N2, int Dim>
static void ssd(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  double & sqrd,
  Eigen::RowVectorXd & B1,
  Eigen::RowVectorXd & B2)
{
  const Eigen::Matrix<double,N1,Dim> A = V1;
  const Eigen::Matrix<double,N2,Dim> B = V2;
  Eigen::Matrix<double,1,N1> W1;
  Eigen::Matrix<double,1,N2> W2;
  igl::simplex_simplex_squared_distance(A,B,sqrd,W1,W2);
  B1 = W1;
  B2 = W2;
}

// Simplices with more corners than this fall back to the dynamically sized
// instantiation. 4 covers points, segments, triangles and tets.
#define GPTOOLBOX_SSD_MAX_STATIC_CORNERS 4

template <int Dim>
static void ssd_dispatch(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  double & sqrd,
  Eigen::RowVectorXd & B1,
  Eigen::RowVectorXd & B2)
{
  const int n1 = (int)V1.rows();
  const int n2 = (int)V2.rows();
#define GPTOOLBOX_SSD_CASE(M1,M2) \
  case (M1)*(GPTOOLBOX_SSD_MAX_STATIC_CORNERS+1)+(M2): \
    return ssd<M1,M2,Dim>(V1,V2,sqrd,B1,B2);
  if(
    n1>=1 && n1<=GPTOOLBOX_SSD_MAX_STATIC_CORNERS &&
    n2>=1 && n2<=GPTOOLBOX_SSD_MAX_STATIC_CORNERS)
  {
    switch(n1*(GPTOOLBOX_SSD_MAX_STATIC_CORNERS+1)+n2)
    {
      GPTOOLBOX_SSD_CASE(1,1) GPTOOLBOX_SSD_CASE(1,2)
      GPTOOLBOX_SSD_CASE(1,3) GPTOOLBOX_SSD_CASE(1,4)
      GPTOOLBOX_SSD_CASE(2,1) GPTOOLBOX_SSD_CASE(2,2)
      GPTOOLBOX_SSD_CASE(2,3) GPTOOLBOX_SSD_CASE(2,4)
      GPTOOLBOX_SSD_CASE(3,1) GPTOOLBOX_SSD_CASE(3,2)
      GPTOOLBOX_SSD_CASE(3,3) GPTOOLBOX_SSD_CASE(3,4)
      GPTOOLBOX_SSD_CASE(4,1) GPTOOLBOX_SSD_CASE(4,2)
      GPTOOLBOX_SSD_CASE(4,3) GPTOOLBOX_SSD_CASE(4,4)
      default: break;
    }
  }
#undef GPTOOLBOX_SSD_CASE
  return ssd<Eigen::Dynamic,Eigen::Dynamic,Dim>(V1,V2,sqrd,B1,B2);
}

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[]
  )
{
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs==2,"nrhs should be == 2");
  mexErrMsgTxt(nlhs<=5,"nlhs should be <= 5");

  MatrixXd V1,V2;
  parse_rhs_double(prhs+0,V1);
  parse_rhs_double(prhs+1,V2);

  mexErrMsgTxt(V1.rows()>0,"V1 must not be empty");
  mexErrMsgTxt(V2.rows()>0,"V2 must not be empty");
  mexErrMsgTxt(V1.cols()==V2.cols(),
    "V1 and V2 must have the same number of columns");
  // igl::simplex_simplex_squared_distance indexes faces with a 64-bit mask.
  mexErrMsgTxt(V1.rows()<=62,"V1 has too many corners");
  mexErrMsgTxt(V2.rows()<=62,"V2 has too many corners");

  double sqrd;
  RowVectorXd B1,B2;
  switch((int)V1.cols())
  {
    case 2: ssd_dispatch<2>(V1,V2,sqrd,B1,B2); break;
    case 3: ssd_dispatch<3>(V1,V2,sqrd,B1,B2); break;
    default: ssd_dispatch<Eigen::Dynamic>(V1,V2,sqrd,B1,B2); break;
  }

  if(nlhs>=5){ prepare_lhs_double(B2,plhs+4); }
  if(nlhs>=4){ prepare_lhs_double(B1,plhs+3); }
  // Closest points, recovered from the barycentric coordinates.
  if(nlhs>=3){ prepare_lhs_double((B2*V2).eval(),plhs+2); }
  if(nlhs>=2){ prepare_lhs_double((B1*V1).eval(),plhs+1); }
  prepare_lhs_double(Eigen::Matrix<double,1,1>(sqrd),plhs+0);

  std::cout.rdbuf(outbuf);
}

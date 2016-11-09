#include <igl/remove_unreferenced.h>
#include <igl/unique.h>
#include <igl/copyleft/cgal/snap_rounding.h>
#include <igl/copyleft/cgal/resolve_intersections.h>
#include <igl/copyleft/cgal/subdivide_segments.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <algorithm>

// mex(  MEXOPTS{:}, MSSE42, STDCPP11,  LIBIGL_BASE{:},EIGEN_INC, CGAL_INC,CGAL_LIB{:},CGAL_FLAGS,  LIBIGL_LIBCGAL, LIBIGL_LIBBOOLEAN, CORK_INC,CORK_LIB{:},  BOOST_INC,BOOST_LIB{:},  'mesh_boolean.cpp');

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  Eigen::MatrixXd V,VI;
  Eigen::MatrixXi E,EI;
  Eigen::VectorXi I;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);

  mexErrMsgTxt(mxGetN(prhs[0]) == 0 || mxGetN(prhs[0]) == 2,
      "Mesh vertex list must be #V by 3 list of vertex positions");
  mexErrMsgTxt(mxGetN(prhs[1])==0 || 2 == mxGetN(prhs[1]),
      "E must segments");
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,E);

  igl::copyleft::cgal::snap_rounding(V,E,VI,EI,I);


  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_index(I,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(EI,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(VI,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::copyleft::cgal::snap_rounding<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif

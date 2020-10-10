// Force header only version
#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif
#include <mex.h>

#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/STR.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/prepare_lhs.h>

#include <igl/xml/serialize_xml.h>
#include <Eigen/Core>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <wordexp.h>
#endif
#include <iostream>

// http://www.alecjacobson.com/weblog/?p=4477
typedef CGAL::Epeck::FT EScalar;
typedef Eigen::Matrix<EScalar,Eigen::Dynamic,Eigen::Dynamic> MatrixXE;
namespace igl
{
  namespace xml
  {
    namespace serialization_xml
    {
      template <> inline void serialize(
        const MatrixXE & obj,
        tinyxml2::XMLDocument* doc,
        tinyxml2::XMLElement* element,
        const std::string& name)
      {
        const std::function<std::string(const MatrixXE::Scalar &) > to_string = 
          [](const MatrixXE::Scalar & v)->std::string
          {
            return
              STR(CGAL::exact(v));
          };
        serialize(obj,name,to_string,doc,element);
      }
      template <> inline void deserialize(
        MatrixXE & obj,
        const tinyxml2::XMLDocument* doc,
        const tinyxml2::XMLElement* element,
        const std::string& name)
      {
        const std::function<void(const std::string &,MatrixXE::Scalar &)> & 
          from_string = 
          [](const std::string & s, MatrixXE::Scalar & v)
          {
            std::stringstream(s)>>v;
          };
        deserialize(doc,element,name,from_string,obj);
      }
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
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs >= 1, "The number of input arguments must be >=1.");
  mexErrMsgTxt(mxIsChar(prhs[0]),"File name should be string");
  string filename;
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
  {
    wordexp_t exp_result;
    wordexp(mxArrayToString(prhs[0]), &exp_result, 0);
    filename = exp_result.we_wordv[0];
  }
#endif



  MatrixXE V;
  MatrixXi F;

  // Read mesh
  igl::xml::deserialize_xml(V,"vertices",filename);
  igl::xml::deserialize_xml(F,"faces",    filename);

  plhs[0] = mxCreateDoubleMatrix(V.rows(),V.cols(), mxREAL);
  double * Vp = mxGetPr(plhs[0]);
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 2:
    {
      prepare_lhs_index(F,plhs+1);
      // fall through
    }
    case 1:
    {
      const int m = V.rows();
      const int n = V.cols();
      plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
      double * Vp = mxGetPr(plhs[0]);
      for(int i = 0;i<m;i++)
      {
        for(int j = 0;j<n;j++)
        {
          Vp[i+m*j] = CGAL::to_double(V(i,j));
        }
      }
      // fall through
    }
    case 0: break;
  }

}

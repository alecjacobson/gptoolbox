//     mex -v -DMEX intersect_other.cpp ...
//       -Isrc -I/opt/local/include/ -I/opt/local/include/eigen3 ...
//       -I/usr/local/igl/libigl/include -L/opt/local/lib ...
//       -lCGAL -lCGAL_Core -lgmp -lmpfr ...
//       -lboost_thread-mt -lboost_system-mt
// 
#ifdef MEX
#  include <mex.h>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

#include <igl/read_triangle_mesh.h>
#include <igl/pathinfo.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/doublearea.h>
#include <igl/writeDMAT.h>
#include <igl/unique_simplices.h>
#include <igl/C_STR.h>
#ifdef MEX
#  define IGL_REDRUM_NOOP
#endif
#include <igl/REDRUM.h>
#ifdef MEX
#  include <igl/matlab/MexStream.h>
#  include <igl/matlab/mexErrMsgTxt.h>
#  include <igl/matlab/parse_rhs.h>
#endif
#include <igl/cgal/intersect_other.h>

#ifdef MEX
#  include "mex.h"
#endif

#include <iostream>
#include <string>

#ifdef MEX
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

#else
int main(int argc, char * argv[])
{
#endif
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::cgal;

  MatrixXd V,U;
  MatrixXi F,G;
  bool first_only = false;

  string prefix;
  bool use_obj_format = false;
#ifdef MEX
  const int NUM_REQ = 4;
  if(nrhs < NUM_REQ)
  {
    mexErrMsgTxt(C_STR("nrhs < "<<NUM_REQ));
  }
  const auto & parse_mesh = [](
      const mxArray *prhs[],
      MatrixXd & V,
      MatrixXi & F)
  {
    parse_rhs_double(prhs,V);
    parse_rhs_index(prhs+1,F);
    mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
    mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  };
  parse_mesh(prhs,V,F);
  parse_mesh(prhs+2,U,G);

  if(nrhs>NUM_REQ)
  {
    int i = NUM_REQ;
    while(i<nrhs)
    {
      if(!mxIsChar(prhs[i]))
      {
        mexErrMsgTxt("Parameter names should be char strings");
      }
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("FirstOnly",name) == 0)
      {
        if((i+1)>=nrhs)
        {
          mexErrMsgTxt(C_STR("Parameter '"<<name<<"' requires argument"));
        }
        i++;
        if(!mxIsLogical(prhs[i]))
        {
          mexErrMsgTxt(C_STR("Parameter '"<<name<<"' requires Logical arg"));
        }
        mxLogical * v = (mxLogical *)mxGetData(prhs[i]);
        first_only = *v;
      }else
      {
        mexErrMsgTxt(C_STR("Unsupported parameter: "<<name));
      }
      i++;
    }
  }
#else
  if(argc <= 2)
  {
    cerr<<"Usage:"<<endl<<"  intersect_other [path to .off/.obj mesh A] [path to .off/.obj mesh B] "
      "[0 or 1 for first only]"<<endl;
    return 1;
  }
  // Apparently CGAL doesn't have a good data structure triangle soup. Their
  // own examples use (V,F):
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/AABB_tree/Chapter_main.html#Subsection_64.3.7
  //
  // Load mesh
  string filenameA(argv[1]);
  if(!read_triangle_mesh(filenameA,V,F))
  {
    //cout<<REDRUM("Reading "<<filename<<" failed.")<<endl;
    return false;
  }
  cout<<GREENGIN("Read "<<filenameA<<" successfully.")<<endl;

  string filenameB(argv[2]);
  if(!read_triangle_mesh(filenameB,U,G))
  {
    //cout<<REDRUM("Reading "<<filename<<" failed.")<<endl;
    return false;
  }
  cout<<GREENGIN("Read "<<filenameB<<" successfully.")<<endl;

  {
    // dirname, basename, extension and filename
    string dirname,b,ext;
    string prefixA,prefixB;
    pathinfo(filenameB,dirname,b,ext,prefixB);
    pathinfo(filenameA,dirname,b,ext,prefixA);
    prefix = dirname + "/" + prefixA + "-" + prefixB;
    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    use_obj_format = ext == "obj";
  }

  if(argc>3)
  {
    //http://stackoverflow.com/a/9748431/148668
    char *p;
    long d = strtol(argv[3], &p, 10);
    if (errno != 0 || *p != '\0')
    {
      cerr<<"first only param should be 0 or 1"<<endl;
    }else
    {
      first_only = d!=0;
    }
  }
#endif
  const auto validate = [](const MatrixXd & V, const MatrixXi & F) -> bool
  {
    // Check that there aren't any combinatorially or geometrically degenerate triangles
    VectorXd A;
    doublearea(V,F,A);
    if(A.minCoeff()<=0)
    {
#ifdef MEX
      mexErrMsgTxt("Geometrically degenerate face found.");
#else
      cerr<<"Geometrically degenerate face found."<<endl;
      return false;
#endif
    }
    VectorXi F12,F23,F31;
    F12 = F.col(0)-F.col(1);
    F23 = F.col(1)-F.col(2);
    F31 = F.col(2)-F.col(0);
    if(
      F12.minCoeff() == 0 || 
      F23.minCoeff() == 0 || 
      F31.minCoeff() == 0)
    {
#ifdef MEX
      mexErrMsgTxt("Combinatorially degenerate face found.");
#else
      cerr<<"Geometrically degenerate face found."<<endl;
      return false;
#endif
    }
    return true;
  };

  if(!validate(V,F) || !validate(U,G))
  {
#ifndef MEX
    // Otherwise should have called mexErr
    return 1;
#endif
  }

  // Now mesh self intersections
  MatrixXi IF;
  {
    intersect_other(V,F,U,G,first_only,IF);
#ifndef MEX
    cout<<"writing pair list to "<<(prefix+"-IF.dmat")<<endl;
    writeDMAT((prefix+"-IF.dmat").c_str(),IF);
#endif
  }

#ifdef MEX
  // Prepare left-hand side
  nlhs = 1;

  // Treat indices as reals
  plhs[0] = mxCreateDoubleMatrix(IF.rows(),IF.cols(), mxREAL);
  double * IFp = mxGetPr(plhs[0]);
  MatrixXd IFd = (IF.cast<double>().array()+1).matrix();
  copy(&IFd.data()[0],&IFd.data()[0]+IFd.size(),IFp);

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);

#else
  return 0;
#endif
}

//   = Compile = 
//   I've noticed that the default gcc/g++ can be a lot of trouble with the mex,
//   cgal, boost constellation. On Mac OS X I recommend using gcc-4.7. I have had
//   bad experiences with gcc-4.3, resulting in bogus runtime errors and
//   assertions.
// 
//     = mexopts.sh =
//     Setting up your mexopts properly on a Mac depends on which version of XCode
//     you have installed. 
// 
//     Open up mexopts.sh. First you need to point the SDKROOT variable to your
//     sdk. This is usually either something like:
// 
//     SDKROOT='/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/'
// 
//     or 
// 
//     SDKROOT='/Developer/SDKs/MacOSX10.7.sdk'
// 
//     Next, you need to remove the -no-cpp-precomp flag and add the
//     -frounding-math flag to your CFLAGS and CXXFLAGS
// 
//     I also found it best to change CC and CXX to point to a g++-4.7 (or
//     g++-mp-4.7), but clang should also work.
// 
//   To mex on a mac, open up matlab and issue:
// 
//     % Alec's Macbook Air, Intel Core i7, 10.7.5, MATLAB 2012a, XCode 4.3.3,
//     % Boost 1.52.0_1, CGAL 4.1, g++-mp-4.7
//     mex -v -DMEX selfintersect.cpp ...
//       -Isrc -I/opt/local/include/ -I/opt/local/include/eigen3 ...
//       -I/usr/local/igl/libigl/include -L/opt/local/lib ...
//       -lCGAL -lCGAL_Core -lgmp -lmpfr ...
//       -lboost_thread-mt -lboost_system-mt
// 
//   If you have a non-standard installation of eigen3 and libigl, you will need
//   to be sure that in the mex commands above you change:
// 
//   /usr/local/igl/libigl ---> path to libigl directory 
//   /opt/local/include/eigen3---> path to eigen3 directory 
//   /opt/local/include/ ---> path to CGAL/Boost/gmp/mpfr include directory 
//   /opt/local/lib/ ---> path to CGAL/Boost/gmp/mpfr lib directory 
// 
//     = libboost = 
//     It seems MATLAB 2012b's libboost does not play nicely with macports
//     libboost-mt. Should compile with matlab's libboost directly: -lboost_*
//     rather than -lboost_*-mt
// 
//     % Alec's iMac, Intel Core i7, 10.7.5, MATLAB 2012b, XCode 4.2, Boost
//     % 1.52.0_1, CGAL 4.1, g++-mp-4.7
// 
//   = Run =
//   You can then run these with the following examples:
// 
//   V1 = [0 0 0;1 0 0;0 1 0;0.5 0.5 0.5;0 0 -0.5;1 1 -0.5];
//   F1 = [1 2 3;4 5 6];
//   [SV,SF] = selfintersect(V1,F1);
// 
#ifdef MEX
#  include <mex.h>
#  include <igl/C_STR.h>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#else
#include <igl/read_triangle_mesh.h>
#include <igl/pathinfo.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#endif
#include <igl/doublearea.h>
#include <igl/unique_simplices.h>
#ifdef MEX
#  define IGL_REDRUM_NOOP
#endif
#include <igl/REDRUM.h>
#ifdef MEX
#  include <igl/matlab/MexStream.h>
#  include <igl/matlab/mexErrMsgTxt.h>
#  include <igl/matlab/parse_rhs.h>
#endif
#include <igl/cgal/remesh_self_intersections.h>

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

  MatrixXd V;
  MatrixXi F;
  igl::cgal::RemeshSelfIntersectionsParam params;

  string prefix;
  bool use_obj_format = false;
#ifdef MEX
  if(nrhs < 2)
  {
    mexErrMsgTxt("nrhs < 2");
  }
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,F);
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  if(nrhs>2)
  {
    int i = 2;
    while(i<nrhs)
    {
      if(!mxIsChar(prhs[i]))
      {
        mexErrMsgTxt("Parameter names should be char strings");
      }
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("DetectOnly",name) == 0)
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
        params.detect_only = *v;
      }else if(strcmp("FirstOnly",name) == 0)
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
        params.first_only = *v;
      }else
      {
        mexErrMsgTxt(C_STR("Unsupported parameter: "<<name));
      }
      i++;
    }
  }
#else
  if(argc <= 1)
  {
    cerr<<"Usage:"<<endl<<"  selfintersect [path to .off/.obj mesh] "
      "[0 or 1 for detect only]"<<endl;
    return 1;
  }
  // Apparently CGAL doesn't have a good data structure triangle soup. Their
  // own examples use (V,F):
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/AABB_tree/Chapter_main.html#Subsection_64.3.7
  //
  // Load mesh
  string filename(argv[1]);
  if(!read_triangle_mesh(filename,V,F))
  {
    //cout<<REDRUM("Reading "<<filename<<" failed.")<<endl;
    return false;
  }
  cout<<GREENGIN("Read "<<filename<<" successfully.")<<endl;
  {
    // dirname, basename, extension and filename
    string dirname,b,ext;
    pathinfo(filename,dirname,b,ext,prefix);
    prefix = dirname + "/" + prefix;
    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    use_obj_format = ext == "obj";
  }
  if(argc>2)
  {
    //http://stackoverflow.com/a/9748431/148668
    char *p;
    long d = strtol(argv[2], &p, 10);
    if (errno != 0 || *p != '\0')
    {
      cerr<<"detect only param should be 0 or 1"<<endl;
    }else
    {
      params.detect_only = d!=0;
    }
  }
#endif
  MatrixXi IF;
  VectorXi J,IM;
  if(F.rows()>0)
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
      return 1;
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
      return 1;
#endif
    }

    // Now mesh self intersections
    {
      MatrixXd tempV;
      MatrixXi tempF;
      remesh_self_intersections(V,F,params,tempV,tempF,IF,J,IM);
      //cout<<BLUEGIN("Found and meshed "<<IF.rows()<<" pair"<<(IF.rows()==1?"":"s")
      //  <<" of self-intersecting triangles.")<<endl;
      V=tempV;
      F=tempF;
#ifndef MEX
      cout<<"writing pair list to "<<(prefix+"-IF.dmat")<<endl;
      writeDMAT((prefix+"-IF.dmat").c_str(),IF);
      cout<<"writing map to F list to "<<(prefix+"-J.dmat")<<endl;
      writeDMAT((prefix+"-J.dmat").c_str(),J);
      cout<<"writing duplicat index map to "<<(prefix+"-IM.dmat")<<endl;
      writeDMAT((prefix+"-IM.dmat").c_str(),IM);
      if(!params.detect_only)
      {
        if(use_obj_format)
        {
          cout<<"writing mesh to "<<(prefix+"-selfintersect.obj")<<endl;
          writeOBJ(prefix+"-selfintersect.obj",V,F);
        }else
        {
          cout<<"writing mesh to "<<(prefix+"-selfintersect.off")<<endl;
          writeOFF(prefix+"-selfintersect.off",V,F);
        }
      }
#endif
    }

  // Double-check output

#ifdef DEBUG
    // There should be *no* combinatorial duplicates
    {
      MatrixXi tempF;
      unique_simplices(F,tempF);
      if(tempF.rows() < F.rows())
      {
        cout<<REDRUM("Error: selfintersect created "<<
          F.rows()-tempF.rows()<<" combinatorially duplicate faces")<<endl;
      }else
      {
        assert(tempF.rows() == F.rows());
        cout<<GREENGIN("selfintersect created no duplicate faces")<<endl;
      }
      F = tempF;
    }
#endif
  }

#ifdef MEX
  // Prepare left-hand side
  switch(nlhs)
  {
    case 5:
    {
      // Treat indices as reals
      plhs[4] = mxCreateDoubleMatrix(IM.rows(),IM.cols(), mxREAL);
      double * IMp = mxGetPr(plhs[4]);
      VectorXd IMd = (IM.cast<double>().array()+1).matrix();
      copy(&IMd.data()[0],&IMd.data()[0]+IMd.size(),IMp);
      // Fallthrough
    }
    case 4:
    {
      // Treat indices as reals
      plhs[3] = mxCreateDoubleMatrix(J.rows(),J.cols(), mxREAL);
      double * Jp = mxGetPr(plhs[3]);
      VectorXd Jd = (J.cast<double>().array()+1).matrix();
      copy(&Jd.data()[0],&Jd.data()[0]+Jd.size(),Jp);
      // Fallthrough
    }
    case 3:
    {
      // Treat indices as reals
      plhs[2] = mxCreateDoubleMatrix(IF.rows(),IF.cols(), mxREAL);
      double * IFp = mxGetPr(plhs[2]);
      MatrixXd IFd = (IF.cast<double>().array()+1).matrix();
      copy(&IFd.data()[0],&IFd.data()[0]+IFd.size(),IFp);
      // Fallthrough
    }
    case 2:
    {
      // Treat indices as reals
      plhs[1] = mxCreateDoubleMatrix(F.rows(),F.cols(), mxREAL);
      double * Fp = mxGetPr(plhs[1]);
      MatrixXd Fd = (F.cast<double>().array()+1).matrix();
      copy(&Fd.data()[0],&Fd.data()[0]+Fd.size(),Fp);
      // Fallthrough
    }
    case 1:
    {
      plhs[0] = mxCreateDoubleMatrix(V.rows(),V.cols(), mxREAL);
      double * Vp = mxGetPr(plhs[0]);
      copy(&V.data()[0],&V.data()[0]+V.size(),Vp);
      break;
    }
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);

#else
  return 0;
#endif
}

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
#  include <igl/matlab/validate_arg.h>
#  include <igl/matlab/parse_rhs.h>
#endif
#include <igl/copyleft/cgal/remesh_self_intersections.h>

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
  using namespace igl::copyleft::cgal;

  MatrixXd V;
  MatrixXi F;
  igl::copyleft::cgal::RemeshSelfIntersectionsParam params;

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
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_logical(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        params.detect_only = *v;
      }else if(strcmp("FirstOnly",name) == 0)
      {
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_logical(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        params.first_only = *v;
      }else if(strcmp("StitchAll",name) == 0)
      {
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_logical(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        params.stitch_all = *v;
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

#include <mex.h>
#include <igl/C_STR.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/read_triangle_mesh.h>
#include <igl/pathinfo.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/doublearea.h>
#include <igl/unique_simplices.h>
#define IGL_REDRUM_NOOP
#include <igl/REDRUM.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>

#include "mex.h"

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

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
  MatrixXi IF;
  VectorXi J,IM;
  if(F.rows()>0)
  {
    // Check that there aren't any combinatorially or geometrically degenerate triangles
    VectorXd A;
    doublearea(V,F,A);
    if(A.minCoeff()<=0)
    {
      mexErrMsgTxt("Geometrically degenerate face found.");
    }
    if(
       (F.array().col(0) == F.array().col(1)).any() ||
       (F.array().col(1) == F.array().col(2)).any() ||
       (F.array().col(2) == F.array().col(0)).any())
    {
      mexErrMsgTxt("Combinatorially degenerate face found.");
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
}

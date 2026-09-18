#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab_format.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/copyleft/tetgen/tetgenio_to_tetmesh.h>
#include <igl/boundary_facets.h>
#ifndef TETLIBRARY
#define TETLIBRARY
#endif
#include <tetgen.h>
#include <Eigen/Core>
#include <cstring>

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs+0,TV);
  parse_rhs_index(prhs+1,TT);
  mexErrMsgTxt(TV.cols()==3,"TV should be #TV by 3");
  mexErrMsgTxt(TT.cols()==4,"TT should be #TT by 4");

  std::string flags = "";
  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Flags",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        flags = mxArrayToString(prhs[++i]);
      }else if(strcmp("Boundary",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        parse_rhs_index(prhs+i+1,TF);
        i++;
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  // Reconstruct a previous tetrahedralization and refine it. Tetgen's
  // post-refinement vertex smoothing pass ('-s', on by default whenever
  // Steiner points are inserted) is known to segfault (inside
  // lawsonflip3d/move_vertex) on some degenerate reconstructed meshes;
  // disable it (smooth_maxiter=0) unless the caller already asked for a
  // specific smoothing behavior via Flags.
  std::string full_flags = flags + "rQ";
  if(flags.find('s')==std::string::npos)
  {
    full_flags += "s0/0";
  }

  tetgenio in,out;
  in.firstnumber = 0;

  in.numberofpoints = TV.rows();
  in.pointlist = new REAL[TV.rows()*3];
  for(int i = 0;i<TV.rows();i++)
  {
    for(int c = 0;c<3;c++)
    {
      in.pointlist[3*i+c] = TV(i,c);
    }
  }

  in.numberoftetrahedra = TT.rows();
  in.numberofcorners = 4;
  in.tetrahedronlist = new int[TT.rows()*4];
  for(int i = 0;i<TT.rows();i++)
  {
    for(int c = 0;c<4;c++)
    {
      in.tetrahedronlist[4*i+c] = TT(i,c);
    }
  }

  if(TF.rows()==0)
  {
    // tetgen's "-r" reconstruction is more robust when told explicitly
    // which trifaces are on the boundary, rather than left to infer it.
    igl::boundary_facets(TT,TF);
  }
  mexErrMsgTxt(TF.cols()==3,"Boundary should be #TF by 3");
  in.numberoftrifaces = TF.rows();
  in.trifacelist = new int[TF.rows()*3];
  for(int i = 0;i<TF.rows();i++)
  {
    for(int c = 0;c<3;c++)
    {
      in.trifacelist[3*i+c] = TF(i,c);
    }
  }

  char * cflags = new char[full_flags.size()+1];
  std::strcpy(cflags,full_flags.c_str());
  try
  {
    ::tetrahedralize(cflags,&in,&out);
  }catch(int e)
  {
    delete[] cflags;
    std::string reason;
    switch(e)
    {
      case 1: reason = "out of memory"; break;
      case 3: reason = "input contains self-intersections"; break;
      case 4: reason = "input has a very small feature size (try -T to set a coarser tolerance)"; break;
      case 5: reason = "input has two very close facets (try -Y)"; break;
      case 10: reason = "invalid input"; break;
      default: reason = "internal tetgen error"; break;
    }
    ::mexErrMsgTxt(("tetgen failed to refine tetrahedralization: "+reason+" (code "+std::to_string(e)+").").c_str());
  }
  delete[] cflags;
  mexErrMsgTxt(out.numberoftetrahedra>0,"Tetgen failed to refine tetrahedralization.");

  Eigen::MatrixXd RV;
  Eigen::MatrixXi RT,RF;
  Eigen::VectorXi RM,RR,PT;
  Eigen::MatrixXi FT,RN;
  int num_regions;
  bool ret = igl::copyleft::tetgen::tetgenio_to_tetmesh(
    out,RV,RT,RF,RM,RR,RN,PT,FT,num_regions);
  mexErrMsgTxt(ret,"Failed to convert tetgen output.");
  if(RF.rows()==0)
  {
    igl::boundary_facets(RT,RF);
  }

  switch(nlhs)
  {
    case 3:
      prepare_lhs_index(RF,plhs+2);
    case 2:
      prepare_lhs_index(RT,plhs+1);
    case 1:
      prepare_lhs_double(RV,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

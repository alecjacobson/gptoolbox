#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/list_to_matrix.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <iostream>
#include <vector>
#include <algorithm>

template <int dim>
void helper(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd A1,A2,B1,B2;
  MatrixXi I;
  parse_rhs_double(prhs+0,A1);
  parse_rhs_double(prhs+1,A2);
  mexErrMsgTxt(A1.rows()==A2.rows(),"A2 must be #A by dim");
  mexErrMsgTxt(A1.cols()==dim,"A2 must be #A by dim");
  mexErrMsgTxt(A2.cols()==dim,"A2 must be #A by dim");

  typedef std::vector<int>::iterator ListIterator;
  typedef 
    CGAL::Box_intersection_d::Box_with_handle_d<double,dim, ListIterator > 
    Box;

  const auto box_up = [](
      const Eigen::MatrixXd & A1, 
      const Eigen::MatrixXd & A2,
      std::vector<Box> & boxes,
      std::vector<int> & list)
  {
    boxes.clear();
    boxes.reserve(A1.rows());
    list.resize(A1.rows());
    std::iota(list.begin(), list.end(), 0);
    for(int i = 0;i<A1.rows();i++)
    {
      Box box;
      std::vector<double> min_point(dim);
      std::vector<double> max_point(dim);
      for(int c = 0;c<dim;c++)
      {
        min_point[c] = A1(i,c);
        max_point[c] = A2(i,c);
      }
      boxes.emplace_back(
        &min_point[0],
        &max_point[0],
        list.begin() + i);
    }
  };
  std::vector<Box> Aboxes;
  std::vector<int> Alist;
  box_up(A1,A2,Aboxes,Alist);


  std::vector<std::vector<int> > vI;
  const auto cb = [&](const Box &a, const Box &b) -> void
  {
    vI.push_back({
        *a.handle(),
        *b.handle()
        });
  };
  if(nrhs == 2)
  {
    CGAL::box_self_intersection_d(Aboxes.begin(),Aboxes.end(),cb);
  }
  else
  {
    parse_rhs_double(prhs+2,B1);
    parse_rhs_double(prhs+3,B2);
    mexErrMsgTxt(B1.rows()==B2.rows(),"B2 must be #B by dim");
    mexErrMsgTxt(B1.cols()==dim,"B1 must be #B by dim");
    mexErrMsgTxt(B2.cols()==dim,"B2 must be #B by dim");
    std::vector<Box> Bboxes;
    std::vector<int> Blist;
    box_up(B1,B2,Bboxes,Blist);
    CGAL::box_intersection_d(
      Aboxes.begin(),Aboxes.end(),
      Bboxes.begin(),Bboxes.end(),
      cb);
  }
  igl::list_to_matrix(vI,I);
  switch(nlhs)
  {
    case 1:
      prepare_lhs_index(I,plhs+0);
    default:break;
  }
}

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=4 || nrhs==2,"nrhs should be == 2 or >= 4");
  const int dim = mxGetN(prhs[0]);
  switch(dim)
  {
    case 4:
      return helper<4>(nlhs,plhs,nrhs,prhs);
    case 3:
      return helper<3>(nlhs,plhs,nrhs,prhs);
    case 2:
      return helper<2>(nlhs,plhs,nrhs,prhs);
    case 1:
      return helper<1>(nlhs,plhs,nrhs,prhs);
    default:
      mexErrMsgTxt("Unsupported dimension.");
      break;
  }


  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

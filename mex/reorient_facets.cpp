#include <mex.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/embree/reorient_facets_raycast.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  using namespace igl::matlab;
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::Matrix<bool,Eigen::Dynamic,1> I;
  mexErrMsgTxt(nrhs >= 2, "Must provide V and F");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);
  mexErrMsgTxt(V.size() == 0 || V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.size() == 0 || F.cols()==3,"F must be #F by 3");


  double rays_total = F.rows()*100;
  double rays_minimum = 10;
  bool facet_wise = false;
  bool use_parity = false;
  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("NumRays",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        rays_total = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("MinRays",name) == 0)
      {
        validate_arg_double(i,nrhs,prhs,name);
        rays_minimum = (double)*mxGetPr(prhs[++i]);
      }else if(strcmp("Facetwise",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        facet_wise = *v;
      }else if(strcmp("UseParity",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        mxLogical * v = (mxLogical *)mxGetData(prhs[++i]);
        use_parity = *v;
      } else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  Eigen::VectorXi C;
  igl::embree::reorient_facets_raycast(
    V,F,rays_total,rays_minimum,facet_wise,use_parity,false,I,C);

  // Conservative in case FF = F
  Eigen::MatrixXi FF(F.rows(),F.cols());
  for(int i = 0;i<I.rows();i++)
  {
    if(I(i))
    {
      FF.row(i) = (F.row(i).reverse()).eval();
    }else
    {
      FF.row(i) = F.row(i);
    }
  }

  mexErrMsgTxt(nlhs<=2,"Too many output parameters");
  switch(nlhs)
  {
    case 2:
      prepare_lhs_logical(I,plhs+1);
      // Fall through
    case 1:
      prepare_lhs_index(FF,plhs+0);
    default:break;
  }
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#include <igl/principal_curvature.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>


// Los que tienen el tipo de objeto definido son input, el resto output??

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd V,PD1,PD2,PV1,PV2;
  MatrixXi F;
  int radius = 2;
  bool useKring = true;
    
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs,V); // Aquí se pasa del prhs a la matriz V de Eigen
  parse_rhs_index(prhs+1,F); // Aquí se pasa del prhs a la matriz F de Eigen
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Radius",name) == 0)
      {
        validate_arg_scalar(i,nrhs,prhs,name);
        validate_arg_double(i,nrhs,prhs,name);
        radius = (int)*mxGetPr(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
      }
      i++;
    }
  }

  
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,radius,useKring);

  switch(nlhs)
  {
      case 4:
          prepare_lhs_double(PV2,plhs+3);
      case 3:
          prepare_lhs_double(PV1,plhs+2);
      case 2:
          prepare_lhs_double(PD2,plhs+1);
      case 1:
          prepare_lhs_double(PD1,plhs+0);
      default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

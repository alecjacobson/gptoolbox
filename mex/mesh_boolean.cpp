#ifdef MEX

#include <igl/copyleft/boolean/mesh_boolean.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>
#define IGL_NO_CORK
#ifndef IGL_NO_CORK
#  include <igl/copyleft/boolean/mesh_boolean_cork.h>
#endif
#include <igl/copyleft/boolean/string_to_mesh_boolean_type.h>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/C_STR.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>

#include <cstring>

enum BooleanLibType
{
  BOOLEAN_LIB_TYPE_LIBIGL = 0,
  BOOLEAN_LIB_TYPE_CORK = 1,
  BOOLEAN_LIB_TYPE_LIBIGL_TRY_CORK_RESOLVE = 2,
  NUM_BOOLEAN_LIB_TYPES = 3
};

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & VA,
  Eigen::MatrixXi & FA,
  Eigen::MatrixXd & VB,
  Eigen::MatrixXi & FB,
  igl::copyleft::boolean::MeshBooleanType & type,
  BooleanLibType & boolean_lib,
  bool & debug
  )
{
  using namespace std;
  using namespace igl;
  using namespace igl::copyleft::boolean;
  using namespace igl::matlab;
  using namespace igl::copyleft::cgal;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 5, "The number of input arguments must be >=5.");

  const auto & parse_mesh = [](
    const mxArray *prhs[], 
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
  {
    const int dim = mxGetN(prhs[0]);
    mexErrMsgTxt(dim==0 || dim == 3,
      "Mesh vertex list must be #V by 3 list of vertex positions");
    mexErrMsgTxt(mxGetN(prhs[1])==0 || dim == mxGetN(prhs[1]),
      "Mesh \"face\" simplex size must equal dimension");
    parse_rhs_double(prhs,V);
    parse_rhs_index(prhs+1,F);
  };
  parse_mesh(prhs,VA,FA);
  parse_mesh(prhs+2,VB,FB);
  mexErrMsgTxt(mxIsChar(prhs[4]), C_STR("Type should be char"));
  const char * type_str = mxArrayToString(prhs[4]);
  mexErrMsgTxt(string_to_mesh_boolean_type(type_str,type),
    C_STR(type_str << " is not a known boolean operation"));

  {
    int i = 5;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("BooleanLib",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * type_name = mxArrayToString(prhs[++i]);
        if(strcmp("libigl",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL;
#ifndef IGL_NO_CORK
        }else if(strcmp("cork",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_CORK;
        }else if(strcmp("libigl-try-cork-resolve",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL_TRY_CORK_RESOLVE;
#endif
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown BooleanLib: "<<type_name));
        }
      }else if(strcmp("Debug",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        debug = (bool)*mxGetLogicals(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
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
  using namespace igl::copyleft::boolean;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd VA,VB,VC;
  MatrixXi FA,FB,FC;
  VectorXi J,I;
  MeshBooleanType type;
  BooleanLibType boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL;
  bool debug = false;
  parse_rhs(nrhs,prhs,VA,FA,VB,FB,type,boolean_lib,debug);
  if(debug)
  {
    cout<<"parsed input."<<endl;
  }
  switch(boolean_lib)
  {
    case BOOLEAN_LIB_TYPE_LIBIGL:
      mesh_boolean(VA,FA,VB,FB,type,VC,FC,J,I);
      break;
#ifndef IGL_NO_CORK
    case BOOLEAN_LIB_TYPE_CORK:
      mesh_boolean_cork(VA,FA,VB,FB,type,VC,FC);
      break;
    case BOOLEAN_LIB_TYPE_LIBIGL_TRY_CORK_RESOLVE:
    {
      // gcc needs an explicit type here, const auto & fails to reveal template
      // in mesh_boolean call (clang relizes it though)
      const std::function<void(
          const Eigen::Matrix<double, -1, 3>&, 
          const Eigen::Matrix<int, -1, 3>&, 
          Eigen::Matrix<double, -1, 3>&, 
          Eigen::Matrix<int, -1, 3>&, 
          Eigen::Matrix<int, -1, 1>&,
          Eigen::Matrix<int, -1, 1>&)> &
        try_cork_resolve = [](
          const Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> & V,
          const Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> & F,
                Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> & CV,
                Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> & CF,
                Eigen::Matrix<Eigen::VectorXi::Scalar,Eigen::Dynamic,1> & J,
                Eigen::Matrix<Eigen::VectorXi::Scalar,Eigen::Dynamic,1> & I)
      {
        typedef Matrix<MatrixXd::Scalar,Dynamic,3> MatrixX3S;
        typedef Matrix<MatrixXi::Scalar,Dynamic,3> MatrixX3I;
        typedef Matrix<MatrixXi::Scalar,Dynamic,2> MatrixX2I;
        typedef Matrix<MatrixXi::Scalar,Dynamic,1> VectorXI;
        const auto & cork_resolve = [](
          const MatrixX3S &V, 
          const MatrixX3I &F, 
          MatrixX3S &CV, 
          MatrixX3I &CF,
          VectorXI &/*J*/,
          VectorXI &/*I*/)
        {
          Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> _1;
          Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> _2;
          mesh_boolean_cork(V,F,_1,_2,MESH_BOOLEAN_TYPE_RESOLVE,CV,CF);
        };
        // OK if &(CV,CF) = &(V,F)
        const auto & libigl_resolve = [](
          const MatrixX3S &V, 
          const MatrixX3I &F, 
          MatrixX3S &CV, 
          MatrixX3I &CF,
          VectorXI &J,
          VectorXI &I)
        {
          using namespace Eigen;
          MatrixX3S SV;
          MatrixX3I SF;
          MatrixX2I SIF;
          igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
          remesh_self_intersections(V,F,params,SV,SF,SIF,J,I);
          for_each(SF.data(),SF.data()+SF.size(),[&I](int & a){a=I(a);});
          {
            VectorXI UIM;
            remove_unreferenced(SV,SF,CV,CF,UIM);
            for_each(I.data(),I.data()+I.size(),[&UIM](int & a){a=UIM(a);});
          }
        };
        const auto & validate = [](
          const MatrixX3S &V, const MatrixX3I &F)->bool
        {
          MatrixX3S SV;
          MatrixX3I SF;
          MatrixX2I SIF;
          VectorXI SJ,SIM,UIM;
          igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
          params.detect_only = true;
          params.first_only = true;
          remesh_self_intersections(V,F,params,SV,SF,SIF,SJ,SIM);
          return SIM.size() == 0;
        };
        cork_resolve(V,F,CV,CF,J,I);
        if(!validate(CV,CF))
        {
          libigl_resolve(CV,CF,CV,CF,J,I);
        }
        // Cork does not keep track of birth parents (to my knowledge)
        J.setConstant(CF.rows(),1,-1);
        I.setConstant(V.rows(),1,-1);
      };
      mesh_boolean(VA,FA,VB,FB,type,try_cork_resolve,VC,FC,J,I);
      break;
    }
#endif
    default:
      assert(false && "Unknown boolean lib");
      break;
  }
  if(debug)
  {
    cout<<"Computed boolean."<<endl;
  }
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 4:
    {
      prepare_lhs_index(I,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_index(J,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(FC,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(VC,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#endif


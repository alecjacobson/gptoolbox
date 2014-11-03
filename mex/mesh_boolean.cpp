#ifdef MEX
// mex -v -largeArrayDims -DMEX ...
//   CXX='/usr/bin/clang++' LD='/usr/bin/clang++' ...
//   LDOPTIMFLAGS='-O ' ...
//   -I/usr/local/igl/libigl/include ...
//   -I/opt/local/include/eigen3 ...
//   -I/opt/local/include ...
//   -L/opt/local/lib ...
//   -I/usr/local/igl/libigl/external/cork/include ...
//   -L/usr/local/igl/libigl/external/cork/lib -lcork ...
//   -lCGAL -lCGAL_Core -lgmp -lmpfr ...
//   -lboost_thread-mt -lboost_system-mt ...
//   -o mesh_boolean mesh_boolean.cpp

#include <igl/boolean/mesh_boolean.h>
#include <igl/boolean/mesh_boolean_cork.h>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
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
  MeshBooleanType & type,
  BooleanLibType & boolean_lib)
{
  using namespace std;
  using namespace igl;
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
  if(strcmp(type_str,"union")==0)
  {
    type = MESH_BOOLEAN_TYPE_UNION;
  }else if(strcmp(type_str,"intersect")==0)
  {
    type = MESH_BOOLEAN_TYPE_INTERSECT;
  }else if(strcmp(type_str,"minus")==0)
  {
    type = MESH_BOOLEAN_TYPE_MINUS;
  }else if(strcmp(type_str,"xor")==0)
  {
    type = MESH_BOOLEAN_TYPE_XOR;
  }else if(strcmp(type_str,"resolve")==0)
  {
    type = MESH_BOOLEAN_TYPE_RESOLVE;
  }else
  {
    mexErrMsgTxt(false,"Unknown type");
  }

  {
    int i = 5;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      const auto requires_arg = 
        [](const int i, const int nrhs, const char * name)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
      };
      const auto validate_char = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsChar(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires char argument"));
      };
      if(strcmp("BooleanLib",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_char(i,prhs,name);
        const char * type_name = mxArrayToString(prhs[i]);
        if(strcmp("libigl",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL;
        }else if(strcmp("cork",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_CORK;
        }else if(strcmp("libigl-try-cork-resolve",type_name)==0)
        {
          boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL_TRY_CORK_RESOLVE;
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown BooleanLib: "<<type_name));
        }
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

  igl::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd VA,VB,VC;
  MatrixXi FA,FB,FC;
  MeshBooleanType type;
  BooleanLibType boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL;
  parse_rhs(nrhs,prhs,VA,FA,VB,FB,type,boolean_lib);
  switch(boolean_lib)
  {
    case BOOLEAN_LIB_TYPE_LIBIGL:
      mesh_boolean(VA,FA,VB,FB,type,VC,FC);
      break;
    case BOOLEAN_LIB_TYPE_CORK:
      mesh_boolean_cork(VA,FA,VB,FB,type,VC,FC);
      break;
    case BOOLEAN_LIB_TYPE_LIBIGL_TRY_CORK_RESOLVE:
    {
      const auto & try_cork_resolve = [](
        const Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> & V,
        const Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> & F,
              Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> & CV,
              Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> & CF)
      {
        typedef Matrix<MatrixXd::Scalar,Dynamic,3> MatrixX3S;
        typedef Matrix<MatrixXi::Scalar,Dynamic,3> MatrixX3I;
        typedef Matrix<MatrixXi::Scalar,Dynamic,2> MatrixX2I;
        typedef Matrix<MatrixXi::Scalar,Dynamic,1> VectorXI;
        const auto & cork_resolve = [](
          const MatrixX3S &V, const MatrixX3I &F, MatrixX3S &CV, MatrixX3I &CF)
        {
          Eigen::Matrix<Eigen::MatrixXd::Scalar,Eigen::Dynamic,3> _1;
          Eigen::Matrix<Eigen::MatrixXi::Scalar,Eigen::Dynamic,3> _2;
          mesh_boolean_cork(V,F,_1,_2,MESH_BOOLEAN_TYPE_RESOLVE,CV,CF);
        };
        // OK if &(CV,CF) = &(V,F)
        const auto & libigl_resolve = [](
          const MatrixX3S &V, const MatrixX3I &F, MatrixX3S &CV, MatrixX3I &CF)
        {
          using namespace Eigen;
          MatrixX3S SV;
          MatrixX3I SF;
          MatrixX2I SIF;
          VectorXI SJ,SIM,UIM;
          RemeshSelfIntersectionsParam params;
          remesh_self_intersections(V,F,params,SV,SF,SIF,SJ,SIM);
          for_each(SF.data(),SF.data()+SF.size(),[&SIM](int & a){a=SIM(a);});
          {
            remove_unreferenced(SV,SF,CV,CF,UIM);
          }
        };
        const auto & validate = [](
          const MatrixX3S &V, const MatrixX3I &F)->bool
        {
          MatrixX3S SV;
          MatrixX3I SF;
          MatrixX2I SIF;
          VectorXI SJ,SIM,UIM;
          RemeshSelfIntersectionsParam params;
          params.detect_only = true;
          params.first_only = true;
          remesh_self_intersections(V,F,params,SV,SF,SIF,SJ,SIM);
          return SIM.size() == 0;
        };
        cork_resolve(V,F,CV,CF);
        if(!validate(CV,CF))
        {
          libigl_resolve(CV,CF,CV,CF);
        }
      };
      mesh_boolean(VA,FA,VB,FB,type,try_cork_resolve,VC,FC);
      break;
    }
    default:
      assert(false && "Unknown boolean lib");
      break;
  }

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
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


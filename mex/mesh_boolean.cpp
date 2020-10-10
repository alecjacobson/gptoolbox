#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/mesh_boolean_type_to_funcs.h>
#include <igl/copyleft/cgal/BinaryWindingNumberOperations.h>
#include <igl/remove_unreferenced.h>
#ifdef WITH_CORK
#  include <igl/copyleft/cork/mesh_boolean.h>
#endif
#include <igl/copyleft/cgal/string_to_mesh_boolean_type.h>

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
  NUM_BOOLEAN_LIB_TYPES = 2
};

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  std::vector<Eigen::MatrixXd> & Vlist,
  std::vector<Eigen::MatrixXi> & Flist,
  igl::MeshBooleanType & type,
  std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) > & wind_func,
  std::function<int(const int, const int)> & keep_func,
  BooleanLibType & boolean_lib,
  bool & debug
  )
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::copyleft::cgal;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const auto & parse_mesh = [](
    const mxArray *pV, 
    const mxArray *pF, 
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
  {
    const int dim = mxGetN(pV);
    mexErrMsgTxt(dim==0 || dim == 3,
      "Mesh vertex list must be #V by 3 list of vertex positions");
    mexErrMsgTxt(mxGetN(pF)==0 || dim == mxGetN(pF),
      "Mesh \"face\" simplex size must equal dimension");
    parse_rhs_double(&pV,V);
    parse_rhs_index(&pF,F);
  };

  int i = 0;
  while(!mxIsChar(prhs[i]))
  {
    if(mxIsCell(prhs[i]))
    {
      mexErrMsgTxt(
        mxIsCell(prhs[i+1]),
        "Cell input for vertices requires cell input for faces, too"); 
      const int k = mxGetNumberOfElements(prhs[i]);
      mexErrMsgTxt(
        k == mxGetNumberOfElements(prhs[i+1]),
        "Face cell input must be same length");
      for(int j = 0;j<k;j++)
      {
        Vlist.emplace_back();
        Flist.emplace_back();
        parse_mesh(
          mxGetCell(prhs[i],j),
          mxGetCell(prhs[i+1],j),
          Vlist[Vlist.size()-1],Flist[Flist.size()-1]);
      }
    }else if(mxIsNumeric(prhs[i]))
    {
      Vlist.emplace_back();
      Flist.emplace_back();
      parse_mesh(prhs[i],prhs[i+1],Vlist[Vlist.size()-1],Flist[Flist.size()-1]);
    }else
    {
      mexErrMsgTxt(false,"Mesh inputs must be numeric");
    }
    i+=2;
  }

  mexErrMsgTxt(mxIsChar(prhs[i]), C_STR("Type should be char"));
  std::string type_str = std::string(mxArrayToString(prhs[i]));
  if(type_str != "")
  {
    mexErrMsgTxt(string_to_mesh_boolean_type(type_str,type),
      C_STR(type_str << " is not a known boolean operation"));
    igl::copyleft::cgal::mesh_boolean_type_to_funcs(type,wind_func,keep_func);
  }
  i++;

  {
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
#ifdef WITH_CORK
        }else if(strcmp("cork",type_name)==0)
        {
          mexErrMsgTxt( type_str != "",
            "Cork does not support arbitrary extraction functions");
          boolean_lib = BOOLEAN_LIB_TYPE_CORK;
#endif
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown BooleanLib: "<<type_name));
        }
      }else if(strcmp("WindingNumberFilter",name) == 0)
      {
        validate_arg_function_handle(i,nrhs,prhs,name);
        mxArray * handle = const_cast<mxArray*>(prhs[++i]);
        wind_func = [handle](const Eigen::RowVectorXi & w)->int
        {
          int fnlhs = 1;
          mxArray *fplhs[1];
          int fnrhs = 2;
          mxArray *fprhs[2];
          // same as preparing lhs
          fprhs[0] = handle;
          igl::matlab::prepare_lhs_double(w,fprhs+1);
          mexCallMATLAB(fnlhs, fplhs, fnrhs, fprhs,"feval");
          // This might be slow
          int res = 0;
          if( mxIsNumeric(fplhs[0]) )
          {
            res = (int)(mxGetScalar(fplhs[0]));
          }else if( mxIsLogical(fplhs[0]))
          {
            res = (bool)*mxGetLogicals(fplhs[0]);
          }else{
            mexErrMsgTxt(
              false, "WindingNumberFilter should return number or logical");
          }
          // clean up
          //mxDestroyArray(fprhs[0]);
          mxDestroyArray(fprhs[1]);
          mxDestroyArray(fplhs[0]);
          return res;
        };
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
  mexErrMsgTxt(type_str != "" || wind_func,
    "If no type provided, then must provided extraction function");
}

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  std::vector<MatrixXd> Vlist;
  std::vector<MatrixXi> Flist;
  MatrixXd VC;
  MatrixXi FC;
  VectorXi J;
  MeshBooleanType type;
  std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) > wind_func;
  std::function<int(const int, const int)> keep_func = 
    igl::copyleft::cgal::KeepInside();

  BooleanLibType boolean_lib = BOOLEAN_LIB_TYPE_LIBIGL;
  bool debug = false;
  parse_rhs(
    nrhs,
    prhs,
    Vlist,
    Flist,
    type,
    wind_func,
    keep_func,
    boolean_lib,debug);
  if(debug)
  {
    cout<<"parsed input."<<endl;
  }
  switch(boolean_lib)
  {
    case BOOLEAN_LIB_TYPE_LIBIGL:
      igl::copyleft::cgal::mesh_boolean(
        Vlist,Flist,wind_func,keep_func,VC,FC,J);
      break;
#ifdef WITH_CORK
    case BOOLEAN_LIB_TYPE_CORK:
    {
      mexErrMsgTxt(Vlist.size() == 2 && Flist.size() == 2,
        "Must provide exactly two meshes for cork.");
      assert(Vlist.size() == 2);
      assert(Flist.size() == 2);
      igl::copyleft::cork::mesh_boolean(
        Vlist[0],Flist[0],
        Vlist[1],Flist[1],
        type,VC,FC);
      break;
    }
#endif
    default:
      mexErrMsgTxt(false,"Unknown boolean lib.");
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

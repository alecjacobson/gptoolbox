#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>

#include <igl/per_vertex_normals.h>
#include <igl/parallel_for.h>
#include <igl/signed_distance.h>
#include <igl/per_edge_normals.h>
#include <igl/matlab/validate_arg.h>
#include <igl/per_face_normals.h>
#include <igl/WindingNumberAABB.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>

#include <Eigen/Core>
#include <iostream>

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  igl::SignedDistanceType & type)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 3, "The number of input arguments must be >=3.");

  const int dim = mxGetN(prhs[0]);


  parse_rhs_double(prhs,P);
  parse_rhs_double(prhs+1,V);
  parse_rhs_index(prhs+2,F);

  mexErrMsgTxt(P.cols()==3 || P.cols()==2,"P must be #P by (3|2)");
  mexErrMsgTxt(V.cols()==3 || V.cols()==2,"V must be #V by (3|2)");
  mexErrMsgTxt(V.cols()==P.cols(),"dim(V) must be dim(P)");
  mexErrMsgTxt(F.cols()==V.cols(),"F must be #F by dim(V)");

  type = SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("SignedDistanceType",name) == 0)
      {
        validate_arg_char(i,nrhs,prhs,name);
        const char * type_name = mxArrayToString(prhs[++i]);
        if(strcmp("pseudonormal",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
        }else if(strcmp("winding_number",type_name)==0)
        {
          type = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
        }else
        {
          mexErrMsgTxt(false,C_STR("Unknown SignedDistanceType: "<<type_name));
        }
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }
}

//void precompute(
//  const Eigen::MatrixXd & V,
//  const Eigen::MatrixXi & F,
//  const igl::SignedDistanceType sign_type,
//  State & state)
//{
//  using namespace igl;
//  using namespace std;
//  using namespace Eigen;

  // This will remember the data structures for subsequent calls (without an
  // calls in between with different (V,F) or type
  static Eigen::MatrixXd g_V;
  static Eigen::MatrixXi g_F;
  static igl::SignedDistanceType g_sign_type = 
    igl::NUM_SIGNED_DISTANCE_TYPE;
  static igl::AABB<Eigen::MatrixXd,3> g_tree;
  static igl::WindingNumberAABB<
    Eigen::RowVector3d,
    Eigen::MatrixXd,
    Eigen::MatrixXi> g_hier;
  static Eigen::MatrixXd g_FN,g_VN,g_EN;
  static Eigen::MatrixXi g_E;
  static Eigen::VectorXi g_EMAP;

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

  MatrixXd P,V,C,N;
  MatrixXi F;
  VectorXi I;
  VectorXd S;
  SignedDistanceType type;
  parse_rhs(nrhs,prhs,P,V,F,type);

  if(F.rows() > 0)
  {
    switch(V.cols())
    {
      case 2:
      {
        // Persistent data not supported for 2D
        signed_distance(P,V,F,type,S,I,C,N);
        break;
      }
      case 3:
      {
        if(g_sign_type != type || g_V != V || g_F != F)
        {
          g_V = V;
          g_F = F;
          g_sign_type = type;
          // Clear the tree
          g_tree.deinit();

          // Prepare distance computation
          g_tree.init(V,F);
          switch(type)
          {
            default:
              assert(false && "Unknown SignedDistanceType");
            case SIGNED_DISTANCE_TYPE_DEFAULT:
            case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
              g_hier.set_mesh(V,F);
              g_hier.grow();
              break;
            case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
              // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
              // [Bærentzen & Aanæs 2005]
              per_face_normals(V,F,g_FN);
              per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,
                g_FN,g_VN);
              per_edge_normals(
                V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,
                g_FN,g_EN,g_E,g_EMAP);
              break;
          }
        }

        N.resize(P.rows(),3);
        S.resize(P.rows(),1);
        I.resize(P.rows(),1);
        C.resize(P.rows(),3);
        //for(int p = 0;p<P.rows();p++)
        igl::parallel_for(P.rows(),[&](const int p)
        {
          const Eigen::RowVector3d q(P(p,0),P(p,1),P(p,2));
          double s,sqrd;
          Eigen::RowVector3d c;
          int i;
          switch(type)
          {
            default:
              assert(false && "Unknown SignedDistanceType");
            case SIGNED_DISTANCE_TYPE_DEFAULT:
            case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
              signed_distance_winding_number(
                g_tree,g_V,g_F,g_hier,q,s,sqrd,i,c);
              break;
            case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
            {
              RowVector3d n(0,0,0);
              signed_distance_pseudonormal(
                g_tree,g_V,g_F,g_FN,g_VN,g_EN,g_EMAP,
                q,s,sqrd,i,c,n);
              N.row(p) = n;
              break;
            }
          }
          I(p) = i;
          S(p) = s*sqrt(sqrd);
          C.row(p) = c;
        },10000);
        break;
      }
    }
  }

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 4:
    {
      prepare_lhs_double(N,plhs+3);
      // Fall through
    }
    case 3:
    {
      prepare_lhs_double(C,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(I,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(S,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

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
#include <igl/matlab_format.h>

#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/barycenter.h>
#include <igl/slice_mask.h>
#include <igl/slice.h>
#include <igl/remove_unreferenced.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign.h>
#include <igl/copyleft/cgal/assign_scalar.h>
// Inputs:
//   S  3 by #S list of corner values corresponding to first three rows of BC
//     output.
// Outputs:
//   FF  #FF by 3 list of triangle indices into rows of BC
//   J  #FF list of indices into cols of S
//   BC
Eigen::VectorXi counts;
void upper_envelope(
  const Eigen::MatrixXd & S,
  Eigen::MatrixXi & FF,
  Eigen::VectorXi & J,
  Eigen::MatrixXd & BC)
{
  // THIS BASE CASE IS A PERFORMANCE OPTIMIZATION
  const int m = S.cols();
  typedef Eigen::MatrixXd::Scalar Scalar;
  // Determine if base case
  typedef Eigen::MatrixXd::Index Index;
  typedef Eigen::Matrix<Scalar,1,Eigen::Dynamic> RowVectorXS;
  typedef Eigen::Matrix<Scalar,3,1> Vector3S;
  const Vector3S Smax = S.rowwise().maxCoeff();
  Eigen::Array<bool,3,Eigen::Dynamic> is_max(3,m);
  // Consider all functions
  for(int j = 0;j<m;j++)
  {
    bool all_max = true;
    // consider each corner
    for(int i = 0;i<3;i++)
    {
      // *Strict* Equality
      is_max(i,j) = S(i,j) == Smax(i);
      all_max = all_max && is_max(i,j);
    }
    if(all_max)
    {
      // Yay! found a base case (btw, there could be more than one candidate)
      FF.resize(1,3); FF<<0,1,2;
      J.resize(1,1); J<<j;
      BC.resize(3,3); BC<<1,0,0,0,1,0,0,0,1;
      counts(0)++;
      return;
    }
  }

  // THIS BASE CASE IS ALSO A PERFORMANCE OPTIMIZATION
  // for the case where there are only two functions interacting, is it really
  // helpful to know `is_max` ? I.e., is that even enough information to
  // determine:
  //   1) which two functions are interacting? and
  //   2) that, in fact, only two functions are interacting?
  Index xj;
  {
    Index _;
    S.maxCoeff(&_,&xj);
  }
  // Alright so xj reaches the max value. 
  // Better lower bound? One per corner?
  Eigen::VectorXi overlap(m);
  {
    int o = 1;
    overlap(0) = xj;
    Scalar low = S.col(xj).minCoeff();
    for(int j = 0;j<m;j++)
    {
      if(j == xj) { continue; }
      if(S.col(j).maxCoeff()>low)
      {
        overlap(o++) = j;
        // update lower bound? Can we shrink the extent?
      }
    }
    overlap.conservativeResize(o);
  }
  counts(overlap.size()-1)++;
  typedef CGAL::Epeck Kernel;
  typedef Kernel::FT ExactScalar;
  //std::cout<<"overlap.size(): "<<overlap.size()<<std::endl;
  if(overlap.size() == 2)
  {
    const int yj = overlap(1);
    Eigen::Matrix<ExactScalar, Eigen::Dynamic, 1>  Dxy(3,1);
    for(int i = 0;i<3;i++)
    {
      Dxy(i) = ExactScalar(S(i,xj))-ExactScalar(S(i,yj));
    }
    // edge where sign is the same 
    int e0 = -1;
    for(int e = 0;e<3;e++)
    {
      if((Dxy((e+1)%3)>0) == (Dxy((e+2)%3)>0))
      {
        e0 = e;
        break;
      }
    }
    assert(e0>=0);
    //
    //   v2
    //    o
    //    |\
    //    | \
    // e1 a  \ e0
    //    |\  \
    //    | \  \
    //    o--b--o
    //  v0   e2  v1
    //
    const int e1 = (e0+1)%3;
    const int e2 = (e0+2)%3;
    BC.resize(5,3);
    BC.topLeftCorner(3,3).setIdentity();
    ExactScalar s20 = -Dxy(e2)/(Dxy(e0)-Dxy(e2));
    ExactScalar s01 = -Dxy(e0)/(Dxy(e1)-Dxy(e0));
    const int a = 3;
    const int b = 4;
    igl::copyleft::cgal::assign_scalar(1.0-s20,BC(a,e2));
    igl::copyleft::cgal::assign_scalar(    s20,BC(a,e0));
                                               BC(a,e1) = 0;
    igl::copyleft::cgal::assign_scalar(1.0-s01,BC(b,e0));
    igl::copyleft::cgal::assign_scalar(    s01,BC(b,e1));
                                               BC(b,e2) = 0;
    FF.resize(3,3);
    FF<<
      a,e0,b,
      b,e1,e2,
      b,e2,a;
    J.resize(3,1);
    const bool flag = Dxy(e0)>0;
    J(0) = flag?xj:yj;
    J(1) = flag?yj:xj;
    J(2) = flag?yj:xj;


    return;
  }
  
  //
  // Ja = {j∈[1,m] | S(a,j) = max(S(a,:))} 
  //   // could be many triangles kissing at a
  // Jb = {j∈Ja | S(b,j) = max(S(b,Ja))}
  //   // could still be a fan around [a,b] edge
  // Jc = {j∈Jb | S(c,j) = max(S(c,Jb))}
  //   // CLAIM: Jc are co-planar and part of upper-envelope at a
  //

  const auto brute_force = [](
    const Eigen::MatrixXd & S,
    Eigen::MatrixXi & FF,
    Eigen::VectorXi & J,
    Eigen::MatrixXd & BC)
  {
    const int m = S.cols();
    // THIS IS O(m²) EVEN IF THE OUTPUT IS O(1)
    // Consider case where upper envelope of m-1 faces is O(m²) and then add an
    // mth function above all of them.

    // General case
    // Create "Height Field" triangle soup
    Eigen::Matrix<Scalar,Eigen::Dynamic,3> HV(3*m,3);
    Eigen::MatrixXi HF(m,3);
    for(int j = 0;j<m;j++)
    {
      for(int i = 0;i<3;i++)
      {
        HF(j,i) = 3*j+i;
        HV(3*j+i,0) = i==1;
        HV(3*j+i,1) = i==2;
        HV(3*j+i,2) = S(i,j);
      }
    }
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::Plane_3 Plane_3;
    typedef Eigen::Matrix<ExactScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXES;
    MatrixXES SV;
    Eigen::MatrixXi SF;
    Eigen::VectorXi SJ;
    {
      Eigen::MatrixXi _1;
      Eigen::VectorXi _2;
      igl::copyleft::cgal::remesh_self_intersections(
        HV,HF,{false,false,true},SV,SF,_1,SJ,_2);
    }
    const auto d = [](const MatrixXES & X) -> Eigen::MatrixXd
    {
      Eigen::MatrixXd Xd;
      igl::copyleft::cgal::assign(X,Xd);
      return Xd;
    };
    //std::cout<<"  "<<igl::matlab_format(d(SV),"SV")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((SF.array()+1).eval(),"SF")<<std::endl;
    MatrixXES SC;
    igl::barycenter(SV,SF,SC);
    // Unfortunately, the triangle with "heighest" barycenter is not necessarily
    // part of the upper envelope.

    // This is O(Mm) where M number of SF, which could be O(m²), so this is O(m³).
    // Yikes!
    //
    std::vector<Point_3> HP;HP.reserve(HV.rows());
    for(int i = 0;i<HV.rows();i++) { HP.emplace_back(HV(i,0),HV(i,1),HV(i,2)); }

    // Consider each output triangle
    Eigen::Array<bool,Eigen::Dynamic,1> keep(SF.rows(),1);
    for(int f = 0;f<SF.rows();f++)
    {
      // centroid
      Point_3 cf(SC(f,0),SC(f,1),SC(f,2));
      keep(f) = true;
      // Consider every input function
      for(int j = 0;j<m;j++)
      {
        // Skip the function of f
        if(j == SJ(f)) { continue; }
        const auto on = CGAL::orientation(HP[HF(j,0)],HP[HF(j,1)],HP[HF(j,2)],cf);
        if(on == CGAL::NEGATIVE || (on == CGAL::COPLANAR && j<SJ(f)))
        {
          keep(f) = false;
          break;
        }
      }
    }
    Eigen::MatrixXi KF;
    igl::slice_mask(SF,keep,1,KF);
    igl::slice_mask(SJ,keep,1,J);
    // Remove unreferenced
    MatrixXES UV;
    {
      Eigen::VectorXi _;
      igl::remove_unreferenced(SV,KF,UV,FF,_);
    }
    //std::cout<<"  "<<igl::matlab_format(d(UV),"UV")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((FF.array()+1).eval(),"UF")<<std::endl;
    // Augh, make the original corners come first
    Eigen::VectorXi I = Eigen::VectorXi::LinSpaced(UV.rows(),0,UV.rows()-1);
    for(int i = 0;i<UV.rows();i++)
    {
      for(int j = 0;j<3;j++)
      {
        if(
          UV(i,0) == ExactScalar(j==1) &&
          UV(i,1) == ExactScalar(j==2))
        {
          std::swap(I(i),I(j));
        }
      }
    }
    // matlibberish mindfuck
    Eigen::VectorXi II(I.size());
    MatrixXES IV(UV.rows(),UV.cols());
    for(int i = 0;i<I.size();i++)
    {
      IV.row(i) = UV.row(I(i));
      II(I(i)) = i;
    }
    for(int f = 0;f<FF.rows();f++)
    {
      for(int c = 0;c<FF.cols();c++)
      {
        FF(f,c) = II(FF(f,c));
      }
    }
    //std::cout<<"  "<<igl::matlab_format((I.array()+1).eval(),"I")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((II.array()+1).eval(),"II")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format(d(IV),"IV")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((FF.array()+1).eval(),"FF")<<std::endl;

    BC.resize(IV.rows(),3);
    for(int i = 0;i<IV.rows();i++)
    {
      const auto alpha = ExactScalar(1.0)-IV(i,0)-IV(i,1);
      igl::copyleft::cgal::assign_scalar(alpha,BC(i,0));
      igl::copyleft::cgal::assign_scalar(IV(i,0),BC(i,1));
      igl::copyleft::cgal::assign_scalar(IV(i,1),BC(i,2));
    }
    //std::cout<<"  "<<igl::matlab_format(BC,"BC")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((FF.array()+1).eval(),"FF")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((J.array()+1).eval(),"J")<<std::endl;
  };

  Eigen::MatrixXd oS;
  igl::slice(S,overlap,2,oS);
    //std::cout<<"  "<<igl::matlab_format(S,"S")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format((overlap.array()+1).eval(),"overlap")<<std::endl;
    //std::cout<<"  "<<igl::matlab_format(oS,"oS")<<std::endl;
  Eigen::VectorXi oJ;
  brute_force(oS,FF,oJ,BC);
  // reindex
  J.resize(oJ.size());
  for(int i = 0;i<oJ.size();i++) { J(i) = overlap(oJ(i)); }
}

void upper_envelope(
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & S,
  Eigen::MatrixXi & FF,
  Eigen::VectorXi & J,
  Eigen::VectorXi & FJ,
  Eigen::SparseMatrix<double> & BC)
{
  // number of vertices
  int n = S.rows();
  counts.setZero(S.cols());
  std::vector<Eigen::Triplet<double> > BCIJV;
  BCIJV.reserve(2*n);
  // original vertices come first
  for(int vi = 0;vi<n;vi++){ BCIJV.emplace_back(vi,vi,1.0); }

  // amortized dynamic number of rows
  const auto append_bottom = [](
     const auto & Y,
     auto & X,
     int & rX)
  {
    assert(X.cols() == Y.cols() && "X.cols() should match Y.cols()");
    // base case I'm not sure is necessary
    if( Y.rows() == 0 ) { return; }
    // add rows if necessary
    if(rX+Y.rows() > X.rows())
    {
      X.conservativeResize((X.rows()+Y.rows())*2,Eigen::NoChange);
    }
    X.block(rX,0,Y.rows(),X.cols()) = Y;
    rX += Y.rows();
  };
  FF.resize(2*F.rows(),F.cols()); int rFF = 0;
  J .resize(2*F.rows());           int rJ = 0;
  FJ.resize(2*F.rows());          int rFJ = 0;

  const int m = S.cols();
  // loop over each triangle
  const int nf = F.rows();
  for(int f = 0;f<nf;f++)
  {
    //std::cout<<f<<":"<<std::endl;
    Eigen::MatrixXd Sf(3,m);
    Sf<<S.row(F(f,0)), S.row(F(f,1)), S.row(F(f,2));
    Eigen::MatrixXi FFf;
    Eigen::VectorXi Jf;
    Eigen::MatrixXd BCf;
    upper_envelope(Sf,FFf,Jf,BCf);
    // reindex FFf
    for(int i = 0;i<FFf.rows();i++)
    {
      for(int j = 0;j<FFf.cols();j++)
      {
        if(FFf(i,j)<3)
        {
          // first three are original corners of F(f,:)
          FFf(i,j) = F(f,FFf(i,j));
        }else
        {
          // otherwise will be appended to BC as new vertices
          FFf(i,j) += n;
        }
      }
    }
    append_bottom(FFf,FF,rFF);
    append_bottom(Eigen::VectorXi::Constant(FFf.rows(),1,f),FJ,rFJ);
    append_bottom(Jf,J,rJ);
    // append to BC
    for(int i = 0;i<BCf.rows();i++)
    {
      for(int j = 0;j<BCf.cols();j++)
      {
        BCIJV.emplace_back(n,F(f,j),BCf(i,j));
      }
      // increment count of number of vertices
      n++;
    }
  }
  // Build BC
  BC.resize(n,S.rows());
  BC.setFromTriplets(BCIJV.begin(),BCIJV.end());
  // Trim dynamic arrays
  FF.conservativeResize(rFF,Eigen::NoChange);
  FJ.conservativeResize(rFJ,Eigen::NoChange);
  J.conservativeResize(rJ,Eigen::NoChange);
}

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  mexUnlock();
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  const bool stitch_all = false;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd S;
  mexErrMsgTxt(nrhs>=3,"nrhs should be at least 3");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);
  parse_rhs_double(prhs+2,S);
  mexErrMsgTxt(V.rows() == S.rows(),"V and S should have same # rows");
  // number of functions
  const int m = S.cols();
  // simplex size
  const int ss = F.cols();

  Eigen::MatrixXd VV;
  Eigen::MatrixXi FF;
  Eigen::SparseMatrix<double> BC;
  Eigen::VectorXi J,FJ;

  upper_envelope(F,S,FF,J,FJ,BC);
  // Could consider outputting this as we go, if ever the bottleneck
  VV = BC*V;

  mexErrMsgTxt(!stitch_all,"Stitch all not yet supported");

  switch(nlhs)
  {
    case 6:
      prepare_lhs_double(counts.cast<double>().eval(),plhs+5);
    case 5:
      prepare_lhs_double(BC,plhs+4);
    case 4:
      prepare_lhs_index(FJ,plhs+3);
    case 3:
      prepare_lhs_index(J,plhs+2);
    case 2:
      prepare_lhs_index(FF,plhs+1);
    case 1:
      prepare_lhs_double(VV,plhs+0);
    default:
      break;
  }



  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

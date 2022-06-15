#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab_format.h>
#include <igl/matlab/prepare_lhs.h>
#include <random>
#include <igl/blue_noise.h>
#include <Eigen/Core>

// Return a random point in the annulus centered at the origin with inner radius
// r and outer radius R
Eigen::RowVector2d random_point_in_annulus(const double r, const double R)
{
  // rejection sampling 
  Eigen::RowVector2d p = Eigen::RowVector2d::Random()*R;
  const double l2 = p.squaredNorm();
  if(l2 < r*r || l2>R*R)
  {
    return random_point_in_annulus(r,R);
  }else
  {
    return p;
  }
}

// Implements "Fast Poisson Disk Sampling in Arbitrary Dimensions" [Robert
// Bridson 2007] to compute blue noise sampling of a given box according to a
// given minimum distance.
//
// Inputs:
//   w  width and 
//   h  height of box to sample
//   r  Poisson disk radius: minimum distance between samples
//   k  number of rejection sampling attempts before acccepting candidate {30}
void blue_noise(
  const double bw,
  const double bh,
  const double r,
  const int k,
  Eigen::MatrixXd & points)
{
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  // Initial point
  Eigen::RowVector2d p0( dis(gen)*bw, dis(gen)*bh );
  // Step 0
  const double h = r/sqrt(2.); // 2 as in 2D
  // grid(x/h,y/h)
  Eigen::MatrixXi grid = Eigen::MatrixXi::Constant(ceil(bw/h)+1,ceil(bh/h)+1,-1);
  int num_points = 0;
  points.resize(bw*bh/r/r*2/3,2);
  std::vector<int> active;
  // Step 1
  const auto add_point = 
    [&points,&active,&grid,&num_points,&h](const Eigen::RowVector2d & p)
  {
    // Expand array if necessary
    if(num_points >= points.rows())
    {
      points.conservativeResize(2*points.rows()+1,points.cols());
    }
    points.row(num_points) = p;
    active.push_back(num_points);
    grid(p(0)/h,p(1)/h) = num_points;
    num_points++;
  };
  add_point(p0);
  points.conservativeResize(num_points,points.cols());
  while(!active.empty())
  {
    // choose a random index in active
    const int a = dis(gen)*active.size();
    bool found_new_sample;
    // generate up to k points near the (now freshly swapped) last element in active
    Eigen::RowVector2d q = points.row(active[a]);
    for(int i = 0;i<k;i++)
    {
      Eigen::RowVector2d p = random_point_in_annulus(r,2*r);
      p += q;
      //p.x = periodic.x ? (p.x+box.w) % box.w : p.x
      //p.y = periodic.y ? (p.y+box.h) % box.h : p.y
      const int gx = floor(p(0)/h);
      const int gy = floor(p(1)/h);
      // check if in grid
      if(gx<0 || gy<0 || gx>=grid.rows() || gy>=grid.cols()){ continue; }
      // check if points too close
      found_new_sample = true;
      for(int nx = std::max(0,gx-2); nx < std::min(  (int)grid.rows(),gx+3);nx++)
      {
        for(int ny = std::max(0,gy-2); ny < std::min((int)grid.cols(),gy+3);ny++)
        {
          if(grid(nx,ny) == -1){ continue; }
          if((p-points.row(grid(nx,ny))).squaredNorm()<r*r)
          { 
            found_new_sample = false; 
          }
          if(!found_new_sample) break;
        }
        if(!found_new_sample) break;
      }
      if(found_new_sample)
      {
        add_point(p);
        break;
      }
    }
    if(!found_new_sample)
    {
      // remove this point rom active list
      std::swap(active[a],active[active.size()-1]);
      active.pop_back();
    }
  }
  points.conservativeResize(num_points,2);
}

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexErrMsgTxt(nrhs>=2,"nrhs should be ≥2");
  Eigen::VectorXi FI;
  Eigen::MatrixXd B,P;
  const bool box2d = mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1;
  if(box2d)
  {
    mexErrMsgTxt(lhs<=1,"nlhs should be ≤1");
    Eigen::MatrixXd V;
    parse_rhs_double(prhs+0,V);
    const int dim = V.cols();
    int k = 30;
    mexErrMsgTxt(dim==2,"dim should be == 2");
    if(nrhs>3)
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
        if(strcmp("K",name) == 0)
        {
          validate_arg_double(i,nrhs,prhs,name);
          validate_arg_scalar(i,nrhs,prhs,name);
          double * v = (double *)mxGetData(prhs[++i]);
          k = *v;
        }else
        {
          mexErrMsgTxt(C_STR("Unsupported parameter: "<<name));
        }
        i++;
      }
    }

    Eigen::RowVectorXd minV = V.array().colwise().minCoeff();
    Eigen::RowVectorXd maxV = V.array().colwise().maxCoeff();
    const double r = (double)*mxGetPr(prhs[1]);
    blue_noise(maxV(0)-minV(0),maxV(1)-minV(1),r,k,P);
    if(P.rows()>0)
    {
      P.array().rowwise() += minV.array();
    }
  }else
  {
    Eigen::MatrixXd V;
    parse_rhs_double(prhs+0,V);
    const int dim = V.cols();
    mexErrMsgTxt(dim==3 || dim==2,"dim should be == 2|3");
    if(dim == 2)
    {
      V.conservativeResize(V.rows(),3);
      V.col(2).setZero();
    }
    Eigen::MatrixXi F;
    parse_rhs_index(prhs+1,F);
    const double r = (double)*mxGetPr(prhs[2]);
    igl::blue_noise(V,F,r,B,FI,P);
    P.conservativeResize(P.rows(),dim);
  }
  

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_double(B,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(FI,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_double(P,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}


#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 

#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point ;
//
// Setup an enriched polyhedron type which stores an id() field in the items
//
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface_mesh; 
typedef Surface_mesh::Halfedge_handle Halfedge_handle ;
typedef Surface_mesh::Vertex_handle   Vertex_handle ;
namespace SMS = CGAL::Surface_mesh_simplification ;
typedef SMS::Edge_profile<Surface_mesh> Profile ;
// The following is a Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Stats
{
  Stats() 
    : collected(0)
    , processed(0)
    , collapsed(0)
    , non_collapsable(0)
    , cost_uncomputable(0) 
    , placement_uncomputable(0) 
  {} 
  
  std::size_t collected ;
  std::size_t processed ;
  std::size_t collapsed ;
  std::size_t non_collapsable ;
  std::size_t cost_uncomputable  ;
  std::size_t placement_uncomputable ; 
} ;
struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor( Stats* s) : stats(s){} 
  // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const&, boost::optional<double> const& )
  {
    ++ stats->collected ;
    //std::cerr << "\rEdges collected: " << stats->collected << std::flush ;
  }                
  
  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(Profile const&          
                 ,boost::optional<double> cost
                 ,std::size_t             initial
                 ,std::size_t             current
                 )
  {
    ++ stats->processed ;
    if ( !cost )
      ++ stats->cost_uncomputable ;
      
    //if ( current == initial )
    //  std::cerr << "\n" << std::flush ;
    //std::cerr << "\r" << current << std::flush ;
  }                
  
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(Profile const&          
                   ,boost::optional<Point>  placement
                   )
  {
    if ( !placement )
      ++ stats->placement_uncomputable ;
  }                
  
  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable( Profile const& )
  {
    ++ stats->non_collapsable;
  }                
  
  // Called AFTER each edge has been collapsed
  void OnCollapsed( Profile const&, Vertex_handle )
  {
    ++ stats->collapsed;
  }                
  
  Stats* stats ;
} ;

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::copyleft::cgal;
  using namespace Eigen;
  MatrixXd V,W;
  MatrixXi F,G;
  bool adaptive = false;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,F);
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  mexErrMsgTxt(
    mxIsDouble(prhs[2]) && mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1,
    "fraction to decimate should be scalar");
  double ratio = * mxGetPr(prhs[2]);
  mexErrMsgTxt((ratio>0 && ratio<1),"Ratio should be in (0,1)");
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Adaptive",name) == 0)
      {
        validate_arg_logical(i,nrhs,prhs,name);
        validate_arg_scalar(i,nrhs,prhs,name);
        adaptive = (bool)*mxGetLogicals(prhs[++i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  Surface_mesh surface_mesh; 
  if(!mesh_to_polyhedron(V,F,surface_mesh))
  {
    mexErrMsgTxt("Failed to convert (V,F) to cgal polyhedron");
  }

  // The items in this polyhedron have an "id()" field 
  // which the default index maps used in the algorithm
  // need to get the index of a vertex/edge.
  // However, the Polyhedron_3 class doesn't assign any value to
  // this id(), so we must do it here:
  int index = 0 ;
  
  for( Surface_mesh::Halfedge_iterator eb = surface_mesh.halfedges_begin()
     , ee = surface_mesh.halfedges_end()
     ; eb != ee
     ; ++ eb
     ) 
    eb->id() = index++;
  index = 0 ;
  for( Surface_mesh::Vertex_iterator vb = surface_mesh.vertices_begin()
     , ve = surface_mesh.vertices_end()
     ; vb != ve
     ; ++ vb
     ) 
    vb->id() = index++;
    
  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
 
  Stats stats ;
  My_visitor vis(&stats) ;
  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  if(adaptive)
  {
    int r = SMS::edge_collapse
             (surface_mesh
             ,stop
             ,CGAL::parameters::visitor      (vis)
             );
  }else
  {
    int r = SMS::edge_collapse
             (surface_mesh
             ,stop
             ,CGAL::parameters::visitor      (vis)
             .get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
             .get_placement(SMS::Midpoint_placement<Surface_mesh>())
             );
  }

  polyhedron_to_mesh(surface_mesh,W,G);

  switch(nlhs)
  {
    case 2:
      prepare_lhs_index(G,plhs+1);
    case 1:
      prepare_lhs_double(W,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}

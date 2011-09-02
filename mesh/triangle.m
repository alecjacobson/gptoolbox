function [TV,TF,TN,VV,VE,VRP,VRD] = triangle(varargin)
  % TRIANGLE interface to the Triangle Program.
  %
  % [TV,TF,TN] = triangle(filename,options)  Triangulates file in filename
  % [TV,TF,TN] = triangle(V,options)  Triangulates convex hull of given points
  % [TV,TF,TN] = triangle(V,F,options)  Refine an existing mesh
  % [TV,TF,TN] = triangle(V,E,H,options)  Triangulates a planar straight line graph
  %
  % Inputs:
  %   filename  triangulates existing .poly or .node file
  %   V  #V by 3 constrained points of triangulation
  %   F  #F by 3 list of triangle indices
  %   E  #E by 2 list of constraint edge indices
  %   H  #H by 2 hole points
  %   options
  %     'Quality'  Quality mesh generation with no angles smaller than 20
  %       degrees. An alternate minimum angle may be specified after the `q'.
  %     'MaxArea'  Imposes a maximum triangle area constraint. A fixed area 
  %       constraint (that applies to every triangle) may be specified, or a 
  %       list of #F area constraints (negative area means triangle is left 
  %       unconstrained)
  %     'ConstrainRefine' only if refining an existing mesh, means that edges
  %       should not be eliminated. Can specify edges after this option,
  %       otherwise all edges of existing mesh are used.
  %     'NoBoundarySteiners' Prohibits the insertion of Steiner points on
  %       the mesh boundary. 
  %     'NoEdgeSteiners' prohibits insertion of Steiner points on any segment, 
  %       including internal segments.
  %     'MaxSteiners' specifies max number of allowed Steiner points.
  % Outputs:
  %   TV  #TV by 2 list of mesh vertices
  %   TF  #TF by 3 list of triangle indices
  %   TN  #TN by 3 list of triangle neighbors
  %   VV  #VV by 2 list of voronoi vertices
  %   VE  #VE by 2 list of voronoi internal edges
  %   VRP  #VRP by 1 list of voronoi infinite ray start points
  %   VRD  #VRP by 2 list of voronoi infinite direction unit vectors
  %
  % Note: In the command line program, Triangle, you can specify that you want
  % to refine an existing mesh AND constrain the triangulation. I'm not sure
  % what this means exactly as far as how the inputs used. It is not yet
  % supported here.
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: delaunay, DelaunayTri
  %

  % default values
  quality = -1;
  max_area = -1;
  no_boundary_steiners = false;
  no_edge_steiners = false;
  max_steiners = -1;
  % possible actions are:
  %   'TriangulatePoly'
  %   'RefineMesh'
  triangulate_poly = false;
  triangulate_existing = false;
  refine_mesh = false;
  neighbors = nargout >= 3;
  voronoi = nargout >= 4;

  % parse inputs
  if nargin >=1 && ischar(varargin{1})
    filename = varargin{1};
    % strip file extension to get prefix
    if ~isempty(regexp(filename,'\.poly$'))
      prefix = regexprep(filename,'\.poly$','');
      triangulate_poly = true;
    elseif ~isempty(regexp(filename,'\.node$'))
      prefix = regexprep(filename,'\.node$','');
    else
      prefix = file_name;
    end
    triangulate_existing = true;
    ii = 2;
  elseif( nargin == 1 || (nargin > 1 && ischar(varargin{2})))
    V = varargin{1};
    ii = 2;
  elseif( nargin == 2 || (nargin > 2 && ischar(varargin{3})))
    refine_mesh = true;
    V = varargin{1};
    MF = varargin{2};
    ii = 3;
  elseif( nargin == 3 || (nargin > 3 && ischar(varargin{4})))
    triangulate_poly = true;
    V = varargin{1};
    PE = varargin{2};
    PH = varargin{3};
    ii = 4;
  end

  % parse options
  while(ii <= nargin)
    if( strcmp(varargin{ii},'Quality') == 1)
      if( (ii+1)<=nargin && ~ischar(varargin{ii+1}))
        ii = ii + 1;
        assert(ii <= nargin);
        quality = varargin{ii};
      else
        quality = 20;
      end
    elseif( strcmp(varargin{ii},'MaxArea') == 1)
      ii = ii + 1;
      assert(ii <= nargin);
      max_area = varargin{ii};
    elseif( strcmp(varargin{ii},'NoBoundarySteiners') == 1)
      no_boundary_steiners = true;
    elseif( strcmp(varargin{ii},'NoEdgeSteiners') == 1)
      no_edge_steiners = true;
    elseif( strcmp(varargin{ii},'MaxSteiners') == 1)
      ii = ii + 1;
      assert(ii <= nargin);
      max_steiners = varargin{ii};
    elseif( strcmp(varargin{ii},'ConstrainRefine') == 1)
      if(refine_mesh)
        triangulate_poly = true;
        % Don't think holes are supported here, if they are then
        % 'ConstrainRefine' should be able to take two arguements
        PH = [];
        if( (ii+1)<=nargin && ~ischar(varargin{ii+1}))
          ii = ii + 1;
          assert(ii <= nargin);
          PE = varargin{ii};
        else
          PE = edges(MF);
        end
      else
        warning('Ignoring ''ConstrainRefine'' option');
      end
    else
      error([varargin{ii} ' is not a valid option']);
    end
    ii = ii + 1;
  end


  % warnings
  if(quality > 34)
    warning( ...
      ['Quality: ' num2str(quality) ' > 34, triangle might not terminate']);
  end

  % command line version
  % get a free temporary prefix
  if ~exist('prefix','var')
    prefix = tempprefix();
  end
  
  % build command line parameters string
  params = '-';

  if ~triangulate_existing
    writeNODE([prefix '.node'],V);
  end

  if triangulate_poly
    params = [params 'p'];
    if ~triangulate_existing
      % print poly file
      writePOLY([prefix '.poly'],[],PE,PH);
    end
  end

  if(refine_mesh)
    params = [params 'r'];
    % print .ele file
    writeELE([prefix '.ele'],MF);
  end

  if(quality >= 0)
    params = [params 'q' num2str(quality)];
  end
  if(max(size(max_area)) > 1)
    params = [params 'a'];
  end
  if(max(size(max_area)) > 1)
    error('Area constraints specified per triangle not yet supported...');
  elseif(max_area >= 0)
    params = [params 'a' num2str(max_area)];
  end
  if(max_steiners >= 0)
    params = [params 'S' num2str(max_steiners)];
  end
  if(no_boundary_steiners && ~no_edge_steiners)
    params = [params 'Y'];
  end
  if(no_edge_steiners)
    params = [params 'YY'];
  end
  if(voronoi)
    params = [params 'v'];
  end
  if(neighbors)
    params = [params 'n'];
  end

  path_to_triangle = '/opt/local/bin/triangle';
  command = [path_to_triangle ' ' params ' ' prefix];
  command
  [status, result] = system( command );
  if(status ~= 0)
     error(result);
  end

  % read outputs from files
  TV = readNODE([prefix '.1.node']);
  TF = readELE([prefix '.1.ele']);
  % Triangle likes to use 1-indexed though .ele reader is 0-indexed
  if(( min(TF(:)) > 1) && (max(TF(:)) > size(TV,1)))
    TF = TF-1;
  end

  if(nargout > 2)
    TN = [];
    VV = [];
    VE = [];
    VRP = [];
    VRD = [];
  end

  if(neighbors)
    TN = readELE([prefix '.1.neigh']);
    % Triangle likes to use 1-indexed though .ele reader is 0-indexed
    if(( min(TN(:)) > 1) && (max(TN(:)) > size(TF,1)))
      TN = TN-1;
    end
  end

  if(voronoi)
    VV = readNODE([prefix '.1.v.node']);
    [VE,VRP,VRD] = readEDGE([prefix '.1.v.edge']);
  end

  % delete temporary output files
  delete([prefix '.1.node']);
  delete([prefix '.1.ele']);
  if(exist([prefix '.1.poly'],'file'))
    delete([prefix '.1.poly']);
  end
  if(exist([prefix '.1.v.node'],'file'))
    delete([prefix '.1.v.node']);
  end
  if(exist([prefix '.1.v.edge'],'file'))
    delete([prefix '.1.v.edge']);
  end

  % Delete temporary input files
  if ~triangulate_existing
    delete([prefix '.node']);
  end 

  if triangulate_poly && ~triangulate_existing
    delete([prefix '.poly']);
  end

  if refine_mesh && ~triangulate_existing
    delete([prefix '.ele']);
  else

end


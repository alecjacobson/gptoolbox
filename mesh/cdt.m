function [VV,TT,FF,TN,IFF] = cdt(varargin)
  % CDT compute the constrained delaunay triangulation of a given mesh using
  % tetgen.
  %
  % [VV,TT,FF,TN,IFF] = cdt(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of mesh positions
  %   F  #F by dim list of face indices
  %   Optional:
  %     'TriangleFlags' followed by a string of additional triangle flags, e.g.
  %       '-q'
  %     'TetgenFlags' followed by a string of additional tetgen flags, e.g.
  %       '-J'
  %     'UseBoundingBox' followed by whether to use an explicit bounding box
  %       rather than tetgen's convex hull
  %     'BoundingBoxPush' followed by a push factor to enlargen the bounding
  %       box around its centroid {1.01}
  %     'BoundingBoxUpsample' followed by a factor f, to upsample the bounding
  %       box faces until min(doublearea(BV,BF))<f*max(doublearea(V,F)) {inf}
  %     'Quiet' followed by true or {false}
  % Outputs:
  %   VV  #VV by dim list of vertices over which TT is defined, should
  %     contain V as a prefix and should be exactly V unless F self-intersects
  %   TT  #TT by dim+1 list of tetrahedra indices
  %   FF  #FF by dim list of F potentially subdivided to be defined over VV,
  %     should be exactly F unless F self-intersects
  %   TN  #TT by dim+1 list of tetrahedra neighbors
  %   IFF  #FF by 1 list of original facets indices corresponding to each facet
  %     in FF
  %

  % parse input
  V = varargin{1};
  if nargin>=2
    F = varargin{2};
  else
    F = [];
  end
  % default values for optional input
  tetgen_flags = '';
  triangle_flags = '';
  use_bounding_box = 0;
  bb_push = 1.01;
  bb_ups = inf;
  quiet = 0;

  % parse any optional input
  ii = 3;
  while(ii <= nargin)
    switch varargin{ii}
    case 'TriangleFlags'
      ii = ii + 1;
      assert(ii<=nargin);
      triangle_flags = varargin{ii};
    case 'TetgenFlags'
      ii = ii + 1;
      assert(ii<=nargin);
      tetgen_flags = varargin{ii};
    case 'UseBoundingBox'
      ii = ii+1;
      assert(ii<=nargin);
      use_bounding_box = 1*(varargin{ii});
    case 'BoundingBoxPush'
      ii = ii+1;
      assert(ii<=nargin);
      bb_push = varargin{ii};
    case 'BoundingBoxUpsample'
      ii = ii+1;
      assert(ii<=nargin);
      bb_ups = varargin{ii};
    case 'Quiet'
      ii = ii+1;
      assert(ii<=nargin);
      quiet = varargin{ii};
      quiet = quiet * 1;
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  if quiet
    fid = fopen('/dev/null','w');
  else
    fid = 1;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Check for cached result, do NOT edit variables until cache is checked,
  % your function code comes later. See below
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear varargin;
  % get a list of current variables in this scope, this is the input "state"
  variables = who;
  % get a temporary file's name
  tmpf = [tempname('.') '.mat'];
  % save the "state" to file, so we can get a md5 checksum
  save(tmpf,'-regexp',sprintf('^%s$|',variables{:}),'-ascii');
  % get md5 checksum on input "state", we append .cache.mat to the check sum
  % because we'll use the checksum as the cache file name
  [s,cache_name] = system(['/sbin/md5 -r ' tmpf ' | awk ''{printf "."$1".cache.mat"}''']);
  % clean up
  delete(tmpf);
  clear s tmpf variables;
  % If the checksum cache file exists then we've seen this input "state"
  % before, load in cached output "state"
  if(exist(cache_name,'file'))
    fprintf(fid,'Cache exists. Using cache...\n');
    % use cache
    load(cache_name);
  % Otherwise this is the first time we've seen this input "state", so we
  % execute the function as usual and save the output "state" to the cache 
  else
    fprintf(fid,'First time. Creating cache...\n');


    dim = size(V,2);
    assert(numel(F) ==0 || dim >= size(F,2), sprintf( ...
      'dimension of V (%d) < simplex size (%d)',dim,size(F,2)));
    switch dim
    case 2
      % we're dealing with edges
      E = F;
      if use_bounding_box
        [BV,BF] = bounding_box([max(V)+eps;min(V)-eps]);
        dir = [sign(bsxfun(@minus,BV,mean(BV)))];
        mag_max = max(sqrt(sum(bsxfun(@minus,BV,mean(BV)).^2,2)));
        BV = bsxfun(@plus,bb_push*mag_max*dir,mean(BV));
        BE = [];
        EmBE = E;
      else
        BV = [];
        BF = [];
        % explicitly add convex hull
        F = delaunay(V);
        BE = outline(F);
        EmBE = setdiff(sort(E,2),sort(BE,2),'rows');
      end
      [VV,TT] = triangle([V;BV],[BE;EmBE;size(V,1)+BF],[],'Flags',triangle_flags);

      % This is a bullshit hack way of determining which edges were subdivided
      % A better way is using boundary markers in triangle like in
      % diffcurves_test.m
      if nargout>=3
      % First determine which edges were not subdivided
      E_indices = 1:size(E,1);
      [maintained] = ismember(sort(E,2),sort(edges(TT),2),'rows');
      EE = E(maintained,:);
      IEE = E_indices(maintained);
      % and those that were
      ES = E(~maintained,:);
      IES = E_indices(~maintained);
      % hopefully #ES is small, let Alec know if it isn't, we should rewrite this
      if size(ES,1)>100
        warning('many (%d) edges to subdivide, this could be slow',size(ES,1));
      end
      % build edge-length adjacency matrix
      C = adjacency_edge_cost_matrix(VV,edges(TT));
      % average edge length for tolerances
      h = avgedge(VV,TT);
      % loop over edges that were subdivided
      for e = 1:size(ES,1)
        % get shortest path from source to dest
        [D,P] = dijk(C,ES(e,2),ES(e,1));
        % check that D is what a straight line path should be

        orig_len = sqrt(sum((V(ES(e,2),:)-V(ES(e,1),:)).^2,2));
        if abs(orig_len-D)>1e-13*h
          warning('path %d to %d may not be straight: |%g - %g| = %g', ...
            ES(e,2),ES(e,1),orig_len,D,abs(orig_len-D));
        end
        d = ES(e,1);
        NE = [];
        % trace back path 
        while d~=ES(e,2)
          NE = [NE;d P(d)];
          d = P(d);
        end
        % append to edges
        EE = [EE;NE];
        IEE = [IEE(:);repmat(IES(e),size(NE,1),1)];
      end
      % rename to match output
      FF = EE;
      end
      if nargout >= 5
        IFF = IEE;
      end
      if nargout >= 4
        % Compute element neighbors
        TN = tt(TT);
      end
    case 3

      % It's really important that we have at least tetgen 1.5
      [status,result] = system([path_to_tetgen ' -help']);
      tetgen_version = ...
        str2double(regexprep(result,'.*Version ([0-9].[0-9]).*','$1'));
      if isnan(tetgen_version) || tetgen_version < 1.5
        error('%s version (%g) < 1.5',path_to_tetgen,tetgen_version);
      end

      FF = F;

      % Q: Should clean up be the caller's responsibility?
      % A: Yes, I think so. We can parse tetgen and give hints.

      %% number of input points
      %n = size(V,1);

      %% Deal with duplicate input points
      %[V,I,J] = remove_duplicate_vertices(V);
      %if n ~= size(V,1)
      %  warning('cdt: %d exactly duplicate points removed from input',n-size(V,1));
      %  n = size(V,1);
      %  % Remap faces
      %  FF = J(F);
      %end

      %if ~isempty(F)
      %  NF = F( (F(:,1) ~= F(:,2)) & (F(:,2) ~= F(:,3)) & (F(:,3) ~= F(:,1)),:);
      %  if size(NF,1) < size(F,1)
      %    warning('cdt: %d faces with combinatorally duplicate vertices removed\n', ...
      %      size(F,1)-size(NF,1));
      %    FF = NF;
      %  end
      %end

      % Try to mesh with all faces included directly
      BM = -ones(size(FF,1),1);

      if use_bounding_box
        [BV,BF] = bounding_box(V);
        % Subdivide
        if ~isempty(FF) && size(FF,2)==3
          max_dblA = max(doublearea(V,FF));
          while min(doublearea(BV,BF))> bb_ups*max_dblA
            [BV,BF] = upsample(BV,BF);
          end
        end
        % push a little a away
        BV = bsxfun(@plus,bb_push*bsxfun(@minus,BV,mean(BV)),mean(BV));
        BF = BF + size(V,1);
        V = [V;BV];
        FF = [FF;BF];
        % true boundary marker
        BM = [BM;ones(size(BF,1),1)];
      end

      prefix = tempname;
      poly_filename = [prefix '.poly'];
      writePOLY_tetgen(poly_filename,V,FF,[],'BoundaryMarkers',BM);


      % -p: mesh PLC
      % -C: check consistency of result
      % -c: retain convex hull, doesn't seem to play nicely with -q
      %% -O0: optimization level
      %% -Y: no steiner points on boundary
      %% -J: Do not jettison duplicate vertices (we should have already handled
      %%   these)
      % -n: output neighbor information
      %% -F: suppress output of face information
      %% -S0: specifies maximum number of added points (may still add schoenhardt
      %   points)
      % There seems to be a secret -o[val] flag that optimizes the max dihedral
      % angle
      delaunay_flags = '-pCn';
      if ~use_bounding_box
        % then retain convex hull. Note: this seems broken, hence the option to
        % use bounding box
        delaunay_flags = [delaunay_flags 'c'];
      end
      %delaunay_flags = '-pq0.75 -a0.0002 -O5YcJnFg';

      command = [path_to_tetgen ' ' delaunay_flags ' ' tetgen_flags ' ' poly_filename];
      fprintf(fid,'%s\n',command);
      if quiet
        [status, result] = system(command);
      else
        [status, result] = system(command,'-echo');
      end
      fprintf(fid,'%s',result);
      if status ~= 0
        fprintf(fid,'%s',result);
      else
        % Print some warnings based on tetgen output
        if strfind(result,'Steiner points on boundary edges')
          fprintf(fid,'%s',result);
          warning('Tetgen added steiner points on boundary');
        end
        if strfind(result,'Jettisoning redundants points.')
          fprintf(fid,'%s',result);
          warning('Tetgen removed duplicate points');
        end

      end

      % See terminatetetgen() in tetgen.h
      switch(status)
      case 0
      case 3
        error('Self intersections. Try selfintersect() or remove self-intersections');
      %% Self-intersection
      %  warning('Throwing away self-intersecting facets');
      %  % Seems tetgen sets status=3 if self-intersecting
      %  % Rerun and detect self-intersecting faces
      %  % -p: mesh PLC, -d: detect all self intersections, -J: don't jettison
      %  % duplicate points (since we want self-intersecting faces defined over
      %  % original vertices)
      %  selfintersection_flags = '-pdJ ';
      %  command = [path_to_tetgen ' ' selfintersection_flags ' ' poly_filename];
      %  fprintf(fid,command);
      %  [status, result] = system(command);
      %  if strfind(result,'No faces are intersecting.')
      %    fprintf(fid,'%s',result);
      %    error('Tetgen unsure whether there are self-intersections');
      %  end
      %  assert(status==0);
      %  % tetgen always writes output to file:
      %  face_filename = [prefix '.1.face'];
      %  node_filename = [prefix '.1.node'];
      %  % faces participating in self-intersections
      %  % Faces will index new nodes, but they better be the input ones
      %  IF = readFACE(face_filename);
      %  % clean up temp files
      %  delete(face_filename);
      %  %VV = readNODE(node_filename);
      %  %assert(size(VV,1)==size(V,1));
      %  %delete(node_filename);

      %  % Segments aren't supported by tetgen the way I'd like them to be (i.e. I
      %  % get seg faults)
      %  %% Get edges of intersecting faces
      %  %IFE = edges(IF);
      %  %IFsegs = mat2cell(IFE,ones(size(IFE,1),1),[2]);

      %  % get orignal faces \ intersecting faces 
      %  % Using sort like this assumes that tetgen ignores orientation
      %  FmIF = setdiff(sort(FF,2),sort(IF,2),'rows');
      %  Facets = [];
      %  Facets.facets = mat2cell(FmIF,ones(size(FmIF,1),1),[size(FmIF,2)]);
      %  Facets.boundary_marker = -ones(numel(Facets.facets),1);
      %  Facets.holes = cell(numel(Facets.facets),1);

      %  % Rewrite poly file
      %  writePOLY(poly_filename,V,Facets,[]);

      %  % shoot for the moon
      %  command = [path_to_tetgen ' ' delaunay_flags ' ' tetgen_flags ' ' poly_filename];
      %  fprintf(fid,command);
      %  [status, result] = system(command);
      %  result
      %  status
      case 4
      % Small feature
        error('Try again with reduced T?');
      case 5
      % close input facets
        error('Try again with -Y');
      case 6
      % close input facets
        error('Input error. Check yo` input');
      otherwise
        error('Tetgen returned status %d != 0',status);
      end

      % tetgen always writes output to file:
      %   xxxx.1.ele  tetrahedra
      %   xxxx.1.node tetrahedra vertices
      ele_filename = [prefix '.1.ele'];
      face_filename = [prefix '.1.face'];
      node_filename = [prefix '.1.node'];
      neigh_filename = [prefix '.1.neigh'];

      %delete(poly_filename);

      VV = readNODE(node_filename);
      TT = readELE(ele_filename);
      TN = readNEIGH(neigh_filename);
      TN(TN>=0) = TN(TN>=0)-1;


      [FF,BF] = readFACE(face_filename);
      %FF = FF-1;
      FF = FF(BF==-1,:);

      %[VV1,TT1] = readMESH(mesh_filename);
      % TODO: handle duplicate input vertices in a reasonable way
      % TODO: Tetgen seems to add steiners when input has coplanar neighboring
      % triangles (aka triangulated planar quads). It seems its enough to do:
      % V = V+eps*rand(size(V));
      if ~quiet
        if size(VV,1) > size(V,1)
          warning('New mesh has extra %d vertices',size(VV,1)-size(V,1));
        elseif size(VV,1) < size(V,1)
          warning('New mesh has %d fewer vertices',size(V,1)-size(VV,1));
        elseif max(abs(VV-V))>eps
          warning('max(abs(VV-V)) = %g > %g',max(abs(VV-V)),eps);
          %% Remap to original vertex list
          %TT = I(TT);
          %VV = V;
        end
      end
      %%% Reverse order so that boundary_faces(TT) meets CCW convention
      %TT = TT(:,[3 2 1 4]);
      %TN = TN(:,[3 2 1 4]);

      %% remove bounding box facets
      %if use_bounding_box
      %  old_FF_count = size(FF,1);
      %  FF = ...
      %    FF(all(reshape(~ismember(VV(FF(:),:),BV,'rows'),size(FF,1),3),2),:);
      %  rm_count = old_FF_count-size(FF,1);
      %  switch rm_count
      %  case 0:11
      %    warning('Only removed %d boundary facets (instead of 12)');
      %  case 12
      %    % alles gut
      %  otherwise
      %    warning('Removed %d boundary facets (instead of 12)',rm_count);
      %  end
      %end
    otherwise
      error(sprintf('Unsupported dimension (%d)',dim));
    end

    if nargout>=3
    [IN] = in_elements(FF,TT);
    if ~all(IN)
       warning('%d faces not in elements',sum(~IN));
    end
    end

    % get list of variables present in this scope at finish of function code,
    % this is the output "state"
    variables = who;
    % save output "state" to file, using md5 checksum cache file name
    save(cache_name,'-regexp',sprintf('^%s$|',variables{:}));

  end
end

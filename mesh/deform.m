classdef deform < handle
  % DEFORM a 2D mesh in a matlab plot window using various deformation methods.
  % This class must extend handle so that it may be updated by changes in the
  % plot windows it creates.
  %
  % this = deform(V,F,C,'ParameterName',ParameterValue',...)
  %
  % Inputs:
  %  V  #V by 2 list of vertex positions
  %  F  #F by 3 list of face indices
  %  C  #C by 2 list of control point positions
  %  OPTIONS:
  %    'DeformMethod' followed by deformation method:
  %       'skinning' followed by W, a #V by #handles matrix of weights
  %       'mls' followed by mlsd, a struct with moving least squares data
  %       'arap' (see 'ArapEnergy' for more options)
  %       'takeo-asap'
  %       'lscm'
  %       'takeo-arap'
  %       'takeo-arap-iterative'
  %       'iterative-laplacian-mesh-editing'
  %       'laplacian-mesh-editing'
  %       'custom' followed by a function that takes the mesh and the boundary
  %         conditions and outputs the new vertex positions. Dummy identity
  %         function:
  %           f = @(V,F,b,bc) V;
  %    'PointHandles' followed by a list of indices into C, for the point
  %      handles. Default value is 1:size(C,1)
  %    'BoneEdges' followed by a #BE by 2 list of indices into C for the bone
  %      handles. Default value is []
  %    'CageEdges' followed by a #CE by 2 list of indices into P specifying 
  %      cage edges for point handles (for display purposes only)
  %    'InterpMode' could be 'LBS' for Linear Blend Skinning or 'DQLBS' for
  %      Dual Quaternion Linear Blend Skinning, default is 'LBS'
  %    'StretchBones' followed by a set of weights #V by #BE. Use stretchable
  %      bones skinning formula instead of LBS, uses extra set of weights per
  %      bone.
  %    'StiffPoints' followed by a list of bone weights #V by #BE, expects
  %      'BoneEdges' to be present and to index into P rather than C. Uses
  %      extra set of weights per bone edge.
  %    'AutoDOF' followed by a method to compute automatic degrees of freedom
  %      for point handles, possible options are:
  %       'none'
  %       'arap'
  %       'pseudoedges'
  %    'PseudoEdges' followed by a #PE by 2 list of pseudo edges for computing
  %      automatic rotations at point handles
  %    'Groups'  followed by G a #V list of groups ids, or a number of groups
  %      to compute based on weights
  %    'Free'
  %      List of "free" handles, that is, LBS handles that are not controled
  %      by the user
  %    'Fixed'
  %      List of "fixed" handles, that is, LBS handles that are *completely*
  %      controlled by the user (AutoDOF methods should not affect these
  %      handles). Note that "fixed" is not necessarily the complement of
  %      "free"
  %    'UseLastFrame' followed by true of false, default is true
  %    'VisualizeWeights' puts a subplot to the right with the selected
  %      handle's weights in a plot using color to indicate value, default is 
  %      off
  %    'VisualizeRotations' puts a subplot to the right with a
  %      visualization of the rotations fit to the vertices in the current
  %      deformation, default is off
  %    'WeightVisualizationWeights' followed by an array the same size as W,
  %      default is to just use W
  %    'ShowContour'  puts a subplot to the right with the selected handle's
  %      weights in a contour plot, default is off
  %    'ContourWeights' followed by an array the same size as W
  %      default is to just use W
  %    'ContourLevels' followed by the number of isolevels in the contour plot,
  %      default is 20
  %    'ShowPixelDiff' followed by another set of weights show the difference
  %      between the current deformation using the current weights and a
  %      deformation by this other set of weights.
  %    'ColorByGroups'  followed by true or false, default is false
  %    'ColorByRotations'  followed by true or false, default is false
  %    'DepthOrdering' list of handle indices in order of depth, default is
  %      1:nh
  %    'Tol'  stopping tolerance for arap optimization, see arap.m
  %    'MaxIter'
  %       max number of local-global iterations used in arap_dof, default is 10
  %    'Image' deform an image rather than mesh (V,F are ignored). Followed by:
  %       im  w by h by {1|3} image
  %       X  x-coordinates of grid over image
  %       Y  y-coordinates of grid over image
  %    'ColorScheme' 
  %      'BBW' Use color scheme similar to [Jacobson et al 2011]
  %      'MLS' Use Color scheme similar to "Image Deformation Using Moving
  %        Least Squares", by [Schaefer et al] 
  %    'Tets' followed by a #T by 4 list of tet indices
  %    'ArapEnergy'
  %      followed by a string specifying which arap energy definition to use.
  %      One of the following:
  %        'spokes'  "As-rigid-as-possible Surface Modeling" by [Sorkine and
  %          Alexa 2007], rotations defined at vertices affecting incident
  %          edges
  %        'elements'  "A local-global approach to mesh parameterization" by
  %          [Liu et al.  2010] or "A simple geometric model for elastic
  %          deformation" by [Chao et al.  2010], rotations defined at
  %          elements (triangles or tets) 
  %        'spokes-and-rims'  Adapted version of "As-rigid-as-possible Surface
  %          Modeling" by [Sorkine and Alexa 2007] presented in section 4.2 of
  %          or "A simple geometric model for elastic deformation" by [Chao et
  %          al.  2010], rotations defined at vertices affecting incident
  %          edges and opposite edges
  %    'Stripes' followed by number of strips and dimension to use for stripes
  %    'ControlGroups' followed by #C by 1 list of group ids 
  %      group control handles, so that when one in a group is selected and
  %      deformed the deformation is applied to all in the group
  % Output:
  %   this  object giving access to plot handles and input variables of
  %     deformation and deformation plots
  %
  % Example:
  %   % given a texture (im,alpha) use a deform class instance (d) to show a
  %   % deformed image
  %   [oim, oalpha] = ...
  %     texture_map([d.new_V d.W(:,2)],d.F,texture_coords(d.V),im,alpha);
  %   
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read only class fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties(GetAccess=public,SetAccess=protected)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Domain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of mesh vertices
    n
    % number of control vertices
    c
    % mesh vertex positions
    V
    % mesh connectivity
    F
    % tet mesh connectivity
    Tets = [];
    % control point positions
    C
    % indices of closest mesh vertices to C
    b
    % point handles, indices into C
    P
    % bone handles, #B by 2 list of bone endpoints, indices into C
    BE = [];
    % cage edges, #CE by 2 list of cage edge endpoints, indices into P
    CE = [];
    % image
    im
    % grid dimensions
    grid_h
    grid_w
    % flag whether using grid
    grid = false;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deformation 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % method used for deformation updates
    deform_method = 'skinning';
    % custom deform function handle
    custom_deform
    % energy for solving arap, default is Sorkine and Alexa
    arap_energy = 'spokes';
    % pose rotations of point handles
    R  
    % pose positions of control vertices
    new_C  
    % pose mesh vertex positions
    new_V
    % plot handle to scatter plot of interpolated handles
    ch  
    %% plot handle to scatter plot of free handles
    %fch  
    % Pose transformatiosn of handles
    TR
    % moving least squares data
    mlsd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skinning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skinning weights
    W = [];
    % Skinning transformations
    L = [];
    % interpolation mode:
    %  LBS
    %  DQLBS
    interp_mode = 'LBS';
    % use stretchable, twistable bones deformation formula
    stretch_bones = false;
    % endpoint weights for stretch, twistable bones skinning, or incident bone
    % weights for stiff points skinning
    WBE = [];
    % use stiff points deformation formula
    stiff_points = false;
    % use automatically computed degrees of freedom at point handles
    auto_dof = 'none';
    % pseudo-edges for auto_dof
    PE = [];
    % use paritioning of mesh into groups to speed up auto_dof
    G = [];
    % list of "free" point handles, indices into P
    free = [];
    % list of "fixed" point handles, indices into P
    fixed = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot handle to main triangle mesh plot
    tsh  
    % handle to grid warp plot
    gsh
    % plot handle to weight visualization plot
    wvsh  
    % plot handle to rotation visualization plot
    rvsh  
    % plot handle to pixel difference visulation plot
    wpdsh  
    % show weight visualization off to right
    visualize_weights = false;
    % show rotation visualization off to right
    visualize_rotations = false;
    % weights used for weight visualization
    WVW
    % show contour plot off to right 
    show_contour = false;
    % weights used for contour plots
    CW
    % number of levels in contour plot
    contour_levels = 10;
    % show pixel difference of deformation with W and W_other off to right
    show_pixel_diff = false;
    % other weights used to compute pixel difference
    W_other = [];
    % color mesh according to partitioning groups
    color_by_groups = false;
    % color mesh according to rotations fit at each vertex cell
    color_by_rotations = false;
    % depth order
    depth_ordering
    % color scheme
    color_scheme = 'BBW'
    % stripes
    stripe_num = 0
    stripe_dim = 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slaves
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % list of slaves
    slaves = [];
    % internal flag used for controlling slaves
    slave_seen = true;
    % main plot axis handle
    ah  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interaction 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag for whether controls are currently being clicked down (dragged) on
    is_down = false;
    % keep track of down position
    down_pos = [];
    % vector from camera plane to axes front
    cp2cf = [];
    % keep track of last two drag positions
    drag_pos = [];
    last_drag_pos = [];
    % keep track of index of selected control point
    ci;
    % Control groups
    CG = [];
    % select fixed vertices in groups
    flood_fill_sel;
    % type of click ('left','right')
    down_type  = '';
    win_max = [];
    win_min = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Other
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_iterations = 10;
    stretch_ratio = 0;
    tol = 0.01;
    % use last frame to initialize deformation updates
    use_last_frame = true;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Private class fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties(Access=private)
    % Precomputation for arap automatic degrees of freedom
    Lall = [];
    % number of subplots in figure window
    number_of_subplots;
    % number of point handles
    np
    % number of bone handles
    nb
    % extra column for each mesh vertex position: depth
    Zdepth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handles to plot elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B_plot_outer
    B_plot_inner
    CE_plot_outer
    CE_plot_inner
    PE_plot
    Wr
    contour_handle
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Public Class methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods(Access=public)

    function this = deform(varargin)
    % See also: deform

      % parse mandatory input
      [this,remaining_inputs] = this.parse_mandatory_inputs(varargin);
      % Parse additional options
      this = this.parse_optional_inputs(remaining_inputs);
      % Set up plots
      this = this.init_plots();
      % set up axes
      this = this.set_axes();
      % Prompt user
      if(this.visualize_weights || this.show_contour)
        fprintf(['\nCLICK a control point to visualize its corresponding ' ...
          'weights.\n']);
      end
      fprintf( ...
        ['DRAG a control point to deform the mesh.\n' ...
        'RIGHT CLICK DRAG a control point to rotate point handles.\n\n']);
      return
    end

    function free_bone_endpoint_indicesC = free_bone_endpoint_indices(this)
      % deform.FREE_BONE_ENDPOINT_INDICES return set of indices into C that are
      % free bone end points
      %
      all_BE = unique(this.BE(:))';
      free_bone_endpoint_indicesC = all_BE(ismember(all_BE,this.free));
    end


    function interpolatedC = interpolated(this)
      % deform.INTERPOLATED return set of indices into C that are NOT in free
      %
      all = 1:size(this.C,1);
      interpolatedC = all(~ismember(all,this.free));
    end

    function this = set_axes(this)
      % deform.SET_AXES enlarges the main plot window's axes to fit
      % current mesh and control points
      %

      % don't set axes if dimension are greater than 2

      if(this.grid)
        X = get(this.gsh,'XData');
        Y = get(this.gsh,'YData');
        cur_V = [X(:) Y(:)];
      else
        cur_V = get(this.tsh,'Vertices');
      end
      cur_V = cur_V(:,1:size(this.V,2));

      % scale window to fit
      this.win_min = min([this.win_min; ...
        cur_V; ...
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],:)]);
      this.win_max = max([this.win_max; ...
        cur_V; ...
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],:)]);
      this.win_max(this.win_min==this.win_max) = this.win_max(this.win_min==this.win_max)+eps;
      [az,el]= view(this.ah);
      axis(this.ah,reshape([this.win_min;this.win_max],1,2*size(cur_V,2)));
      view(this.ah,az,el);

      % update slave deform objects
      if this.is_down
        % BFS to mark all reachable slaves
        slaves = this.slaves;
        slaves_Q = this.slaves;
        while numel(slaves_Q) > 0
          that = slaves_Q(end);
          slaves_Q = slaves_Q(1:(end-1));
          if (this ~= that) && that.slave_seen
            that.slave_seen = false;
            slaves = [slaves that.slaves];
            slaves_Q = [slaves_Q that.slaves];
          end
        end
        % be sure there are no repeats
        slaves = setdiff(unique(slaves),[this]);
        if ~isempty(slaves)
          for that = slaves
            that.slave_seen = true;
            %assert(this ~= that);
            if ~that.is_down
              that.set_axes();
            end
          end
        end
      end
    end

    function this = append_slaves(this,thats)
      % APPEND_SLAVES  append a deform class as a slave to this class. When a
      % user changes controls of this class the slave classes will
      % automatically follow
      %
      % this = append_slaves(this,thats)
      %
      % Inputs:
      %   this  master deform object 
      %   thats list of slave deform objects
      % Outputs:
      %   this  copy of input this with updated slaves
      %
      this.slaves = [this.slaves thats];
    end

    function this = set_new_C(this,new_C)
      this.new_C = new_C;
      this.update_positions();
    end

    function this = update_positions(this)
      % update mesh positions
      %
      title(this.ah,'Updating...');
      drawnow

      % Set initial guess to current deformed vertex positions
      guess = [];
      if(this.use_last_frame)
        guess = this.new_V;
      end

      % Force one-ring around each control point to translate with control point
      b1 = this.b;
      bc1 = this.new_C;
      clamp_one_rings = false;
      if clamp_one_rings
        indices = 1:size(this.V,1);
        % clamp to left and right of two point handles
        left=indices(this.V(:,1)<this.C(1,1));
        right=indices(this.V(:,1)>this.C(2,1));
        b1 = [ left right];
        bc1 = [...
          this.V(left,:) + ...
          repmat(this.new_C(1,:)-this.C(1,:),numel(left),1); ...
          this.V(right,:) + ...
          repmat(this.new_C(2,:)-this.C(2,:),numel(right),1)];
        %A = adjacency_list(this.F,this.b);
        %for ii = 1:size(A,1)
        %  b1 = [b1 A{ii}];
        %  bc1 = [bc1; ...
        %    this.V(A{ii},:) + ...
        %    repmat(this.new_C(ii,:) - ...
        %    this.V(this.b(ii),:),size(A{ii},2),1)];
        %end
      end

      switch this.deform_method
      case 'mls'
        this.new_V = MLSD2DTransform(this.mlsd,this.new_C')';
      case 'arap'
        if isempty(this.Tets)
          this.new_V = arap( ...
            this.V,this.F,b1,bc1, ...
            'V0',guess,'MaxIter',this.max_iterations, ...
            'Energy',this.arap_energy, ...
            'Groups', this.G,'Tol',this.tol);
        else
          this.new_V = arap( ...
            this.V,this.Tets,b1,bc1, ...
            'V0',guess,'MaxIter',this.max_iterations, ...
            'Energy',this.arap_energy, ...
            'Groups', this.G,'Tol',this.tol);
        end
      case 'takeo-arap-iterative'
        this.new_V = takeo_arap( ...
          this.V,this.F,b1,bc1, ...
          'V0',guess,'MaxIter',this.max_iterations, ...
          'Groups', this.G,'Tol',this.tol);
      case 'takeo-arap'
        this.new_V = takeo_arap(this.V,this.F,b1,bc1);
      case 'takeo-asap'
        this.new_V = takeo_asap(this.V,this.F,b1,bc1);
      case 'lscm'
        this.new_V = lscm(this.V,this.F,b1,bc1);
      case 'laplacian-mesh-editing'
        this.new_V = ...
          iterative_laplacian_mesh_editing(this.V,this.F,b1,bc1);
      case 'iterative-laplacian-mesh-editing'
        this.new_V = iterative_laplacian_mesh_editing( ...
          this.V,this.F,b1,bc1, ...
          'V0',guess,'MaxIter',this.max_iterations, 'Tol',this.tol);
      case 'skinning'
        % update mesh positions into this.new_V
        assert(~clamp_one_rings);
        this = skinning_update(this);
      case 'custom'
        this.new_V = this.custom_deform(this.V,this.F,b1,bc1);
      otherwise
        error(['''' this.deform_method '''is not a valid deformation method.']);
      end

      % replace nan entries with original entries
      this.new_V(isnan(this.new_V)) = this.V(isnan(this.new_V));
      %this.new_V = this.V;

      % update plots with this.new_V
      this = this.update_display_positions();
      % update subplot visualizations
      this = this.update_visualizations();

      % update slaves
      this = this.update_slaves();
      title(this.ah,'');

    end

  end% public methods
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Strictly private methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MATLAB <2011 needs this to be public
  %methods(Access=private)
  methods(Access=public)
    function [this,remaining_inputs] = parse_mandatory_inputs(this,inputs)
      % parse mandatory inputs
      %

      % vertex positions of mesh
      this.V = inputs{1};
      % face indices of mesh
      this.F = inputs{2};
      % control vertices of skeleton
      this.C = inputs{3};
      % pass on remaining inputs
      remaining_inputs = inputs(4:end);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set defaults for parameters which depend on input
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % set default point handles
      this.P = 1:size(this.C,1);
      % Be sure that control vertices are in 2D
      if(size(this.C,2) == 3 && size(this.V,2) == 2)
        warning('Discarding z-coordinates in C');
        this.C = this.C(:,1:2);
      end
      % control groups
      this.CG = 1:size(this.C,1);
      % main plot axis dimensions, min and mix
      win_min = repmat(inf,1,size(this.V,2));
      win_max = repmat(-inf,1,size(this.V,2));
    end

    function this = parse_optional_inputs(this,remaining_inputs)
      % parse remaining optional inputs
      %

      ii = 1;
      % number of remaining inputs
      num_ri = numel(remaining_inputs);
      while(ii <= numel(remaining_inputs))
        switch remaining_inputs{ii}
        case 'DeformMethod'
          ii = ii + 1;
          assert(ii<=num_ri);
          assert(ischar(remaining_inputs{ii}));
          this.deform_method = remaining_inputs{ii};
          % some deformation methods are followed by extra arguments
          switch this.deform_method
          case 'skinning'
            ii = ii + 1;
            assert(ii<=num_ri);
            assert(isfloat(remaining_inputs{ii}));
            this.W = remaining_inputs{ii};
            % Weights used for contours
            if isempty(this.CW)
              this.CW = this.W;
            end
            % Weights used for weight visualization
            if isempty(this.WVW)
              this.WVW = this.W;
            end
          case 'mls'
            ii = ii + 1;
            assert(ii<=num_ri);
            assert(isstruct(remaining_inputs{ii}));
            this.mlsd = remaining_inputs{ii};
            % Weights used for contours
            if isempty(this.CW)
              this.CW = this.mlsd.w';
            end
            % Weights used for weight visualization
            if isempty(this.WVW)
              this.WVW = this.mlsd.w';
            end
          case 'custom'
            ii = ii + 1;
            assert(ii<=num_ri);
            assert(isa(remaining_inputs{ii}, 'function_handle'));
            this.custom_deform = remaining_inputs{ii};
          end
        case 'PointHandles'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.P = remaining_inputs{ii};
        case 'BoneEdges'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.BE = remaining_inputs{ii};
        case 'CageEdges'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.CE = remaining_inputs{ii};
        case 'InterpMode'
          ii = ii + 1;
          assert(ii<=num_ri);
          if( strcmp(remaining_inputs{ii},'LBS')) 
            this.interp_mode = 'LBS';
          elseif( strcmp(remaining_inputs{ii},'DQLBS')) 
            this.interp_mode = 'DQLBS';
          else
            error('InterpMode must be either LBS or DQLBS');
          end
        case 'StretchBones'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.WBE = remaining_inputs{ii};
          this.stretch_bones = true;
        case 'StiffPoints'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.WBE = remaining_inputs{ii};
          this.stiff_points= true;
        case 'AutoDOF'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.auto_dof = remaining_inputs{ii};
        case 'PseudoEdges'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.PE = remaining_inputs{ii};
        case 'Groups'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.G = remaining_inputs{ii};
        case 'Free'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.free = remaining_inputs{ii};
        case 'Fixed'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.fixed = remaining_inputs{ii};
        case 'UseLastFrame'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.use_last_frame = remaining_inputs{ii};
        case 'VisualizeWeights'
          this.visualize_weights = true;
        case 'VisualizeRotations'
          this.visualize_rotations = true;
        case 'WeightVisualizationWeights'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.WVW = remaining_inputs{ii};
        case 'ShowContour'
          this.show_contour = true;
        case 'ContourLevels'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.contour_levels = remaining_inputs{ii};
        case 'ContourWeights'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.CW = remaining_inputs{ii};
        case 'ShowPixelDiff'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.Wother = remaining_inputs{ii};
          this.show_pixel_diff = true;
        case 'ColorByGroups'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.color_by_groups = remaining_inputs{ii};
        case 'ColorByRotations'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.color_by_rotations = remaining_inputs{ii};
        case 'DepthOrdering'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.depth_ordering = remaining_inputs{ii};
        case 'Tol'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.tol = remaining_inputs{ii};
        case 'MaxIter'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.max_iterations = remaining_inputs{ii};
        %case 'StretchRatio'
        %  ii = ii + 1;
        %  assert(ii<=num_ri);
        %  this.stretch_ratio = remaining_inputs{ii};
        case 'Image'
          this.grid = true;
          ii = ii + 1;
          assert(ii<=num_ri);
          this.im = remaining_inputs{ii};
          ii = ii + 1;
          assert(ii<=num_ri);
          X = remaining_inputs{ii};
          ii = ii + 1;
          assert(ii<=num_ri);
          Y = remaining_inputs{ii};
          this.grid_h = size(X,1);
          this.grid_w = size(X,2);
          % positions of vertices in grid
          this.V = [X(:) Y(:)];
        case 'ColorScheme'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.color_scheme = remaining_inputs{ii};
        case 'Tets'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.Tets = remaining_inputs{ii};
        case 'ArapEnergy'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.arap_energy = remaining_inputs{ii};
        case 'Stripes'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.stripe_num = remaining_inputs{ii};
          ii = ii + 1;
          assert(ii<=num_ri);
          this.stripe_dim = remaining_inputs{ii};
        case 'ControlGroups'
          ii = ii + 1;
          assert(ii<=num_ri);
          this.CG = remaining_inputs{ii};
        otherwise
          error(['''' remaining_inputs{ii} ''' is not a valid parameter']);
        end
        ii = ii + 1;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set values that may depend on extra parameters
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % checking parameters based on deformation method given
      switch this.deform_method
        case 'skinning'
          if isempty(this.W)
            % use shepard weights as a default
            warning(['''skinning'' chosen as deformation method but no ' ...
              'weights provided. Using shepard weights...']);
          end
          m = numel(this.P) + size(this.BE,1);
          if size(this.W,2) < (numel(this.P) + size(this.BE,1))
            error('Too few weights given');
          elseif (size(this.W,2) > m)
            warning('Too many weights given, assuming remaing weights are ''free''');
            % put all point handles positions at front of C and append zeros
            % for extra weights
            old_nc = size(this.C,1);
            this.C = [ this.C; zeros(size(this.W,2)-m,size(this.C,2))];
            % add extra weights to free handles
            this.free = unique([this.free old_nc+(1:(size(this.W,2)-m))]);
            % rearrange extra weights to come in between real point handles and
            % bones
            this.W = this.W(:, [ ...
              1:numel(this.P) ... % point weights
              (m+1):size(this.W,2) ... % extra "free" weights
              (numel(this.P)+1):m ... % bone weights
              ]);
            % remap fixed indices to new locations
            this.fixed(this.fixed > numel(this.P)) = ...
              this.fixed(this.fixed > numel(this.P)) + (size(this.W,2)-m);
            % point handles for each control position
            this.P = [this.P old_nc+(1:(size(this.W,2)-m))];
          else
            assert(size(this.W,2) == (numel(this.P) + size(this.BE,1)));
          end
        case 'takeo-arap'
          % force max iterations to be 1
          this.max_iterations = 1;
      end

      % if endpoint weights exist then there better be bones and the right
      % number of them
      assert(isempty(this.WBE) || size(this.WBE,2) == size(this.BE,1))
      % number of point handles
      this.np = numel(this.P);
      % number of bone handles
      this.nb = size(this.BE,1);
      % pose positions of control points
      this.new_C = this.C;
      % rotations (angles) stored at each control point
      switch size(this.V,2)
      case 2
        this.R = zeros(this.np,1);
      case 3
        this.R = repmat(eye(3,3),[1 1 this.np]);
      end

      if isempty(this.W)
        this.W = shepard(this.V,this.C,this.P,this.BE,[]);
      end

      if size(this.V,2) == 2
        % depth ordering
        this.depth_ordering = 1:(size(this.W,2));
      end

      if numel(this.G) == 1
        % given number of groups rather than group ids
        this.G = partition(this.W,this.G);
      end

      % number of mesh vertices
      this.n = size(this.V, 1);
      % number of control vertices
      this.c = size(this.C,1);
      assert(size(this.V,2) == size(this.C,2));
      % compute distance from every vertex in the mesh to every control vertex
      D = permute(sum((repmat(this.V,[1,1,this.c]) - ...
        permute(repmat(this.C,[1,1,this.n]),[3,2,1])).^2,2),[1,3,2]);
      % use distances to determine closest mesh vertex to each control vertex
      % Cv(i) is closest vertex in V to ith control vertex in C
      [XXX,this.b] = min(D);

      % if number of unique closest mesh vertices is less than total number,
      % then we have contradictory boundary conditions
      if(~all(size(unique(this.b)) == size(this.b)))
        warning('Multiple control vertices snapped to the same domain vertex');
      end

      % if using last frame then prepare new_V with rest positions
      if(this.use_last_frame)
        this.new_V = this.V;
      end
    end

    function this = init_plots(this)
      % Initialize main plot and subplots as needed
      %

      % Set up figure with enough subplots for any additional visualizations
      % clear current figure
      clf
      % compute number of subplot depending on 
      this.number_of_subplots = 1 + ...
        this.visualize_weights + ...
        this.visualize_rotations + ...
        this.show_contour + ...
        this.show_pixel_diff;
      current_subplot = 1;
      if(this.number_of_subplots>1)
        % only use subplots if we have more than one
        subplot(1,this.number_of_subplots,1);
      end
      this = this.init_main_plot();

      % Contour plot
      if(this.show_contour)
        assert(~this.grid);
        current_subplot = current_subplot + 1;
        assert(all(size(this.CW) == size(this.W)));
        % Put contour plot and isolines in other subplot
        subplot(1,this.number_of_subplots,current_subplot);
        this = this.init_contour_plot();
      end
      % pixel diff plot
      if(this.show_pixel_diff)
        assert(~this.grid);
        current_subplot = current_subplot + 1;
        % subplot for weight visualization
        subplot(1,this.number_of_subplots,current_subplot);
        assert(all(size(this.Wother) == size(this.W)));
        hold on;
        % plot rest mesh
        this.wpdsh = ...
          trisurf(this.F,this.V(:,1),this.V(:,2),zeros(size(this.V,1),1), ...
          'FaceColor','interp');
        colorbar;
        view(2);
        axis equal
        hold off;
      end
      % weight visualization plot
      if(this.visualize_weights)
        assert(~this.grid);
        current_subplot = current_subplot + 1;
        % subplot for weight visualization
        subplot(1,this.number_of_subplots,current_subplot);
        assert(all(size(this.CW) == size(this.W)));
        hold on;
        % plot the original mesh
        this.wvsh = ...
          trisurf(this.F,this.V(:,1),this.V(:,2),zeros(size(this.V,1),1), ...
          'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
        colormap(get(this.wvsh,'Parent'),jet(10));
        view(2);
        axis equal
        hold off;
      end
      % rotation visualization plot
      if(this.visualize_rotations)
        assert(~this.grid);
        current_subplot = current_subplot + 1;
        % subplot for weight visualization
        subplot(1,this.number_of_subplots,current_subplot);
        hold on;
        % plot the original mesh
        this.rvsh = ...
          trisurf(this.F,this.V(:,1),this.V(:,2),zeros(size(this.V,1),1), ...
          'FaceColor','interp','EdgeColor','interp');
        % initial values: 0
        set(this.tsh,'CData',zeros(size(this.V,1),1));
        % set color span
        caxis([-pi/2 pi/2]);
        view(2);
        axis equal
        axis manual
        hold off;
      end
    end %init_plots

    function this = init_main_plot(this)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Main plot
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % compute depth offsets based on weights
      depth_offset = repmat(linspace(0,-1e5,size(this.W,2)),size(this.V,1),1);
      if size(this.V,2) == 2
        this.Zdepth = sum(depth_offset(:,this.depth_ordering).*this.W,2);
        this.Zdepth(isnan(this.Zdepth)) = -1e6;
      end
      % save axis handle
      this.ah = gca;

      % plot the original mesh
      if(this.grid)
        this.gsh = warp( ...
          reshape(this.V(:,1),this.grid_h,this.grid_w), ...
          reshape(this.V(:,2),this.grid_h,this.grid_w), ...
          reshape(this.Zdepth,this.grid_h,this.grid_w), ...
          this.im);
      else
        if size(this.V,2) == 2
          this.tsh = trisurf(this.F,this.V(:,1),this.V(:,2),this.Zdepth, ...
            'FaceColor',[0.5 1 0.5]);
        else 
          this.tsh = trisurf(this.F,this.V(:,1),this.V(:,2),this.V(:,3), ...
            'FaceColor',[0.5 1 0.5]);
        end
        % clicking on mesh counts as clicking on handle
        click_on_mesh = true;
        if click_on_mesh
          set(this.tsh,'ButtonDownFcn',@this.oncontrolsdown);
        end
        if(this.color_by_groups)
          set(this.tsh,'EdgeColor','interp');
          set(this.tsh,'FaceColor','interp');
          if ~isempty(this.G)
            set(this.tsh,'CData',this.G);
          else
            % if groups are empty then use a "group for each vertex"
            [s,randC] = sort(rand(size(this.V,1),1));
            set(this.tsh,'CData',randC);
          end
        elseif(this.color_by_rotations)
          set(this.tsh,'EdgeColor','interp');
          set(this.tsh,'FaceColor','interp');
          set(this.tsh,'CData',zeros(size(this.V,1),1));
          caxis(this.ah,[-pi/2 pi/2]);
        end
        if this.stripe_num > 0
          if this.visualize_weights
            warning('Stripes will not appear while visualizing weights');
          end
          set(this.tsh, ...
            'FaceColor','interp','FaceLighting','phong', ...
            'CData',this.V(:,this.stripe_dim));
          colormap(this.ah,flag(this.stripe_num));
        end
      end
      % 2D view
      view(2);
      axis equal
      % matlab's depth ordering is broken, this forces handles to display on
      % top of grid
      if(this.grid)
        axis manual
      end
      % plot bones
      hold on;
      if(this.nb > 0)
        % plot thick lines for bones (outline of lines)
        this.B_plot_outer = plot( ...
          [this.C(this.BE(:,1),1) this.C(this.BE(:,2),1)]', ...
          [this.C(this.BE(:,1),2) this.C(this.BE(:,2),2)]', ...
          '-k', ...
          'LineWidth',5);
        % plot thin lines for bones (innerline of lines)
        this.B_plot_inner = plot( ...
          [this.C(this.BE(:,1),1) this.C(this.BE(:,2),1)]', ...
          [this.C(this.BE(:,1),2) this.C(this.BE(:,2),2)]', ...
          '-b', ...
          'LineWidth',2);
      end
      if(size(this.PE,1) > 0)
        % plot lines for pseudo edges
        this.PE_plot = plot( ...
          [this.C(this.P(this.PE(:,1)),1) this.C(this.P(this.PE(:,2)),1)]', ...
          [this.C(this.P(this.PE(:,1)),2) this.C(this.P(this.PE(:,2)),2)]', ...
          '--r', ...
          'LineWidth',5);
      end
      if(size(this.CE,1) > 0)
        % plot lines for cage edges
        this.CE_plot_outer = plot( ...
          [this.C(this.P(this.CE(:,1)),1) this.C(this.P(this.CE(:,2)),1)]', ...
          [this.C(this.P(this.CE(:,1)),2) this.C(this.P(this.CE(:,2)),2)]', ...
          '-k', ...
          'LineWidth',5);
        this.CE_plot_inner = plot( ...
          [this.C(this.P(this.CE(:,1)),1) this.C(this.P(this.CE(:,2)),1)]', ...
          [this.C(this.P(this.CE(:,1)),2) this.C(this.P(this.CE(:,2)),2)]', ...
          '-', ...
          'Color', [1 1 0.2], ...
          'LineWidth',2);
      end

      % plot the control points (use 3D plot and fake a depth offset by pushing
      % control points up in z-direction)
      if size(this.V,2) == 2
        this.ch = scatter3( ...
          this.C([this.interpolated() this.free_bone_endpoint_indices()],1), ...
          this.C([this.interpolated() this.free_bone_endpoint_indices()],2), ...
          0.1+0*this.C([this.interpolated() this.free_bone_endpoint_indices()],1), ... 
          'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
          'LineWidth',2,'SizeData',500, ...
          'ButtonDownFcn',@this.oncontrolsdown);

      else
        this.ch = scatter3( ...
          this.C([this.interpolated() this.free_bone_endpoint_indices()],1), ...
          this.C([this.interpolated() this.free_bone_endpoint_indices()],2), ...
          this.C([this.interpolated() this.free_bone_endpoint_indices()],3), ... 
          'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
          'LineWidth',2,'SizeData',100, ...
          'ButtonDownFcn',@this.oncontrolsdown);
      end
      %this.fch = scatter3( ...
      %  this.C(this.free,1),this.C(this.free,2),0.1+0*this.C(this.free,1), ... 
      %  'o','MarkerFaceColor',[0.95 0.9 0.5], ...
      %  'MarkerEdgeColor',[0.2 0.2 0.2], ...
      %  'LineWidth',2,'SizeData',100, ...
      %  'ButtonDownFcn',@this.oncontrolsdown);
      hold off;

      switch this.color_scheme
      case 'MLS'
        set(this.ch, ...
        'Marker','s', ...
        'MarkerFaceColor','b','MarkerEdgeColor','b', ...
        'LineWidth',2,'SizeData',100);
        set(gca, 'visible', 'off');
        set(gcf,'Color','white');
      case 'BBW'
        %do nothing
      otherwise
        error(['''' this.color_scheme '''is not a valid color scheme.']);
      end
      if size(this.V,2) == 3
        axis manual;
      end
    end

    function this = init_contour_plot(this)
      % precompute contour information for each handle
      this.Wr = [];
      fprintf('Precomputing contour information:\n');
      for(ii = 1:(this.np+this.nb))
        fprintf( ...
          [ '  for handle ' num2str(ii) ...
            ' out of ' num2str(this.np+this.nb) '...\n']);
        [Xr,Yr,Wrii] = triinterp(this.V,this.F,this.CW(:,ii),true);
        if(isempty(this.Wr))
          this.Wr = zeros([size(Wrii) (this.np+this.nb)]);
        end
        this.Wr(:,:,ii) = Wrii;
      end
      % display contour
      hold on;
      [XXX,this.contour_handle] = contourf(Xr,Yr,this.Wr(:,:,1));
      this.contour_levels = 20;
      set(this.contour_handle,'LevelStepMode','Manual');
      set(this.contour_handle, ...
        'LevelStep',(max(this.CW(:))-min(this.CW(:)))/this.contour_levels);
      axis equal
      axis manual
      hold off;
      % Find all edges in mesh, note internal edges are repeated
      E = sort( ...
        [ this.F(:,1) this.F(:,2); ...
          this.F(:,2) this.F(:,3); ...
          this.F(:,3) this.F(:,1)],2);
      % determine uniqueness of edges
      [u,m,n] = unique(E,'rows');
      % determine counts for each unique edge
      counts = accumarray(n(:), 1);
      % extract edges that only occurred once
      O = u(counts==1);
      % display outline of mesh
      line(this.V([O(:);O(1)],1),this.V([O(:);O(1)],2),'Color','k','LineWidth',1);
    end

    % Callback for mouse down on controls points
    function oncontrolsdown(this,src,ev)
      this.is_down = true;
      % tell window that drag and up events should be handled by controls
      set(gcf,'windowbuttonmotionfcn',@this.oncontrolsdrag);
      set(gcf,'windowbuttonupfcn',    @this.oncontrolsup  );
      set(gcf,'KeyPressFcn',          @this.onkeypress    );
      % gather type of click
      if(strcmp('normal',get(gcf,'SelectionType')))
        % left-click
        this.down_type = 'left';
      else
        % other (right) click
        this.down_type = 'right';
      end
      % get current mouse position, and remember old one
      this.down_pos = get(gca,'currentpoint');
      % get position of click halfway between front and back
      middle_down_pos = mean(this.down_pos(:,1:size(this.V,2)));
      % keep track of vector from camera plane and front plane of axes
      cp = get(gca,'CameraPosition');
      ct = get(gca,'CameraTarget');
      % vector from camera position to camera target
      cp2ct = ct-cp;
      % length of vector from camera position to camera target
      cp2ct_len = sqrt(sum(cp2ct.^2));
      % height of triangle between [cp,ct,down_pos] treating cp2ct as base
      h = doublearea([cp;ct;this.down_pos],[1 2 3])/2/cp2ct_len;
      % vector from cp to perpendicular bisector of base
      this.cp2cf = sqrt(sum((this.down_pos(1,:)-cp).^2)-h.^2).*cp2ct./cp2ct_len;
      % down position should have same dimension as V
      this.down_pos = this.down_pos(1,1:size(this.V,2));
      % initialize last drag position to down position
      this.last_drag_pos = this.down_pos;
      % initialize drag position to down position
      this.drag_pos = this.down_pos;
      % keep track of control point positions at mouse down
      if size(this.V,2) == 2
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],:) = ...
          [get(this.ch,'XData')' get(this.ch,'YData')'];
      else
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],:) = ...
          [get(this.ch,'XData')' get(this.ch,'YData')' get(this.ch,'ZData')'];
      end
      %this.new_C(this.free,:) = ...
      %  [get(this.fch,'XData')' get(this.fch,'YData')'];
      % get index of closest control point
      [XXX,this.ci] =  ...
        min(sum((this.new_C(:,1:size(this.V,2)) - ...
          repmat(middle_down_pos,size(this.new_C,1),1)).^2,2));
      this.ci = [this.ci;setdiff(find(this.CG == this.CG(this.ci)),this.ci)];
      this = this.update_visualizations();
    end

    % Callback for mouse drag on control points
    function oncontrolsdrag(this,src,ev)
      % keep last drag position
      this.last_drag_pos = this.drag_pos;
      % get current mouse position
      this.drag_pos = get(gca,'currentpoint');
      % handle dragging outside of axes
      camera_depth = get(gca,'CameraPosition') - get(gca,'CameraTarget');
      drag_axis_depth = this.drag_pos(1,:) - this.drag_pos(2,:);
      if sum((camera_depth - drag_axis_depth).^2) < eps
        % "correct" drag pos front position
        this.drag_pos(1,:) = this.drag_pos(1,:) + this.cp2cf;
      end
      % only keep front
      this.drag_pos = this.drag_pos(1,1:size(this.V,2));

      if(strcmp('left',this.down_type))
        % move selected control point by drag offset
        this.new_C(this.ci,:) = ...
          bsxfun(@plus,this.new_C(this.ci,:),this.drag_pos-this.last_drag_pos);
      else
        % control groups not suppported when specifying rotations
        this.ci = this.ci(1);
        % store rotation at *point* handle if selected
        [found, iP] = ismember(this.ci,this.P);
        if(found)
          x_extent = max(this.V(:,1))-min(this.V(:,1));
          angle = 2*pi*(this.drag_pos(1)-this.last_drag_pos(1))/x_extent;
          switch size(this.V,2)
          case 2
            this.R(iP) = this.R(iP) + angle;
          case 3
            this.R(:,:,iP) = ...
              this.R(:,:,iP) * axisangle2matrix(drag_axis_depth,angle);
          end
        end
      end
      this = this.update_positions();
    end

    % Callback for mouse release of control points
    function oncontrolsup(this,src,ev)
      % Tell window to handle drag and up events itself
      set(gcf,'windowbuttonmotionfcn','');
      set(gcf,'windowbuttonupfcn','');
      this = this.set_axes();
      this.is_down = false;
    end

    function onkeypress(this,src,ev)
      % default is to not update on keystrok
      need_to_update = false;
      switch ev.Character
      case 'r'
        % reset to bind positions, identity rotations
        this.new_C = this.C;
        switch size(this.V,2)
        case 2 
          this.R = zeros(this.np,1);
        case 3
          this.R = repmat(eye(3,3),[1 1 this.np]);
        end
        need_to_update = true;
      case 'u'
        % force update
        need_to_update = true;
      end

      % only update positions if needed
      if need_to_update
        % simulate mouse down, up
        old_is_down = this.is_down;
        this.is_down = true;
        this = this.update_positions();
        this.is_down = old_is_down;
      end
    end

    function this = update_visualizations(this)
      % Update visualization to show weights if selected *point* handle
      %

      if ~isempty(this.ci)
        % try to find ci in list of point handles
        [found, iP] = ismember(this.ci(1),this.P);
        if(found)
          % change weights in contour plot
          if(this.show_contour)
            set(contour_handle,'ZData',this.Wr(:,:,iP));
            set(contour_handle,'LevelStep', ...
              (max(this.CW(:,iP))-min(this.CW(:,iP)))/contour_levels);
          end
          % change weights in weight visualization
          if(this.visualize_weights)
            set(this.wvsh,'CData',this.WVW(:,iP));
            V = get(this.wvsh,'Vertices');
            V(:,3) = max(max(V(:,1:2))-min(V(:,1:2))).*this.WVW(:,iP);
            set(this.wvsh,'Vertices',V);
            [az,el]= view(get(this.wvsh,'Parent'));
            axis(get(this.wvsh,'Parent'),'tight');
            view(az,el);
          end
        end
      end

      % change rotations in weight visualization
      if(this.visualize_rotations)
        % fit rotations to current positions
        CSM = covariance_scatter_matrix(this.V,this.F);
        dim = size(this.V,2);
        % dim by dim by n list of covariance matrices
        S = CSM*repmat(this.new_V,dim,1);
        S = permute(reshape(S,[size(CSM,1)/dim dim dim]),[2 3 1]);
        [R,~] = fit_rotations(S);
        % convert matrices to signed angles
        r = asin(squeeze(R(1,2,:)));
        %r = squeeze(SS(2,2,:)./SS(1,1,:))';
        %r = (r.^(1/4))*pi-pi/2;
        %%r = (r-min(r))/(max(r)-min(r));
        set(this.rvsh,'CData',r);
      end
    end

  end%private methods

  methods(Access=protected)

    function this = skinning_update(this)
      % update mesh positions using skinning, using appropriate method based on
      % interp_mode field
      %
      if(strcmp(this.interp_mode,'LBS'))
        this = this.lbs_update();
      elseif(strcmp(this.interp_mode,'DQLBS'))
        % USING DUAL QUATERNION SKINNING
        % number of handles
        m = numel(this.P)+size(this.BE,1);
        % dimension (2 or 3)
        dim = size(this.C,2);
        
        % update transformations based on given automatic degrees of freedom method
        old_R = this.R;
        if ~strcmp(this.auto_dof,'none')
          this = this.auto_dof_update();
          X = repmat([1 0],m,1);
          RX = stacktimes(this.TR(1:dim,1:dim,:),X);
          X = [X zeros(m,1)];
          RX = [RX zeros(m,1)];
          [AX,AN] = axisanglebetween(X,RX,[0 0 1]);
          switch size(this.V,2)
          case 2
            this.R = this.R + AX(1:this.np,3).*AN(1:this.np);
          case 3
            error('Not supported');
          end
        end

        % convert angles around z-axis to quaternions
        [T,AX,AN,Sm,O] = ...
          skinning_transformations(this.C,this.P,this.BE,this.new_C,this.R);
        % restore rotation
        this.R = old_R;

        % Extract scale
        this.TR = zeros(dim,dim+1,m);
        this.TR(1:dim,1:dim,:) = Sm;
        Sm = reshape(Sm,[dim dim*m])';
        this.TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
        % Perform scale as linear blend skinning, before translations and
        % rotations
        [this.new_V] = lbs(this.V(:,1:dim),this.TR,this.W);
        % convert rotations to axis angle
        Q = axisangle2quat(AX,AN);
        % quattrans2udq expect 3D translations, so pad with zeros
        T = [T zeros(size(T,1),1)];
        % convert quaternions and translations into dualquaternions
        DQ = quattrans2udq(Q,T);
        % compute deformation using dual quaternion skinning
        this.new_V = dualquatlbs(this.new_V,DQ,this.W);
      end

    end

    function this = lbs_update(this)
      % update mesh positions using skinning
      %

      % get transformations stored at each point and bone handle
      this.TR = ...
        skinning_transformations(this.C,this.P,this.BE,this.new_C,this.R);
      % update transformations based on given automatic degrees of freedom method
      this = this.auto_dof_update();

      dim = size(this.C,2);
      if(this.stretch_bones)
        % stretchable twistable bones skinning
        this.new_V = zeros(size(this.V));
        % point handles
        this.new_V = ...
          this.new_V + ...
          lbs(this.V(:,1:dim),this.TR(:,:,1:this.np),this.W(:,1:this.np));
        % stretchable bones
        this.new_V = ...
          this.new_V + ...
          stretch_bones_lbs( ...
            this.V(:,1:dim), ...
            this.C, ...
            this.BE, ...
            this.new_C, ...
            this.W(:, ...
            this.np+(1:this.nb)), ...
            this.WBE);
      elseif(this.stiff_points)
        % stiff points lbs
        this.new_V = ...
          stiff_points_lbs( ...
            this.V(:,1:dim), ...
            this.C, ...
            this.P, ...
            this.BE, ...
            this.new_C, ...
            this.W, ...
            this.WBE);
      else
        % regular linear blend skinning
        [this.new_V] = lbs(this.V(:,1:dim),this.TR,this.W);
        if(this.show_pixel_diff)
          [this.other_V] = lbs(this.V(:,1:dim),this.TR,this.Wother);
        end
      end
    end
      
    function [this] = auto_dof_update(this)
      % compute degrees of freedom automatically based on control point
      % positions

      % handle case where all handles are fixed
      indices = 1:(this.np+this.nb);
      if(all(ismember(indices,this.fixed)))
        return;
      end

      % displacements of *point* handles
      D = this.new_C(this.P,:)-this.C(this.P,:);

      switch this.auto_dof
      case 'arap'
        % Compute as rigid as possible *point* handle degrees of freedom
        % doesn't work with bones
        [this.L,this.new_C] = ...
           arap_dof( ...
             this.V,this.F,this.W, ...
             this.C,this.P,this.BE,this.new_C, ...
             'L0',this.L, ...
             ... %'CovarianceScatterMatrix',this.Lall, ...
             'Groups',this.G, ...
             'Free',this.free, ...
             'Fixed',this.fixed,this.TR(:,:,this.fixed), ...
             'Tol', this.tol, ...
             'MaxIter', this.max_iterations);
        dim = size(this.V,2);
        % place L into TR in appropriate place
        this.TR(:,:,1:(this.np+this.nb)) = ...
          permute(reshape(this.L,[this.np+this.nb dim dim+1]),[2 3 1]);
      case 'pseudoedges'
        switch size(this.V,2)
        case 2
          Q = axisangle2quat(repmat([0 0 1],this.np,1),this.R);
        case 3
          Q = dcm2quat(this.R);
        end
        % compute rotations using pseudoedges
        this.TR(:,:,1:this.np) = ...
          pseudoedge_dof( ...
            this.C(this.P,:), ...
            this.PE,this.new_C(this.P,:)-this.C(this.P,:), ...
            Q);
      case 'none'
        % do nothing
      otherwise
        error(['''' this.autodof ''' is not a valid AutoDOF value']);
      end
    end

    function this = update_display_positions(this)
      % update display positions
      %

      % update mesh positions in main plot
      if(this.grid)
        set(this.gsh,'XData',reshape(this.new_V(:,1),this.grid_h,this.grid_w));
        set(this.gsh,'YData',reshape(this.new_V(:,2),this.grid_h,this.grid_w));
        if this.color_by_rotations
          assert(false);
        end
      else
        if size(this.V,2) == 2
          set(this.tsh,'Vertices',[this.new_V(:,1:2) this.Zdepth]);
        else
          set(this.tsh,'Vertices',this.new_V);
        end
        % update color by rotations
        if this.color_by_rotations
          % fit rotations to current positions
          R = fit_rotations(this.V,this.F,this.new_V);
          % convert matrices to signed angles
          r = asin(squeeze(R(1,2,:)));
          set(this.tsh,'CData',r);
        end
      end

      % update pixel difference plot
      if(this.show_pixel_diff)
        pixel_diff = normrow(this.other_V - this.new_V);
        set(this.wpdsh,'CData',pixel_diff);
      end

      % update control point positions
      set(this.ch,'XData', ...
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],1));
      set(this.ch,'YData', ...
        this.new_C([this.interpolated() this.free_bone_endpoint_indices()],2));
      if size(this.V,2) == 3
        set(this.ch,'ZData', ...
          this.new_C([this.interpolated() this.free_bone_endpoint_indices()],3));
      end
      %set(this.fch,'XData',this.new_C(this.free,1));
      %set(this.fch,'YData',this.new_C(this.free,2));
      % update bone plots
      if(this.nb > 0)
        set(this.B_plot_outer,{'XData'}, num2cell([ ...
          this.new_C(this.BE(:,1),1) ...
          this.new_C(this.BE(:,2),1)],2));
        set(this.B_plot_outer,{'YData'}, num2cell([ ...
          this.new_C(this.BE(:,1),2) ...
          this.new_C(this.BE(:,2),2)],2));
        set(this.B_plot_inner,{'XData'}, num2cell([ ...
          this.new_C(this.BE(:,1),1) ...
          this.new_C(this.BE(:,2),1)],2));
        set(this.B_plot_inner,{'YData'}, num2cell([ ...
          this.new_C(this.BE(:,1),2) ...
          this.new_C(this.BE(:,2),2)],2));
      end
      % update pseudo edge plots
      if(size(this.PE,1)>0)
        set(this.PE_plot,{'XData'}, num2cell([ ...
          this.new_C(this.P(this.PE(:,1)),1) ...
          this.new_C(this.P(this.PE(:,2)),1)],2));
        set(this.PE_plot,{'YData'}, num2cell([ ...
          this.new_C(this.P(this.PE(:,1)),2) ...
          this.new_C(this.P(this.PE(:,2)),2)],2));
      end
      % update cage edge plots
      if(size(this.CE,1)>0)
        set(this.CE_plot_outer,{'XData'}, num2cell([ ...
          this.new_C(this.P(this.CE(:,1)),1) ...
          this.new_C(this.P(this.CE(:,2)),1)],2));
        set(this.CE_plot_outer,{'YData'}, num2cell([ ...
          this.new_C(this.P(this.CE(:,1)),2) ...
          this.new_C(this.P(this.CE(:,2)),2)],2));
        set(this.CE_plot_inner,{'XData'}, num2cell([ ...
          this.new_C(this.P(this.CE(:,1)),1) ...
          this.new_C(this.P(this.CE(:,2)),1)],2));
        set(this.CE_plot_inner,{'YData'}, num2cell([ ...
          this.new_C(this.P(this.CE(:,1)),2) ...
          this.new_C(this.P(this.CE(:,2)),2)],2));
      end
    end


    function this = update_slaves(this)
      %
      % update slave deform objects
      if this.is_down
        % BFS to mark all reachable slaves
        slaves = this.slaves;
        slaves_Q = this.slaves;
        while numel(slaves_Q) > 0
          that = slaves_Q(end);
          slaves_Q = slaves_Q(1:(end-1));
          if (this ~= that) && that.slave_seen
            that.slave_seen = false;
            slaves = [slaves that.slaves];
            slaves_Q = [slaves_Q that.slaves];
          end
        end
        % be sure there are no repeats
        slaves = setdiff(unique(slaves),[this]);
        if ~isempty(slaves)
          for that = slaves
            that.slave_seen = true;
            %assert(this ~= that);
            if ~that.is_down
              % control assume master controls correspond to beginning of slave
              % controls
              % Only control as many handles of slave as possible
              nc = min(size(this.new_C,1),size(that.new_C,1));
              switch size(this.V,2)
              case 2
                nr = min(size(this.R,1),    size(that.R,1));
                that.R(1:nr) = this.R(1:nr);
              case 3
                nr = min(size(this.R,3),    size(that.R,3));
                that.R(:,:,1:nr) = this.R(:,:,1:nr);
              end
              that.new_C(1:nc,:) = this.new_C(1:nc,:);
              that.update_positions();
            end
          end
        end
      end
    end
  end% protected

end% classdef deform

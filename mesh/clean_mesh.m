function [SV,SF,SVJ] = clean(V,F,varargin)
  % CLEAN  Clean up a given mesh (V,F)
  %
  % [SV,SF,SVJ] = clean(V,F)
  % [SV,SF,SVJ] = clean(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh positions
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'MaxIter' followed by maximum number of iterations to try {100}
  %     'MinAngle' followed by minimum angle threshold in radians for
  %       degenerate triangles {0.001}
  %     'MinArea' followed by minimum area threshold for small triangles {1e-7}
  %     'MinDist' followed by minimum distance threshold for duplicate vertices
  %       {1e-7}
  %     'Quiet' followed by true or {false}
  %     'SelfIntersections' followed by one of the following:
  %        'mesh'  subdivide along self-intersection contour
  %        'remove'  remove self-intersecting face pairs
  %        'remove-first'  remove first of each self-intersecting faces pair
  %        'ignore'  ignores them
  %     'Single' followed by whether to cast to single precision {false}
  %     'SmallTriangles' followed by one of the following:
  %        'collapse'  collapse small triangles along shortest edge (possibly
  %          creating self-intersections and duplicate facets)
  %        'remove'  simply remove small triangles (creating holes)
  % Outputs:
  %   SV  #SV by 3 list of new mesh positions
  %   SF  #SF by 3 list of new triangle indices
  %   SVJ  #V by 1 list of indices mapping V to SV
  %
  %
  % Example:
  %  [SV,SF,SVJ] = clean(V,F);
  %  [DV,DT,DF,DN] = cdt(SV,SF);
  %  % remap final vertices and final tets to work with original faces
  %  IM = [SVJ(:)' max(SVJ)+1:size(DV,1)];
  %  RIM = 1:size(DV,1);
  %  RIM(IM) = 1:numel(IM);
  %  % Extension of original vertices
  %  VDV = DV(IM,:);
  %  % final tets and faces defined over extension of original vertices
  %  VDF = RIM(DF);
  %  VDT = RIM(DT);
  %
  %  % Remove only degeneracies and duplicates
  %  [SV,SF] = clean(V,F,'MinDist',0,'MinArea',0,'MinAngle',0, ...
  %    'SelfIntersections','ignore','SmallTriangles','remove');
  %

  % Default parameters
  % around 1e-7 Tetgen stops complaining (i.e. seg faulting)
  min_dist = 1e-7;
  min_area = 1e-7;
  min_angle = 0.001;
  use_single = false;
  % Small triangles will contribute ~0 to winding number. Perhaps it's OK to
  % just delete them. Actually this logic is flawed because clean up is only
  % used for the constraint facets, not the winding number computation
  %small_triangles = 'collapse';
  small_triangles = 'remove';
  self_intersections = 'mesh';
  max_iter = 100;
  quiet = false;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','MinAngle','MinArea','MinDist','Quiet','SelfIntersections', ...
      'Single', 'SmallTriangles'}, ...
    {'max_iter','min_angle','min_area','min_dist','quiet', ...
      'self_intersections', 'use_single','small_triangles'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end


  if quiet
    fid = fopen('/dev/null','w');
  else
    fid = 1;
  end

  % Q: Is this guaranteed to converge?
  % A: No.
  SVJbefore = 1:size(V,1);
  SVJ = 1:size(V,1);
  SV = V;
  SF = F;
  iter = 0;
  tic;
  while iter < max_iter
    SV = double(single(SV)); 

    % GEOMETRICALLY DUPLICATE VERTICES (b.c. selfintersect maintains manifolds)
    if min_dist >= 0
      [SV,SVI,SVJ] = remove_duplicate_vertices(SV,min_dist);
      num_dups = size(SVJ,1)-size(SVI,1);
      % remap faces
      SF = SVJ(SF);
      % remap old indices
      SVJ = SVJ(SVJbefore);
      fprintf(fid,'Removed %d geometrically duplicate vertices\n',num_dups);
    end
    [SF,SFI,SFJ] = remove_duplicate_simplices(SF);
    fprintf(fid,'Removed %d combinatorially duplicate facets\n', ...
      size(SFJ,1)-size(SFI,1));
    % COMBINATORIALLY DEGENERATE FACETS
    NSF = SF((SF(:,1) ~= SF(:,2))&(SF(:,2) ~= SF(:,3))&(SF(:,3) ~= SF(:,1)),:);
    if size(NSF,1) < size(SF,1)
      fprintf(fid,'Removed %d facets with combinatorally duplicate vertices\n', ...
        size(SF,1)-size(NSF,1));
      SF = NSF;
    end

    SVbefore = SV;
    SFbefore = SF;
    SVJbefore = SVJ;

    % SMALL TRIANGLES
    switch small_triangles
    case 'collapse'
      SF = collapse_small_triangles(SVbefore,SFbefore,min_area);
    case 'remove'
      SF = SFbefore(doublearea(SVbefore,SFbefore)>min_area,:);
    otherwise
      error(['Unsupported SmallTriangles argument: ' small_triangles]);
    end
    fprintf(fid,'%sd %d small triangles\n', ...
      small_triangles,size(SFbefore,1)-size(SF,1));

    % SMALL ANGLE TRIANGLES
    sm_angle = min(internalangles(SVbefore,SF),[],2)<min_angle;
    SF = SF(~sm_angle,:);
    fprintf(fid,'Removed %d small angle triangles\n',sum(sm_angle));

    % COMBINATORIALLY DUPLICATE FACETS
    [~,SFI,SFJ] = unique(sort(SF,2),'rows','stable');
    SF = SF(SFI,:);
    fprintf(fid,'Removed %d combinatorially duplicate facets\n', ...
      size(SFJ,1)-size(SFI,1));

    if strcmp(self_intersections,'ignore')
      break;
    end
    % SELF-INTERSECTIONS
    [SVtemp,SFtemp,IF] = selfintersect(SVbefore,SF, ...
      'DetectOnly',~strcmp(self_intersections,'mesh'));
    % number of intersection pairs
    num_inters = size(IF,1);
    switch self_intersections
    case 'mesh'
      fprintf(fid,'Added %d vertices and %d faces to mesh self-intersections\n', ...
        size(SVtemp,1) - size(SVbefore,1), ...
        size(SFtemp,1) - size(SF,1));
      SV = SVtemp;
      SF = SFtemp;
    case 'remove'
      offending = ismember(1:size(SF,1),IF(:));
      SF = SF(~offending,:);
      fprintf(fid,'Removed %d self-intersecting faces\n',sum(offending));
    case 'remove-first'
      % Q: Is there an optimal (smallest) number of facets to remove that would
      % illiminate all self-intersections?
      % A: What about sorting each pair by the number of self-intersections each
      % facet participates in. Then removing the first.
      offending = ismember(1:size(SF,1),IF(:,1));
      SF = SF(~offending,:);
      fprintf(fid,'Removed %d self-intersecting faces\n',sum(offending));
    case 'remove-optimal'
      % Q: Would this be more efficient if IF was sparse and we zeroed-out
      % entries?
      rmoff = [];
      while ~isempty(IF)
        % find number of occurences of each offending facet
        C = sparse(IF(:),1,1);
        [~,ii] = max(C);
        % remove insections involving this facet
        IF = IF(~any(IF==ii,2),:);
        rmoff = [rmoff(:);ii];
      end
      SF = SF(~ismember(1:end,rmoff),:);
      fprintf(fid,'Removed %d self-intersecting faces\n',numel(rmoff));
    otherwise
      error(['Unsupported SelfIntersections argument: ' self_intersections]);
    end
    if num_inters == 0
      break;
    end

    % TODO: when repeating we should only check for self-intersections on dirty
    % faces
    iter = iter + 1;
  end
  toc


  % No combinatorially degenerate faces
  assert(~any(SF(:,2)==SF(:,3)));
  assert(~any(SF(:,3)==SF(:,1)));
  assert(~any(SF(:,1)==SF(:,2)));

end

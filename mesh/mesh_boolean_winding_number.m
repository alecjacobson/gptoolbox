function [W,H,J] = mesh_boolean_winding_number(V,F,U,G,operation,varargin)
  % MESH_BOOLEAN_WINDING_NUMBER Compute boolean csg operations on oriented
  % meshes. Uses winding number to determine boolean extraction. This is
  % theoretically prone to errors for "bad" intersections, but will handle
  % "dirty" input meshes (e.g. open boundaries) to some degree.
  %
  % [W,H] = mesh_boolean(V,F,U,G,operation)
  % 
  % Inputs:
  %   V  #V by 3 list of vertex positions of first mesh
  %   F  #F by 3 list of triangle indices into V
  %   U  #U by 3 list of vertex positions of second mesh
  %   G  #G by 3 list of triangle indices into U
  %   operation  followed by operation to perform as a string, one of: 'union',
  %     'intersect', 'minus', 'xor', or 'resolve'
  %     Optional:
  %       'Thresholds' 2-long list of winding number cut offs to use to
  %         determine inside/outside: all faces of F (after intersection) with
  %         winding number w.r.t. (U,G) less than thresholds(1) are considered
  %         inside.
  % Outputs:
  %   W  #W by 3 list of vertex positions of boolean result mesh
  %   H  #H by 3 list of triangle indices into W
  %   J  #H list of indices into [FA;FB] of facet birth parents
  %    
  % See also: mesh_boolean
  %

  VU = [V;U];
  FG = [F;size(V,1)+G];
  [W,SFG,~,J,IM] = selfintersect(VU,FG);
  [W,UIM] = remove_unreferenced(W,IM(SFG));
  SF = IM(SFG(J<=size(F,1),:));
  JF = J(J<=size(F,1));
  SG = IM(SFG(J>size(F,1),:));
  JG = J(J>size(F,1),:);
  SF = UIM(SF);
  SG = UIM(SG);

  th = [0.5 0.5];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Thresholds'}, ...
    {'th'});
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
  A = barycenter(W,SF);
  B = barycenter(W,SG);
  wUGA = winding_number(U,G,A);
  wVFB = winding_number(V,F,B);
  switch operation
  case 'union'
    SF = SF(wUGA<th(1),:);
    SG = SG(wVFB<th(2),:);
    JF = JF(wUGA<th(1),:);
    JG = JG(wVFB<th(2),:);
  case 'intersect'
    SF = SF(wUGA>=th(1),:);
    SG = SG(wVFB>=th(2),:);
    JF = JF(wUGA>=th(1),:);
    JG = JG(wVFB>=th(2),:);
  case 'minus'
    SF = SF(wUGA<th(1),:);
    SG = fliplr(SG(wVFB>=th(2),:));
    JF = JF(wUGA<th(1),:);
    JG = JG(wVFB>=th(2),:);
  case 'xor'
    SF = [SF(wUGA<th(1),:); fliplr(SF(wUGA>=th(1),:))];
    SG = [SG(wVFB<th(2),:); fliplr(SG(wVFB>=th(2),:))];
    JF = [JF(wUGA<th(1),:); JF(wUGA>=th(1),:)];
    JG = [JG(wVFB<th(2),:); JG(wVFB>=th(2),:)];
  case 'resolve'
    % Don't do anyting. Just snap vertices...
  otherwise
    error(['Unsupported operation: ' operation]);
  end
  % Combine meshes and map to unique vertices
  H = [SF;SG];
  J = [JF;JG];
  [W,IM] = remove_unreferenced(W,H);
  H = IM(H);
end

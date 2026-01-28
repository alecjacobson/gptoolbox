function [V,F,I] = capsule_mesh(n,a,varargin)
  % [V,F,I] = capsule_mesh(n,a)
  %
  % Inputs:
  %   n  number of faces along latitude of each hemisphere
  %   a  length of cylindrical portion
  %   Optional:
  %     'R' followed by scalar radius of capsule {1}
  %     'Stacks' followed by number of stacks in cylindrical part
  %       {ceil(2*n/(2*pi*R*a))}
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into rows of V
  %   I  #F list I(f) = 1,2 or 3 indicating bottom hemisphere, cylindrical part,
  %     top hemisphere, respectively.
  %
  % See also: stadium

  R = 1;
  s = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'R','Stacks'},{'R','s'});
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

  if isempty(s)
    s = ceil( (2*n)/(2*pi*R)* a);
  end
  [CV,CF] = cylinder_mesh(R,2*n,'Stacks',s);
  CV = CV.*[1 1 a];
  [X,Y,Z] = sphere(2*n);
  [SF,SV] = surf2patch(X,Y,Z,'triangles');
  SV = SV.*R;
  % Eww
  [SV,~,J] = remove_duplicate_vertices(SV,1e-7);
  SF=J(SF);
  SF = SF(SF(:,1)~=SF(:,2) & SF(:,2)~=SF(:,3) & SF(:,3)~=SF(:,1),:);
  SBC = barycenter(SV,SF);
  T = SBC(:,3)>0;
  TV = SV+[0 0 a];
  TF = SF(T,:);
  BV = SV;
  BF = SF(~T,:);
  V = [BV;CV;TV];
  F = [BF;size(BV,1)+[CF;size(CV,1)+TF]];
  % Eww
  [V,~,J] = remove_duplicate_vertices(V,1e-7);
  F=J(F);
  [V,~,~,F] = remove_unreferenced(V,F);
  I = [repmat(1,size(BF,1),1);repmat(2,size(CF,1),1);repmat(3,size(TF,1),1)];
end

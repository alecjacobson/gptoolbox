function [d1,d2,E] = boundary_matrix(varargin)
  % [d1,d2] = boundary_matrix(F)
  %   or
  % [d1,d2] = boundary_matrix(PI,PC)
  %
  % https://jeremykun.com/2013/04/10/computing-homology/
  % 
  % Inputs:
  %   F  #F by 3 list of triangles into some vertex set 1:max(F(:))
  %    or
  %   PI  #PI stream of polygon indices into rows of some V
  %   PC  #polys+1 list of cumulative sum of dual face valences
  % Outputs:
  %   d1  #V by #E sparse, oriented boundary matrix
  %   d2  #E by #F sparse, oriented boundary matrix
  % 
  % see also: cells, polygons_to_triangles
  %
  switch numel(varargin)
  case 1
    F = varargin{1};
    n0 = max(F(:));
    n2 = size(F,1);
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    K = repmat(1:size(F,1),3,1)';
  case 2
    PI = varargin{1};
    PC = varargin{2};
    n0 = max(PI(:));
    n2 = size(PC,1)-1;
    diffPC = diff(PC);
    I1 = PI;
    % order index per polygon
    PJ = [(1:numel(PI))' - reshape(repelem(PC(1:end-1),diffPC),[],1)];
    % increment by one mod size
    PJ = mod(PJ,reshape(repelem(diffPC,diffPC),[],1))+1;
    I2 = PI(reshape(repelem(PC(1:end-1),diffPC),[],1) + PJ);
    allE = [I1 I2];
    % all halfedges. Halfedge emanating from corresponding vertex in PI
    K = repelem((1:n2)',(PC(2:end)-PC(1:end-1))');
  end

  [E,~,EMAP] = unique(sort(allE,2),'rows');
  n1 = size(E,1);
  d1 = sparse(E,repmat(1:n1,2,1)',repmat([-1 1],n1,1),n0,n1);
  d2 = sparse(EMAP,K,2*(E(EMAP,1)==allE(:,1))-1,n1,n2);
end

function [V,F] = wedding_cake(P,E,h)
  % WEDDING_CAKE Generate a "wedding cake" mesh from a curve and a list of
  % heights.
  % 
  % [V,F] = wedding_cake(P,E,h)
  %
  % Inputs:
  %   P  #h list of #P{i} by 2 list of positions
  %   E  #h list of #E{i} by 3 list of edge indices into P{i}
  %   h  #h list of level heights
  % Outputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %
  
  assert(numel(P) == numel(E));
  assert(numel(P) == numel(h));
  % Loop over curves
  % running mesh
  V = [];
  F = [];
  tot_h = 0;
  for v = 1:numel(P)+1
    if v<=numel(P)
      vP = P{v};
      vE = E{v};
    else
      vP = [];
      vE = [];
    end
    if v ~= 1
      prevP = P{v-1};
      prevE = E{v-1};
      prevI = size(V,1)-2*size(prevP,1)+(1:size(prevP,1));
    else
      prevP = [];
      prevE = [];
      prevI = [];
    end

    %F = [F;size(V,1)-size(prevP,1)+BF];
    % If on any layer but last, then extrude upward
    if v<=numel(P)
      % wall
      [vV,vF] = extrude(vP,vE,'Cap',false);
      vV = bsxfun(@times,[1 1 h(v)],vV);
      vV = bsxfun(@plus,[0 0 tot_h],vV);
      tot_h = tot_h + h(v);
      %  tsurf(vF,vV)
      %  statistics(vV,vF,'Fast',true)
      %  input(num2str(v));
      % append
      F = [F;size(V,1)+vF];
      vI = size(V,1)+size(vP,1)+(1:size(vP,1));
      V = [V;vV];
    end

    % base or level with previous layer
    BP = [prevP;vP];
    I = [prevI vI];
    BE = [prevE;size(prevP,1)+vE];
    [BV,BF] = triangle(BP,BE,[],'Flags','-YY','Quiet');
    % Assumes we've already dealt with self-intersections in curves
    assert(size(BV,1) == size(BP,1));
    % absolute value is a little unsafe here
    w = winding_number(BP,BE,barycenter(BV,BF));
    % Once inside
    BF = BF(w>0.5 & w<1.5,:);
    % If BV == BP, then this is unnecessary
    %[BV,I] = remove_unreferenced(BV,BF);
    %BF = I(BF);
    % flip bottom layer
    if v == 1
      BF = fliplr(BF);
    end
    % remap and append
    F = [F;I(BF)];
    assert(~any(isnan(V(:))))
  end

  % cleanup
  [V,~,I] = remove_duplicate_vertices(V,1e-9);
  F = I(F);

  %V = [AV;BV];
  %F = [AF;size(AV,1)+BF];

  %  [V,F] = triangle(V,O,[],'Flags',sprintf('-q -a%g',max_area));
  %  w = winding_number(S,E,barycenter(V,F))/(2*pi);
  %  F = F(w>0.5,:);
  %  [V,I] = remove_unreferenced(V,F);
  %  F = I(F);

  %end
  
end

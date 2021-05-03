function [PC,CC,DC,IC] = trim_with_spline(PA,CA,PB,CB,tol)
  % TRIM_WITH_SPLINE Trim a given spline (PA,CA) with a "solid" spline (PB,CB)
  %
  % [PC,CC,DC,IC] = trim_with_spline(PA,CA,PB,CB,tol)
  %
  % Inputs:
  %   PA  #PA by dim list of control point locations
  %   CA  #CA by 4 list of indices into PA of cubic Bezier curves
  %   PB  #PB by dim list of control point locations
  %   CB  #CB by 4 list of indices into PB of cubic Bezier curves
  %   tol  tolerance for intersection {1e-7}
  % Outputs:
  %   PC  #PC by dim list of control point locations
  %   CC  #CC by 4 list of indices into PC of cubic Bezier curves
  %   DC  #CC list of flags whether inside B
  %   IC  #CC list of indices into CA
  %

  function [PA,CA,DA,IA] = trim_with_spline_helper(PA,CA,PB,CB)

    [A1,A2] = box_each_element(PA,CA);
    [B1,B2] = box_each_element(PB,CB);
    I = box_intersect(A1,A2,B1,B2);
    T = cell(size(CA,1),1);
    % consider each intersection, gather splits
    for ii = 1:size(I,1)
      ca = I(ii,1);
      cb = I(ii,2);
      Tii = cubic_cubic_intersect(PA(CA(ca,:),:),PB(CB(cb,:),:),tol);
      if ~isempty(Tii)
        T{ca} = [T{ca};Tii(:,1)];
      end
    end
    % conduct all splits
    IA = (1:size(CA,1))';
    for ca = 1:size(CA,1)
      if ~isempty(T{ca})
        Tca = sort(T{ca});
        [PAa,CAa] = cubic_subdivide(PA(CA(ca,:),:),Tca);
        CAa = reshape([CA(ca,1) CAa(2:end-1)+size(PA,1)-2 CA(ca,end)],size(CAa));
        IA([ca size(CA,1)+(1:size(CAa,1)-1)],:) = IA(ca);
        CA([ca size(CA,1)+(1:size(CAa,1)-1)],:) = CAa;
        PA = [PA;PAa(3:end,:)];
      end
    end

    %%fprintf('trim_with_spline...\n');
    %% loop order is imporant 
    %IA = (1:size(CA,1))';
    %for cb = 1:size(CB,1)
    %  %fprintf('  %04d/%04d:\n',cb,size(CB,1));
    %  for ca = 1:size(CA,1)
    %    %progressbar(ca,size(CA,1));
    %    T = cubic_cubic_intersect(PA(CA(ca,:),:),PB(CB(cb,:),:),tol);
    %    if ~isempty(T)
    %      [PAa,CAa] = cubic_subdivide(PA(CA(ca,:),:),T(:,1));
    %      CAa = reshape([CA(ca,1) CAa(2:end-1)+size(PA,1)-2 CA(ca,end)],size(CAa));
    %      IA([ca size(CA,1)+(1:size(CAa,1)-1)],:) = IA(ca);
    %      CA([ca size(CA,1)+(1:size(CAa,1)-1)],:) = CAa;
    %      PA = [PA;PAa(3:end,:)];
    %    end
    %  end
    %end

    % midpoints of each cubic
    %fprintf('midpoints...\n');
    M = zeros(size(CA,1),2);
    for ca = 1:size(CA,1)
      M(ca,:) = cubic_eval(PA(CA(ca,:),:),0.5);
    end
    %fprintf('spline_winding_number...\n');
    DA = abs(spline_winding_number(PB,CB,M))>0.5;
  end

  m = size(CA,1);
  embed = @(P,C) reshape(P(C,:),size(C,1),8);
  CB8 = embed(PB,CB);
  % In either direction
  RB8 = embed(PB,fliplr(CB));
  CA8 = embed(PA,CA);
  % find perfect matches (only considering that every control point is exactly
  % the same; not finding partial co-incidence )
  [idx,dist] = rangesearch([CB8;RB8],CA8,0);
  perfect = cellfun(@(i) ~isempty(i),idx);
  JN = find(~perfect);
  [PA,CN,DN,IN] = trim_with_spline_helper(PA,CA(JN,:),PB,CB);
  CA = [CN;CA(perfect,:)];
  JP = find(perfect);
  IA = [JN(IN);JP];
  DA = [DN;true(numel(JP),1)];
  assert(size(CA,1) == numel(IA));
  assert(max(IA) <= m);
  assert(size(CA,1) == numel(DA));

  % rename
  PC = PA;
  CC = CA;
  DC = DA;
  IC = IA;
end

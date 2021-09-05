function [EF,EI,uE,EMAP] = edge_flaps(F)
  % [EF,EI,uE,EMAP] = edge_flaps(F)
  %
  % Inputs:
  %   F  #F by 3 list of face indices
  % Outputs:
  %   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
  %     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
  %     e=(j->i)
  %   EI  #E by 2 list of edge flap corners (see above).
  %   uE  #uE by 2 list of edge indices into V.
  %   EMAP #F*3 list of indices into uE, mapping each directed edge to unique
  %     unique edge in uE
  %
  E = [F(:,2:3);F(:,[3 1]);,F(:,1:2)];
  sE = sort(E,2);
  [uE,~,EMAP] = unique(sE,'rows');
  I = (1:size(uE,1))';
  [B1,E1] = ismember(uE,E,'rows');
  [B2,E2] = ismember(uE,fliplr(E),'rows');
  J1 = mod(E1(B1)-1,size(F,1))+1;
  J2 = mod(E2(B2)-1,size(F,1))+1;
  C1 = floor((E1(B1)-1)/size(F,1))+1;
  C2 = floor((E2(B2)-1)/size(F,1))+1;
  I1 = I(B1);
  I2 = I(B2);
  K1 = repmat(1,numel(I1),1);
  K2 = repmat(2,numel(I2),1);
  EF = full(sparse([I1;I2],[K1;K2],[J1;J2],size(uE,1),2));
  EI = full(sparse([I1;I2],[K1;K2],[C1;C2],size(uE,1),2));
end

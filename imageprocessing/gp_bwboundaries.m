function [B,L,N,P] = gp_bwboundaries(BW,CONN)
  % GP_BWBOUNDARIES Find boundaries of regions and holes in a logical image
  %
  % B = gp_bwboundaries(BW)
  % [B,L,N] = gp_bwboundaries(BW,CONN)
  %
  % Inputs:
  %   BW  h by w logical BWage
  % Outputs:
  %   B  cell array containing where each element B{b} contains a #B{c} by 2
  %     list of row and column indices of pixels in the bth boundary
  %   L  h by w BWage of ids of regions and holes (0 for "unlabeled" false
  %     pixels touching the border)
  %   N  number of regions (first N elements of B are regions, the rest are
  %     "holes")
  %
  % See also: gp_bwlabel
  %

  if nargin<2
    CONN = 8;
  end

  A4 = fd_laplacian(size(BW))~=0;
  A8 = ((A4*1*A4)>1) + A4;
  A = {[],[]};
  switch CONN
  case 4
    A{1} = A8;
    A{2} = A4;
  case 8
    A{1} = A4;
    A{2} = A8;
  end
  
  LN = gp_bwlabel(BW);
  LP = gp_bwlabel(~BW);
  % Combine
  L = LN + (LN==0).*(LP+max(LN(:)));
  
  
  % BFS from outside
  O = reshape(full(sparse([2:size(BW,1)-1 1:size(BW,1):numel(BW)-size(BW,1) size(BW,1):size(BW,1):numel(BW) numel(BW)-size(BW,1)+1:numel(BW)-1],1,1)),size(BW));
  
  K = {O&BW,false(size(BW))};
  outside = [];
  inside = [];
  
  outside = union(outside,setdiff(unique(L.*O.*~BW),0));
  
  while true
    G = ismember(L,union(inside,outside));
    %BWshow(G);
    %pause
    K{1} = K{1} | (reshape((A{1}*1*G(:))~=0,size(BW)) & BW & ~G);
    %BWshow(matrixnormalize(K{1}-K{2}))
    %pause
    inside = union(inside,setdiff(unique(L.*(K{1}|K{2})),0));
    if isempty(setdiff(L,union(inside,outside)))
      break;
    end
    G = ismember(L,union(inside,outside));
    K{2} = K{2} | (reshape((A{2}*1*G(:))~=0,size(BW)) & ~G);
    outside = setdiff(union(outside,L.*((~BW & reshape((A{1}*1*G(:))~=0,size(BW))))),0);
    %BWshow(matrixnormalize(K{1}-K{2}))
    %pause
    if isempty(setdiff(L,union(inside,outside)))
      break;
    end
  end
  
  B = {{},{}};
  for c = 1:max(L(:))
    for e = 1:2
      BWc = K{e} & L==c;
      [Lc,ncc] = gp_bwlabel(BWc,8);
      for cc = 1:ncc
        BWcc = Lc == cc;
        % Super inefficient search
        Acc = A{3-e};
        Acc(~BWcc,:) = 0;
        Acc(:,~BWcc) = 0;
        Acc = Acc - diag(diag(Acc));
        [s,t] = find(Acc,1,'first');
        P = [];
        p = s;
        while sum(Acc(:))>0
          P = [P;p];
          Acc(p,:) = 0;
          [~,p] = max(Acc(:,p));
        end
        [LI,LJ] = ind2sub(size(BW),P);
        if ~isempty(LI)
          % 
          B{e} = {B{e}{:} [LI([1:end 1]) LJ([1:end 1])]}';
        end
      end
    end
  end
  N = numel(B{1});
  B = {B{1}{:} B{2}{:}}';
  
  L = L.*~ismember(L,L.*(O&~BW));
  [~,~,L] = unique(L);
  L = reshape(L,size(BW));
  
end

function [B,L,N] = gp_bwboundaries(BW,CONN)
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

  A4 = fd_laplacian(size(BW))~=0;
  A8 = (A4*1*A4)>1;
  switch CONN
  case 4
    AP = A8;
    AN = A4;
  case 8
    AP = A4;
    AN = A8;
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
    K{1} = K{1} | (reshape((AP*1*G(:))~=0,size(BW)) & BW & ~G);
    %BWshow(matrixnormalize(K{1}-K{2}))
    %pause
    inside = union(inside,setdiff(unique(L.*(K{1}|K{2})),0));
    if isempty(setdiff(L,union(inside,outside)))
      break;
    end
    G = ismember(L,union(inside,outside));
    K{2} = K{2} | (reshape((AN*1*G(:))~=0,size(BW)) & ~G);
    outside = setdiff(union(outside,L.*((~BW & reshape((AP*1*G(:))~=0,size(BW))))),0);
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
        [LI,LJ] = find(Lc == cc);
        if ~isempty(LI)
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

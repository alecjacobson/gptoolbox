function [ears,ear_opp,flops,flop_opp] = find_ears(F)
  % FIND_EARS  Find all ears (faces with two boundary edges) in a given mesh
  % 
  % [ears,ear_opp,flops,flop_opp] = find_ears(F)
  %
  % Inputs:
  %   F  #F by 3 list of triangle mesh indices
  % Outputs:
  %   ears  #ears list of indices into F of ears
  %   ear_opp  #ears list of indices indicating which edge is non-boundary
  %     (connecting to flops)
  %   flop  #ears list of indices into F of faces incident on ears
  %   flop_opp  #ears list of indices indicating which edge is non-boundary
  %     (connecting to ears)
  % 
  % Known bugs: (V,F) mmust be manifold near ears.

  % Must be manifold (actually only ears need to be)
  [Fp, Fi] = tt(F);
  % Definition: an ear has 2 boundary **edges**
  ears = find(sum(Fp==-1,2)==2);
  [flops,ear_opp] = max(Fp(ears,:),[],2);
  flop_opp = Fi(sub2ind(size(Fp),ears,ear_opp));
end


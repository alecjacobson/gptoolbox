function T = triinterp(V,F,S,U)
  % TRIINTERP
  %
  % T = triinterp(V,F,S,U)
  %
  % Given a scalar field S defined over a triangle mesh with vertices V and
  % faces F, linearly interpolate the data for new positions U. NaNs are
  % returned if position in U is not in or on a mesh triangle
  %
  % Inputs:
  %   V vertex position #V x 3 or #V x 2
  %   F list of faces #F x 3
  %   S scalar field defined over V, #V x 1
  %   U vertex positions #U x 3 or #V x 2
  % Outputs:
  %   T scalar field defined over U, #U x 1
  %

  assert(size(S,2) == 1);

  if(size(V,2) == 2)
    V = [V zeros(size(V,1),1)];
  end
  if(size(U,2) == 2)
    U = [U zeros(size(U,1),1)];
  end

  % gather triangle corners
  VF1 = V(F(:,1),:);
  VF2 = V(F(:,2),:);
  VF3 = V(F(:,3),:);

  % gather triangle areas
  dblA = sqrt(sum(cross(VF2-VF1,VF3-VF1,2).^2,2));
  % gather unnormalized barycenteric coordinates, dblAk(i,j) is the barycentric
  % coordinate of corner k of face i of jth vertex in U

  VF1 = permute(repmat(VF1,[1,1,size(U,1)]),[1,3,2]);
  VF2 = permute(repmat(VF2,[1,1,size(U,1)]),[1,3,2]);
  VF3 = permute(repmat(VF3,[1,1,size(U,1)]),[1,3,2]);
  UU = permute(repmat(U,[1,1,size(F,1)]),[3,1,2]);
  dblA1 = sqrt(sum(cross(VF3-UU,VF2-UU,3).^2,3));
  dblA2 = sqrt(sum(cross(VF1-UU,VF3-UU,3).^2,3));
  dblA3 = sqrt(sum(cross(VF2-UU,VF1-UU,3).^2,3));

  dblA = repmat(dblA,[1,size(U,1)]);
  in = abs(dblA1 + dblA2 + dblA3 - dblA)<1e-10;

  % normalize and mask
  w1 = in.*(dblA1./dblA);
  w2 = in.*(dblA2./dblA);
  w3 = in.*(dblA3./dblA);

  % gather S value at triangle corners
  SF1 = S(F(:,1));
  SF2 = S(F(:,2));
  SF3 = S(F(:,3));
  SF1 = repmat(SF1,1,size(U,1));
  SF2 = repmat(SF2,1,size(U,1));
  SF3 = repmat(SF3,1,size(U,1));

  T = sum(w1.*SF1+w2.*SF2+w3.*SF3,1);
  
  T(~logical(sum(in,1))) = NaN;
end

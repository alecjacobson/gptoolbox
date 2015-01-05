function E = arap_energy(V,F,U,R,L)
  % ARAP_ENERGY  Evaluate the arap energy as described in "As-rigid-as-possible
  % Surface Modeling" by Sorkine and Alexa
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  %   U  #V by dim list of pose domain positions
  %   Optional:
  %     R  dim by dim by #F list of rotations
  %     L  #V by #V laplacian
  % Output:
  %   E  scalar energy value
  %

  % Should take optional 'Energy' parameter like arap.m
  assert(false,'Use arap_gradient instead')

  % fit rotations if not given
  if ~exist('R','var')
    assert(false)
    R = fit_rotations(V,F,U);
  end

  if ~exist('L','var')
    if(size(F,2) == 3)
      L = cotmatrix(V,F);
    elseif(size(F,2) == 4)
      L = cotmatrix3(V,F);
    else
      error('Invalid face list');
    end
  end

  E = 0;
  EF = edges(F);
  for e = 1:size(EF,1)
    i = EF(e,1);
    j = EF(e,2);
    % add (i,j)s contribution to energy at i
    E = E + L(i,j)*...
      sum(((U(i,:) - U(j,:)) - (V(i,:)-V(j,:))*(R(:,:,i)')).^2,2);
    % add (i,j)s contribution to energy at j
    E = E + L(j,i)*...
      sum(((U(j,:) - U(i,:)) - (V(j,:)-V(i,:))*(R(:,:,j)')).^2,2);
  end
end

function [U,phi,vals,K,M] = stvk_modal_derivatives(V,T,nummodes,varargin)
  reg_method = 'null-orthogonal';
  Aeq = [];
  nu = 0.45;
  young = 100;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Regularization','Aeq', 'young', 'nu'}, ...
    {'reg_method','Aeq', 'young', 'nu'});
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
  nu = nu.*ones(size(T, 1), 1);
  young = young .* ones(size(T, 1), 1);
  lambda = young.*nu./((1+nu).*(1-2.*nu));
  mu = .5.*young./(1+nu);
  Aeq(:,[1:3:end 2:3:end 3:3:end]) = Aeq;
  [U,phi,vals,K,M] = STVKPolys3D(V,T,nummodes, reg_method, Aeq, mu, lambda);
  K = K([1:3:end 2:3:end 3:3:end],[1:3:end 2:3:end 3:3:end]);
  M = M([1:3:end 2:3:end 3:3:end],[1:3:end 2:3:end 3:3:end]);
  U = permute(reshape(U,[size(V,2) size(V,1) size(U,2)]),[2 1 3]);
  phi = permute(reshape([phi{cellfun(@length,phi)>0}],size(V,2),size(V,1),[]),[2 1 3]);

  % Modified code from Paul Kry
  function [U, phi, vals, K, M] = STVKPolys3D(V, T, nummodes, reg_method, Aeq, mu, lambda)
  % Compute basis and polynomial coefficients
  %
  %   Issues and observations:
  %
  %   Does a lot of things better than the previous attempt, and exploits
  %   sparsity throughout (i.e., might not be too horrible for larger meshes).
  %
  %   Modal derivatives are supposed to be mass normalized (how?) but are
  %   currently just normalized.  See text after equation 15 of Barbic 2005.
  %
  %   Currently missing weighted PCA to compute the basis (see equation 16 of
  %   Barbic 2005).
  %
  %   Does not do the unpinned modal derivative solve properly, but the
  %   results seem somewhat convincing just the same.  See longer discussion
  %   in appendix D of Barbic 2005.
  %
  %   Stiffness symmetry could be exploited in making a more efficient
  %   coefficient and poly computation code generation.  
  %
  %   Could likewise make polynomial coefficients a function of Lame
  %   parameters as shown in Appendix B of Barbic 2005.  This may benefit
  %   from taking an approach to coefficient computaiton that more closely
  %   follows the appendix.
  %
  %   computeBandArea is very slow for larger models.

  
      fileNameRoot = "beam";
      fileName = "../../models/3D/" + fileNameRoot + ".mesh";
      
      pinned = false;
      computeModalDerivatives = true;
      ShowAnim = 0;
      MakeVideo = 0;  % ShowAnim must be true if you want to make a video
      vidFileName = "beamModes3D-A.mp4";
      GeneratePolys = 0;
      collectStiffnessPoly = 0;
      
      % These calls use gptoolbox
      %[V,T,F] = readMESH( fileName );
      [F,J,K] = boundary_faces(T);
      
      % Use a different loader for .msh files, and again another for tetgen
      % created files.
      
      p = V;
      t = T;
      M = kroneye(massmatrix(p,t,'barycentric'),size(V,2));
      %isM = diag(sparse(1./sqrt(diag(M))));
      faces = F;
      
      %[B, vol] = computeBandArea( p, t );
      % Why negative?
      vol = -volume(p,t);
      B = grad(V,T);
      B([1:3:end 2:3:end 3:3:end],:)=B;
      B = kroneye(B,3);
  
      
      rho = 1;    % density
%       nu = 0.35;  % Poisson ratio
%       E = 2e4;    % Young's modulus
%       mu = E / 2 / (1 + nu);                      % Lamé parameter
%       lambda = E * nu / (1 + nu) / (1 - 2 * nu);  % Lamé parameter
%       
      if ~exist('stvk_d3','file')
        F = sym( 'F', [3, 3] );     % symbolic deformation gradient
        s_lambda = sym( 's_lambda' );     % symbolic deformation gradient
        s_mu = sym( 's_mu' );     % symbolic deformation gradient
        assume( F, 'real' );
        E = 1/2 * (F'*F - eye(3));  % symbolic strain
        psi = 1/2 * s_lambda * (trace(E))^2 + s_mu * trace( E*E ); % STVK energy density
        % Compute first second and third derivative wrt deformaiton gradient
        dpsidF = sym( zeros( 9, 1 ) );
        for i = 1:numel(F)
            dpsidF(i) = diff( psi, F(i) );
        end
        d2psidF2 = sym( zeros( 9, 9 ) );
        for i = 1:numel(F)
            d2psidF2(:,i) = diff( dpsidF, F(i) );
        end   
        d3psidF3 = sym( zeros( 9, 9, 9 ) );
        for i = 1:numel(F)
            d3psidF3(:,:,i) = diff( d2psidF2, F(i) );
        end
        % Build stiffness matrix
        matlabFunction(dpsidF,'File','stvk_d.m','Vars',{F,s_lambda,s_mu});
        matlabFunction(d2psidF2,'File','stvk_d2.m','Vars',{F,s_lambda,s_mu});
        matlabFunction(d3psidF3,'File','stvk_d3.m','Vars',{F,s_lambda,s_mu});
      end
  
      C0 = stvk_d2(eye(3),lambda(1),mu(1));
      n = size(t,1);
      bigC = sparse(9*n,9*n);
      Cstack = zeros(9*n, 9);
      for el = 1:n
          Cstack((el-1)*9+1:el*9, :) = stvk_d2(eye(3),lambda(el),mu(el));
%          bigC(el*9-8:el*9,el*9-8:el*9) = - vol(el) * C0;
      end
      i = repmat((1:n*9)', [1, 9]);
      j = repmat(repmat(1:9, [9, 1]), [n, 1]) + repelem((0:n-1), 9)'*9;
      CC = sparse(i, j, Cstack, 9*n, 9*n);
      bigC = diag(sparse(-repelem(vol, 9))) * CC;
%       bigC2 = diag(sparse(-repelem(vol, 9))) * repdiag(sparse(C0), numel(vol));
%        bigC2 = kron(diag(sparse(-vol)),C0);
      K = B'*bigC*B;
      K = 0.5*(K+K');
       
      % Build lumped mass matrix
      nv = size(p,1);
      %masses = zeros( nv, 1 );
      %for el = 1:n
      %   masses( t(el,:) ) = masses( t(el,:) ) + 1 / 4 * rho * vol(el);
      %end
      
      if ( pinned )
          % work with a smaller matrix of free DOFs
          freeVerts = find( p(:,1) > (min(p(:,1))+0.01) );
          ii = reshape( [freeVerts'*3-2; freeVerts'*3-1; freeVerts'*3], [], 1 );
          modeoffset = 0; 
      else
          ii = 1:numel(p);     % Use all vert DOFs
          if isempty(Aeq)
            modeoffset = 6;  % skip the first 6 modes because of rigid modes
          else
            modeoffset = 0;
          end
      end
       
      %dofmasses = zeros( numel(p), 1 );
      %dofmasses(1:3:end) = masses;
      %dofmasses(2:3:end) = masses;
      %dofmasses(3:3:end) = masses;
      %M = sparse( 1:numel(p), 1:numel(p), dofmasses );
       
      % GEVD, using lumped mass
      sigma = max(max(abs(K(ii,ii))));
      Kiisigma = K(ii,ii)/sigma;
      if isempty(Aeq)
        [Uii, Dii] = eigs( Kiisigma , M(ii,ii), nummodes+modeoffset, 'sm' ); 
      else
        neq = size(Aeq,1);
        [Uii, Dii] = eigs( [Kiisigma Aeq';Aeq sparse(neq,neq)], blkdiag(M(ii,ii),sparse(neq,neq)), nummodes+modeoffset, 'sm' ); 
        Uii = Uii(1:numel(ii),:);
      end
      vals = diag(Dii) * sigma;
          
      U = zeros( numel(p), nummodes+modeoffset );
      U(ii,:) = Uii; % let U be on all verts, including those pinned (if pinned)
  
      if computeModalDerivatives
        if pinned
          K_dec = decomposition(K(ii,ii));
          K_solve = @(rhs) K_dec \ rhs;
        else
          if isempty(Aeq)
            reg_method = 'null-orthogonal';
            switch reg_method
            case 'none'
              K_dec = decomposition(K(ii,ii));
              K_solve = @(rhs) K_dec \ rhs;
            case 'null-orthogonal'
              Q = [K(ii,ii) U(ii,:);U(ii,:)' sparse(size(U,2),size(U,2))];
              dec = decomposition(Q);
              K_solve = @(rhs) speye(numel(ii),numel(ii)+size(U,2)) * (dec \ [rhs; zeros(size(U,2),1)]);
            case 'tikonov'
              K_dec = decomposition(K(ii,ii) + 1e-4*M(ii,ii));
              K_solve = @(rhs) K_dec \ rhs;
            end
          else
            Q = [K(ii,ii) Aeq';Aeq sparse(neq,neq)];
            dec = decomposition(Q);
            K_solve = @(rhs) speye(numel(ii),numel(ii)+neq) * (dec \ [rhs; zeros(neq,1)]);
          end
        end
  
          BU = B*U;
           
          Hel = zeros( 9, 9, 9, n);
          for el = 1:n
           Hel( :, :, :, el) = stvk_d3(eye(3),lambda(el),mu(el));
          end
          for i = 1+modeoffset:nummodes+modeoffset
              for j = i:nummodes+modeoffset
                  % multiply out by the two different eigenvectors
                  %rhs = zeros( size(K,1), 1 );
                  %for el=1:n
                  %    Bel = B(el*9-8:el*9,:);
                  %    for k = 1:numel(F)
                  %        rhs = rhs + vol(el) * ( U(:,j)' * Bel' * H0(:,:,k) * Bel * U(:,i) ) * Bel(k,:)';
                  %    end
                  %end
                  BUi = reshape(BU(:,i),9,[]);
                  BUj = reshape(BU(:,j),9,[]);
                  BUir = repmat(reshape(BUi, [9, 1,  1,  153810]), [1 1 9 1]);

                  H0 = stvk_d3(eye(3),lambda(1),mu(1));
                  HelBUir = permute((squeeze(pagemtimes(Hel, BUir))), [1 3 2]);
%                   vH0BUi = permute(pagemtimes(H0,BUi),[3 1 2]).*permute(vol,[2 3 1]);
                  
                  vH0BUi = permute(HelBUir,[3 1 2]).*permute(vol,[2 3 1]);
                  
                  BUjvH0BUi = permute(sum(permute(BUj,[1 3 2]).*vH0BUi),[2 3 1]);
                  BTBUjvH0BUi = B'*BUjvH0BUi(:);
                  rhs = BTBUjvH0BUi;
  
                  % TODO: use LDLT of K(ii,ii) as we need it a bunch of times
                      Uij = zeros( numel(p), 1 );
                      Uij(ii) = K_solve( rhs(ii) ); % Equation 14 in Section 5.1 of Barbic 2005
                  %fprintf('%g ',norm( K * Uij - rhs, inf)/norm(rhs,inf));
                  %if ( pinned ) 
                  %else % see equaiton 29 of Barbic 2005
                  %    reg = 100; %How much regularizaiton for the free system?
                  %    %Uij = (K - reg*speye(size(K))) \ rhs; % hmm... see Appendix D of Barbic 2005
                  %    Uij = (K - reg*speye(size(K))) \ rhs; % hmm... see Appendix D of Barbic 2005
                  %end
                  %phi{i,j} = Uij / norm(Uij); % arbitrary normalization, should be mass-normalized (how?)
                  % Mass-normalized.
                  Uij = Uij/sqrt(Uij'*M*Uij);
                  phi{i,j} = Uij;
                  %fprintf('%g\n',[ Uij'*K*Uij])
              end
          end
      end
          
      if ShowAnim 
          fov = 30;
          cameraPosition = [ 11, -9, 3 ]*1.3;
          faceColor = [0.5 0.5 0.5];
          faceAlpha = 1;
  
          fh = figure(1); 
          if MakeVideo == 1
              v = VideoWriter(vidFileName,'MPEG-4');
              open(v);
          end
          clf;
          rh = 1/nummodes;        % ratio of plot height for plots
          rv = 1/(nummodes+1);    % ratio of plot width for plots
          for j = 1:nummodes
              rx = 0 + rh*(j-1);
              ry = 1-rv;
              sx = rh;
              sy = rv;
              ax(j) = axes('Position',[rx ry sx sy]);
              
              axis off;
              camva( fov );
              camlight('right', 'infinite');
              camlight('headlight');
              lighting flat;
              camproj('perspective');
              shading flat; 
              campos(cameraPosition);
              camlight('headlight'); % a second light?
              %     set(axesHandle,'CameraViewAngleMode','Manual');
              %     set(gca, 'Projection','perspective');
              pha(j) = patch('vertices', p, 'faces', faces, 'edgecol', [0.2, 0.2, 0.2], 'facecol', 'flat', 'FaceVertexCData', faceColor , 'FaceAlpha', faceAlpha, 'EdgeColor','none'); 
              axis equal
              axis manual;
              axis( [-6 6 -3 3 -3 3] )
          end
          if computeModalDerivatives
              for i = 1:nummodes
                  for j = i:nummodes
                      rx = 0 + rh*(j-1);
                      ry = 1-rv*(1+i);
                      sx = rh;
                      sy = rv;
                      axphi{i,j} = axes('Position',[rx ry sx sy]);
                      axis off;
                      camva( fov );
                      camlight('right', 'infinite');
                      camlight('headlight');
                      lighting flat;
                      camproj('perspective');
                      shading flat; 
                      campos(cameraPosition);
                      camlight('headlight'); % a second light?
                      %     set(axesHandle,'CameraViewAngleMode','Manual');
                      %     set(gca, 'Projection','perspective');
                      phphi{i,j} = patch('vertices', p, 'faces', faces, 'edgecol', [0.2, 0.2, 0.2], 'facecol', 'flat', 'FaceVertexCData', faceColor , 'FaceAlpha', faceAlpha, 'EdgeColor','none'); 
                      axis equal;
                      axis manual;
                      axis( [-6 6 -3 3 -3 3] )
                  end
              end
          end
          
          mag = 8;
          mag2 = 8;
          for theta = 1:4:720
              for j = 1:nummodes
                  mode = j + modeoffset; % offset will add 3 to avoid rigid modes 
                  pha(j).Vertices = p + mag*reshape(U(:,mode)',3,[])'*cosd(theta);
              end
              if computeModalDerivatives
                  for i = 1:nummodes
                      for j = i:nummodes
                          modei = i + modeoffset;
                          modej = j + modeoffset;
                          phphi{i,j}.Vertices = p + mag2*reshape(phi{modei,modej}',3,[])'*cosd(theta);
                      end
                  end
              end
              if MakeVideo == 1
                  figFrame = getframe(fh);
                  writeVideo(v,figFrame);
              else
                  drawnow
              end
          end
  
          if MakeVideo == 1
              close(v);
          end
      end
          
      % TODO: compute mass PCA on the vectors to get a complete basis!!
      basis = U(:,1+modeoffset:end);
      
      if GeneratePolys
          q = sym('q', [nummodes, 1] );
          assume( q, 'real' );
          psifun = matlabFunction(psi);
          dpsidFfun = matlabFunction(dpsidF);
          d2psidF2fun = matlabFunction(d2psidF2);
  
          Uq = basis * q;
          ppUq = reshape(p', [], 1) + Uq;
  
          [~,quarticTerms] = coeffs((sum(q)+1)^4);
          VCoeff = zeros( 1, numel(quarticTerms) );
          [~,cubicTerms] = coeffs((sum(q)+1)^3);
          RCoeff = zeros( numel(q), numel(cubicTerms) );
          [~,quadraticTerms] = coeffs((sum(q)+1)^2);
          KCoeff = zeros( numel(q)^2, numel(quadraticTerms) ); % symmetric, should just do half of them?
  
          for el = 1:n
              elBUiq = B(el*9-8:el*9,:) * ppUq;  % % displacement too!  not just BUq
              F1 = elBUiq(1);
              F2 = elBUiq(2);
              F3 = elBUiq(3);
              F4 = elBUiq(4);
              F5 = elBUiq(5);
              F6 = elBUiq(6);
              F7 = elBUiq(7);
              F8 = elBUiq(8);
              F9 = elBUiq(9);
  % 
  %             F1 = BUiq(el*9-8);
  %             F2 = BUiq(el*9-7);
  %             F3 = BUiq(el*9-6);
  %             F4 = BUiq(el*9-5);
  %             F5 = BUiq(el*9-4);
  %             F6 = BUiq(el*9-3);
  %             F7 = BUiq(el*9-2);
  %             F8 = BUiq(el*9-1);
  %             F9 = BUiq(el*9-0);
  
              BU = B(el*9-8:el*9,:) * basis;
  
              % reduced energy
              Vel = vol(el) * psifun( F1, F2, F3, F4, F5, F6, F7, F8, F9 ); 
              [C,T] = coeffs(Vel);
              assert( all( T == quarticTerms ) ); 
              VCoeff(1,:) = VCoeff(1,:) + double(C); 
  
              % reduced force
              Rel = - vol(el) * BU' * dpsidFfun( F1, F2, F3, F4, F5, F6, F7, F8, F9 ); % see this as dV/dF dF/dq where F = BUq
              for i = 1:numel(Rel)
                  [C,T] = coeffs(Rel(i));
                  assert( all( T == cubicTerms ) ); 
                  RCoeff(i,:) = RCoeff(i,:) + double(C); 
              end    
  
              if ( collectStiffnessPoly)
                  % reduced stiffness
                  Kel = - vol(el) * BU' * d2psidF2fun( F1, F2, F3, F4, F5, F6, F7, F8, F9 ) * BU;
                  for i = 1:numel(Kel)
                      [C,T] = coeffs(Kel(i));
                      %assert( all( T == quadraticTerms ) ); 
                      KCoeff(i,:) = KCoeff(i,:) + double(C); 
                  end    
              end
                  
          end
  
          Vfun = VCoeff * quarticTerms';
          ccode( Vfun, 'File', 'Vfun.c' );
          Rfun = RCoeff * cubicTerms';
          ccode( Rfun, 'File', 'Vfun.c' );
          Kfun = KCoeff * quadraticTerms';
          ccode( Kfun, 'File', 'Kfun.c' );
  
          % Might have slightly better efficiency to assemble all into one array
          % and code gen once?  Is Appendix B a more efficient approach?  There
          % should be some common subexperessions in these polynomials.
          % How would autodiff handle this?
      end
  end
end

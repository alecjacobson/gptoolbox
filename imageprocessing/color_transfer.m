function O = color_transfer(M,I)
  % O = color_transfer(M,I)
  %
  % Example-based Video Color Grading
  %
  % Inputs:
  %   M - reference image (e.g. a frame from Amélie )
  %   I - input image (e.g. a frame from the video to be graded)
  % Outputs:
  %   O - output image (same size as I)
  %

  function [W] = luminance_band_weights(L,k,overlap)
    % L is #pixels list of luminance values
    % k is the number of bands
    % overlap is a percentage of overlap

    %  To allow for this, we separate
    % each frame into three luminance bands with an equal number of
    % pixels. To avoid discontinuities, we make these bands overlap and
    % smoothly transition from one to the other. We achieve this by making
    % bands overlap by 10%, and setting pixels weights to 1 inside their
    % corresponding band with a cubic fallof outside. Finally, the per-pixel
    % weights across all the bands are normalized to 1

    [L_sorted,sort_idx] = sort(L(:));
    W = zeros(numel(L),k);
    % Pure band edges would be
    band_edges = linspace(0,1,k+1);
    band_edges = band_edges(1:end)+[-1;1]*overlap/2;
    band_edges = reshape(band_edges,[],1);
    U = zeros(numel(band_edges),k);
    for i = 1:k
      U((i-1)*(2)+(2:3),i) = 1;
    end
    band_edges(1) = 0;
    band_edges(end) = 1;
    U([1 end]) = 1;
    t = linspace(0,1,numel(L))';
    W = interp1(band_edges,U,t,'pchip');
    W(sort_idx,:) = W;
    W = W./sum(W(:));
  end

  function [mu,sigma] = weighted_mean_covariance(X,W)
    % X is #pixels by #bands
    % W is #pixels by #bands
    % mu is 1 by #bands
    % sigma is #bands by #bands
    mu = sum(X.*W,1)./sum(W,1);
    X_m = X-mu;
    sigma = (X_m'*(W.*X_m))./sum(W,1);
  end


  function O_ok = color_transfer_ok(M_ok,I_ok)
    I_L = I_ok(:,:,1);
    M_L = M_ok(:,:,1);
    O_L = histogram_matching(I_L,M_L,0.1);
    % number of bands
    k = 3;
    w = 0.1;
    I_W = luminance_band_weights(I_L,k,w);
    % #pixels by #bands
    M_W = luminance_band_weights(M_L,k,w);
    % compute weighted mean and weighted covariance for each band of the chrominance channels
    M_ab = reshape(M_ok(:,:,2:3),[],2);
    I_ab = reshape(I_ok(:,:,2:3),[],2);
    I_mu = zeros(k,2);
    I_sigma = zeros(2,2,k);
    M_mu = zeros(k,2);
    M_sigma = zeros(2,2,k);
    T = zeros(2,2,k);
    t = zeros(k,2);
    O_ab = zeros([size(I_ab) k]);
    for i = 1:k
      [I_mu(i,:),I_sigma(:,:,i)] = weighted_mean_covariance(I_ab,I_W(:,i));
      [M_mu(i,:),M_sigma(:,:,i)] = weighted_mean_covariance(M_ab,M_W(:,i));
      %imshow(matrixnormalize(reshape(I_W(:,i),size(I_ok,1),size(I_ok,2))));
      %pause

      sigma_m = M_sigma(:,:,i);
      sigma_i = I_sigma(:,:,i);
      % Linear Monge-Kantorovitch
      % T σᵢ Tᵀ = σₘ
      % T = √σᵢ⁻¹ √( √σᵢ σₘ √σᵢ ) √σᵢ⁻¹
      %
      % where √σ = Pᵀ √D P via eigendecomposition
      %T(:,:,i) = sqrtm(sigma_m) / sqrtm(sigma_i);
      % This is more accurate in terms of σ_o == σₘ but not seeing the
      % difference in the images.
      sqrt_sigma_i = sqrtm(sigma_i);
      T(:,:,i) = sqrt_sigma_i \ sqrtm(sqrt_sigma_i * sigma_m * sqrt_sigma_i) / sqrt_sigma_i;
      
      %T(:,:,i)-T(:,:,i)'

      t(i,:) = M_mu(i,:) - I_mu(i,:) * T(:,:,i);

      O_ab(:,:,i) = I_ab*T(:,:,i) + t(i,:);
      %[O_mu(i,:),O_sigma(:,:,i)] = weighted_mean_covariance(O_ab(:,:,i),I_W(:,i));
      %[I_mu(i,:);M_mu(i,:);O_mu(i,:)]
      %[I_sigma(:,:,i);M_sigma(:,:,i);O_sigma(:,:,i)]
    end
    %for i = 3
    %  clf;
    %  hold on;
    %  sct(M_ab,'filled','SizeData',matrixnormalize(M_W(:,i))*1+1e-7);
    %  sct(I_ab,'filled','SizeData',matrixnormalize(I_W(:,i))*1+1e-7);
    %  sct(O_ab(:,:,i),'filled','SizeData',matrixnormalize(I_W(:,i))*1+1e-7);
    %  sct(I_mu(i,:),'b','filled','SizeData',100);
    %  sct(M_mu(i,:),'r','filled','SizeData',100);
    %  sct(O_mu(i,:),'g','filled','SizeData',100);
    %  % Draw ellipses for M_sigma(:,:,i) and I_sigma(:,:,i)
    %  th = linspace(0,2*pi,100);
    %  [V,D] = eig(M_sigma(:,:,i));
    %  X = V*sqrt(D)*[cos(th);sin(th)];
    %  plot(X(1,:)+M_mu(i,1),X(2,:)+M_mu(i,2),'-r','LineWidth',2);
    %  [V,D] = eig(I_sigma(:,:,i));
    %  X = V*sqrt(D)*[cos(th);sin(th)];
    %  plot(X(1,:)+I_mu(i,1),X(2,:)+I_mu(i,2),'-b','LineWidth',2);
    %  [V,D] = eig(O_sigma(:,:,i));
    %  X = V*sqrt(D)*[cos(th);sin(th)];
    %  plot(X(1,:)+O_mu(i,1),X(2,:)+O_mu(i,2),'-g','LineWidth',2);
    %  hold off;
    %end

    O_ab = sum(O_ab .* permute(I_W,[1 3 2]),3) ./ sum(I_W,2);
    O_ok = cat(3,O_L,reshape(O_ab,size(I_ok,1),size(I_ok,2),2));

    
    %imshow(oklab2rgb(cat(3,O_L,I_ok(:,:,2:3))));
    %imshow(oklab2rgb(O_ok));
    %size(I_ok)
    %size(I_W)
    %imshow([oklab2rgb(I_ok) repmat(reshape(matrixnormalize(I_W(:,3)),size(I_ok,1),size(I_ok,2)),[1 1 3])]);
    %clf;
    %hold on;
    %plot(sort(O_L(:)),linspace(0,1,numel(O_L))','-');
    %plot(sort(I_L(:)),linspace(0,1,numel(I_L))','-');
    %plot(sort(M_L(:)),linspace(0,1,numel(M_L))',':');
    %drawnow

  end

  M = im2double(M);
  I = im2double(I);
  M_ok = rgb2oklab(M);
  I_ok = rgb2oklab(I);
  O_ok = color_transfer_ok(M_ok,I_ok);
  O = oklab2rgb(O_ok);
end

function [denoised, noise] = denoise(noisy_signal,prior,std_noise)
  % Given a noisy image and a similar prior image and the standard deviation of
  % the noise (assumes gaussian distribution) denoise the iamge by building a
  % steerable freq pyramid and using probabilistic algorithms to build a coring
  % function (highest 5 bands only) to remap the image to the original. Also
  % prints and shows plots of the estimators.
  %
  % [denoised, noise] = denoise(noisy_signal,prior,std_noise)
  %
  % Inputs:
  %   noisy_signal  image with noise
  %   prior  similar image without noise 
  %   std_noise  the standard deviation of the noise distribution
  % Output:
  %   denoised  the denoised result image
  %   noise  the generated noise image
  %
  
  % create noise image
  noise = std_noise.*randn(size(noisy_signal));

  % build steerable pyramids for each image
  [prior_pyr,prior_ind]=buildSFpyr(prior,3,3);
  [noisy_signal_pyr,noisy_signal_ind]=buildSFpyr(noisy_signal,3,3);
  [noise_pyr,noise_ind]=buildSFpyr(noise,3,3);

  % reconstruct right away (for debugging)
  %before = reconSFpyr(noisy_signal_pyr,noisy_signal_ind);

  % open a new figure for plotting the estimators
  %screen_size = get(0, 'ScreenSize');
  estimator_figure = figure;
  %set(estimator_figure, 'Position', [0 0 screen_size(3) screen_size(4) ] );
  set(estimator_figure, 'Position', [0 0 1000 200 ] );

  for band_index = 1:5
    % extract top five bands of the prior and noise pyramids
    noise_band = pyrBand(noise_pyr,noise_ind,band_index);
    prior_band = pyrBand(prior_pyr,prior_ind,band_index);
    noisy_signal_band = pyrBand(noisy_signal_pyr,noisy_signal_ind,band_index);

    % small value to add to bins to avoid mayhem
    %EPSILON = 0.1;
    EPSILON = 0.000001;
    % histogram the coefficients 
    histogram_bins = -0.3:0.01:0.3;
    histogram_prior = EPSILON+ ...
      hist(reshape(prior_band,prod(size(prior_band)),1),histogram_bins);
    histogram_noise = EPSILON+ ...
      hist(reshape(noise_band,prod(size(noise_band)),1),histogram_bins);
    % compute estimate of bayesian estimator using formula from simoncelli96c.pdf:
    %        ∫ Pn(y-z) z Px(z) dz
    % b(y) = ---------------------
    %         ∫ Pn(y-z) Px(z) dz 
    %
    % where Pn is the pdf of noise and Px of image. We estimate Pn with
    % histogram_noise and Px with histogram_prior. Then both top and bottom are
    % convolutions ( z is histogram_bins)
    %
    bayesian_estimator = ...
      conv2(histogram_noise,histogram_bins.*histogram_prior,'same') ./ ...
      conv2(histogram_noise,histogram_prior,'same');

    % plot  the estimator with a straight line for comparison
    subplot(1,5,band_index);
    plot(histogram_bins,bayesian_estimator,histogram_bins,histogram_bins);
    grid on
    %legend('bayesian estimator','straight line for comparison');
    
    if(band_index==1)
      title('hifi residual');
    else
      title(['hifi orientation ' num2str(band_index-1)])
    end

    % (re)map signal coefficients using bayesian_estimator 
    denoised_band = interp1(histogram_bins,bayesian_estimator,noisy_signal_band);

    % add back to pyramid
    noisy_signal_pyr = ...
      setPyrBand(noisy_signal_pyr,noisy_signal_ind,denoised_band,band_index);

  end

  % print the plots
  estimator_figure;
  print('-depsc','bayesian_estimator_plots')

  % reconstruct image from modified pyramid
  denoised = reconSFpyr(noisy_signal_pyr,noisy_signal_ind);

end

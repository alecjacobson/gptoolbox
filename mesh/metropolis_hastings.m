function [S,F] = metropolis_hastings(unnorm_distr,next_sample,num_samples,x0)
% Uses the Metropolis-Hastings algorithm to generate a sequence of random samples from an
% unknown distribution p assuming one knows a function f whcih is proportional to p.
%
% [S,F] = metropolis_hastings(unnorm_distr,next_sample,num_samples,x0)
%
% Inputs:
%       unnorm_distr function handle returning the value of the known
%           function which is proportional to the desired distribution
%
%           [f,data_ud] = unnorm_distr(x,data_d)
%
%           Inputs:
%               x #dim vector of sample being considered
%               data_ud empty on first call, or `data_ud` output from prev call
%           Outputs:
%               f scalar value of f
%               data_ud  persistent callback data
%
%       next_sample function handle returning a candidate next sample from 
%           the current sample x; for example, Gaussian centered at x
%
%           [x1,data_ns] = next_sample(x,data_ns)
%
%           Inputs:
%               x #dim vector of previous sample
%               data_ns empty on first call, or `data_ns` output from prev call
%           Outputs:
%               x1 next candidate sample
%               data_ns  persistent callback data
%
%       num_samples integer number of samples in output (this will be more
%           than the total number of considered samples)
%
%       x0 #dim vector of first element in sequence of samples
%
% Outputs:
%       S #num_samples by #dim sequence of samples
%       F #num_samples vector of f evaluated at each row of S
%
%
% Example (run in script to include function definitions):
%
%
% [S,F] = metropolis_hastings(@unnorm_distr,@next_sample,100000,0.1);
% histogram(S) % Should look like a graph of f (absolute value)
% 
% 
% function [x1,data] = next_sample(x0,data)
%     % common choice: normal centered at x0
%     x1 = normrnd(x0,0.1);
% end
% 
% function [f,data] = unnorm_distr(x,data)
%     % f is any function without need for normalization
%     f = max(1 - abs(x),1e-5);
% end

% Data structs
data_ud = [];
data_ns = [];

% Initialize outputs to zero
S = zeros(num_samples,size(x0,2));
F = zeros(num_samples,1);

% Initialize x0 and f0
S(1,:) = x0;
[f0,data_ud] = unnorm_distr(x0,data_ud);

sample_num = 2; % Start from 2
while sample_num < num_samples
    % Compute next candidate sample
    [x1,data_ns] = next_sample(x0,data_ns);
    % Compute f1 at new candidate
    [f1,data_ud] = unnorm_distr(x1,data_ud);
    % Generate random value between 0 and 1
    r = rand(1,1);
    if r<(f1/f0) 
       % If accepted, update current sample with candidate sample and
       % add to running vectors S and F
       x0 = x1;
       f0 = f1;
       S(sample_num,:) = x1;
       F(sample_num,:) = f1;
       sample_num = sample_num + 1; % samples counter
    end
end

end
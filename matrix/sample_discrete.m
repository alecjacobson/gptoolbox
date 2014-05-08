function [S,D,J] = sample_discrete(D,n)
  % SAMPLE_DISCRETE  Sample a discrete distribution given by D n times. This
  % employs the "Alias Table Method". Ensuring that the running time is O(n +
  % m)
  %
  % Inputs:
  %   D  m list of probabilities summing to 1
  %   n  number of samples to output
  % Outputs:
  %   S  n list of indices into D
  %   D  m list of alias table heights
  %   J  m list of "upper half" reindices
  %   

  % Assume column vector
  D = D(:);
  % https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
  m = numel(D);
  assert(abs(sum(D) - 1)<1e-15);
  % We want to sort based on above or below 1/n
  D = D*m;
  % now we can sort based on above or below 1
  % Below queue: http://www.alecjacobson.com/weblog/?p=3933
  % and preallocate up to 2*m (extra m to be safe, I don't think it's needed)
  B = find(D<1);
  Bl = numel(B);
  B(end+1:2*m) = 0;
  % After queue, and preallocate up to m
  A = find(D>=1);
  Al = numel(A);
  A(end+1:m) = 0;
  % Upper half list (initialize with this index because above includes "equal")
  J = (1:m)';

  % While below and above still have elements
  while Bl>0 && Al>0
    % pop below
    b = B(Bl);
    Bl = Bl - 1;
    % Look above
    a = A(Al);
    % a going to be b's upper half
    J(b) = a;
    % lob off enough to fit D(b) to 1
    %D(a) = D(a) - (1.0 - D(b));
    % Reduce rounding error (see comments above)
    D(a) = (D(a) + D(b)) - 1.0;
    if D(a) < 1
      % pop off above
      Al = Al - 1;
      % push onto below
      Bl = Bl+1;
      B(Bl) = a;
    % else still on above
    end
  end

  R1 = floor(rand(n,1)*m)+1;
  R2 = rand(n,1);
  I = (1:m)';
  % select those that are in upper half
  upper = R2>=D(R1);
  S = I(R1);
  S(upper) = J(R1(upper));
  
end

function a = extract_arg(j,func,varargin)
  % a = extract_arg(j,func,…)
  %
  % Inputs:
  %   j  index of the argument to extract (1-based)
  %   func  function handle to call
  %   …  additional arguments to pass to the function
  % Outputs:
  %   a  the j-th output of calling func(…)
  %

  switch j
    case 7
      [~,~,~,~,~,~,a] = func(varargin{:});
    case 6
      [~,~,~,~,~,a] = func(varargin{:});
    case 5
      [~,~,~,~,a] = func(varargin{:});
    case 4
      [~,~,~,a] = func(varargin{:});
    case 3
      [~,~,a] = func(varargin{:});
    case 2
      [~,a] = func(varargin{:});
    case 1
      a = func(varargin{:});
  end
end


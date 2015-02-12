function b1 = betti(F,varargin)
  % BETTI Comput the first betti number of simplicial mesh
  % All of the implemented methods are slow. The smith-normal-form seems like
  % the robust thing to do, but is cubic in the size of the mesh. Perhaps
  % preprocessing with topology preserving edge-(facet?)-collapses would
  % drastically reduce run-time.
  % 
  % b1 = betti(F)
  % b1 = betti(F,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   F  #F by simplex-size mesh indexing some vertex set
  %   Optional:
  %     'Method'  followed by one of the following:
  %       'smith-normal-form' robust implementation using mod 2 integer
  %         arithmatic (actually logicals and XOR for performance). Correct
  %         (robust), slow.
  %       'rank'  Sparse rank computation on boundary map (incidence matrices)
  %       'full-rank' Full rank 〃
  %       'qr'  Use qr to compute rank
  %       {'eigenanalysis'} Find number of zeros in eigenvalues of discrete
  %         Laplacian.
  % Outputs:
  %   b1  first betti number (number of loops)
  %      
  %   

  % default values
  method = 'eigenanalysis';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method'}, ...
    {'method'});
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

  n = max(F(:));

  % For MANIFOLD! triangle meshes this should just use Euler formula like
  % `statistics.m` does.

  % There seems to be two competing views on the "boundary map" or incidence
  % matrix: whether it is oriented (directed) or not:
  %
  % Undirected:
  % https://www.cs.duke.edu/courses/fall06/cps296.1/Lectures/sec-IV-3.pdf
  %
  % Directed:
  % http://comptop.stanford.edu/u/programs/plex/plexintro.pdf
  %
  % Q: Is the Smith Normal Form computation in intergers modulo 2 effectively
  % treating these incidence matrices as signed (directed)?

  E = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  [sE,sJ] = sort(E,2);
  [uE,IA,EMAP] = unique(sE,'rows');


  switch method
  case 'smith-normal-form'
    % Undirected:
    % https://www.cs.duke.edu/courses/fall06/cps296.1/Lectures/sec-IV-3.pdf
    % boundary matrices
    D1 = sparse(  uE,repmat(1:size(uE,1),2,1)',1,n,size(uE,1));
    D2 = sparse(EMAP,repmat(1:size(F,1),3,1)',1,size(uE,1),size(F,1));
    S1 = snf(D1);
    [rB0,~] = find(S1,1,'last');
    rZ1 = size(S1,2) - rB0;
    S2 = snf(D2);
    [rB1,~] = find(S2,1,'last');
    b1 = rZ1 - rB1;
  case {'eigenanalysis','rank','full-rank','qr'}
    % Directed:
    % sJ reveals whether edge was flipped from E to sE (and thus also uE)
    D1 = sparse(  uE,repmat(1:size(uE,1),2,1)',repmat([-1 1],size(uE,1),1),n,size(uE,1));
    D2 = sparse(EMAP,repmat(1:size(F,1),3,1)',sJ(:,1)*2-3,size(uE,1),size(F,1));

    switch method
    case {'rank','full-rank','qr'}
      switch method
      case 'rank'
        dB0 = find(spnull(D1'),1,'last');
        dB1 = find(spnull(D2),1,'last');
        % somehow this is not finding same rank as rank(full())... off by
        % one...
        dB0 = dB0-1;
        dB1 = dB1-1;
      case 'full-rank'
        dB0 = rank(full(D1));
        dB1 = rank(full(D2));
      case 'qr'
        % Q: Is rank(S1) just rank(D1)-1 ?
        % Seems to be as easy as:
        [Q,R,E] = qr(D1');
        [~,dB0] = find(R',1,'last');
        [Q,R,E] = qr(D2');
        [~,dB1] = find(R',1,'last');
      end
      dZ1 = size(uE,1)-dB0;
      b1 = dZ1 - dB1;
    case 'eigenanalysis'
      % http://comptop.stanford.edu/u/programs/plex/plexintro.pdf
      % "discrete Laplacian"
      L1 = D1'*D1 + D2*D2';
      k = min(10,size(uE,1)-1);
      while true
        [EV,ED] = eigs(-L1,k,'sm');
        b1 = sum(diag(abs(ED))<eps);
        if b1 < k
          break;
        end
        k = k*2;
      end
    end
  otherwise
    error('Unknown method: %s',method);
  end

end

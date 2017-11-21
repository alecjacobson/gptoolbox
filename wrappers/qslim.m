function [C,W,G,V2W] = qslim(V,F,t,varargin)
  % QSLIM Very simple wrapper for qslim
  %
  % [C,W,G,V2W] = qslim(V,F,t,'ParameterName',ParameterValue)
  %
  % Inputs:
  %   V  #V by 3 input mesh vertex positions
  %   F  #F by 3 input mesh triangle indices (1-indexed)
  %   t  target number of faces {0 for all collapse}
  %   Optional:
  %     'QSlimFlags' followed by qslim flags
  % Outputs:
  %   C  #collapse list of edge collapse structs with fields
  %     .a index into V of vertex collapsed to
  %     .b index into V of vertex collapsed from
  %     .da  1 by 3 displacement of vertex a
  %     .rm list of indices in F of faces removed by this collapse
  %     .mod list of indices in F of faces modified by this collapse
  %   W  #W by 3 collapsed mesh vertex positions
  %   G  #G by 3 collapsed mesh triangle indices (1-indexed)
  %   V2W  #V list of indices into W
  %   
  % See also: readLOG, perform_edge_collapse, path_to_qslim
  %

  qslim_flags = '';
  if nargin <3
    t = 0;
  end

  v = 1;
  while v<=numel(varargin)
    switch varargin{v}
    case 'QSlimFlags'
      assert((v+1)<=numel(varargin));
      v = v+1;
      qslim_flags = varargin{v};
    otherwise
      error(['Unsupported parametername: ' varargin{v}]);
    end
    v = v+1;
  end

  prefix = tempprefix();
  input = [prefix '.smf'];
  output = [prefix '.log'];

  % Write input to file
  writeSMF(input,V,F);

    % prepare command string
    command = sprintf('%s -t %d -M log -q %s %s >%s', ...
      path_to_qslim,t,qslim_flags,input,output);
  try
    [status,result] = system(command);
    if status ~= 0
      error(result);
    end

    [VV,FF,C] = readLOG(output);
    % should match input
    assert(size(VV,1)==size(V,1));
    assert(size(FF,1)==size(F,1));
    assert(isequal(F,FF));

    delete(input);
    delete(output);
  catch ME
    fprintf('%s\n',command);
    rethrow(ME);
  end

  % This is expensive so only perform if required by output
  if nargout>1
    [W,G,V2W] = perform_edge_collapse(V,F,C,t);
  end

end

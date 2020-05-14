function varargout = cache(fun,varargin)
  % CACHE Rather than call an expensive function again, load its output from the
  % last time it was called with the same input.
  %
  % Replace
  %   [out1,out2, ...] = myfun(in1,in2, ...)
  % with
  %   [out1,out2, ...] = cache(@myfun,in1,in2, ...)
  % 
  % Everytime `[out1,out2, ...] = cache(@myfun,in1,in2, ...)` is called, it will
  % look for a cached .mat file containing the output data `[out1,out2,...]`. If
  % it exists, then it's loaded and returned. If not, then `[out1,out2, ...] =
  % myfun(in1,in2, ...)` and the outputs are cached into a .mat file.
  %
  % There's an extremely unlikely possibility of cache hash collision (md5).
  %
  % If the number of outputs changes a new cache will be generated.
  %
  tmpf = [tempname('.') '.mat'];
  NARGOUT = nargout;
  save(tmpf,'fun','varargin','NARGOUT');
  % get md5 checksum on input "state" to use as filename for cache
  [s,cache_name] = system([ ...
    'dd  if=' tmpf ' 2>/dev/null skip=128 bs=1 | /usr/local/bin/md5sum | ' ...
    'awk ''{printf "cache."$1".mat"}''']);
  % clean up
  delete(tmpf);
  clear s tmpf;
  % If the checksum cache file exists then we've seen this input "state"
  % before, load in cached output "state"
  if(exist(cache_name,'file'))
    warning('Loading %s...',cache_name);
    % use cache
    load(cache_name);
  % Otherwise this is the first time we've seen this input "state", so we
  % execute the function as usual and save the output "state" to the cache 
  else
    % Deal with unknown variable number of outputs
    eval(['[' sprintf('varargout{%d} ',1:nargout) '] = fun(varargin{:});']);
    save(cache_name,'varargout');
  end
end

function C = cache_test(A,B)
  % CACHE_TEST Dummy program that computes C = A+B, but demonstrates how to use
  % md5 caching of function results on input parameters
  %
  % C = cache_test(A,B)
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Check for cached result, do NOT edit variables until cache is checked,
  % your function code comes later. See below
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get a list of current variables in this scope, this is the input "state"
  variables = who;
  % get a temporary file's name
  tmpf = [tempname('.') '.mat'];
  % save the "state" to file, so we can get a md5 checksum
  save(tmpf,'-regexp',sprintf('^%s$|',variables{:}),'-ascii');
  % get md5 checksum on input "state", we append .cache.mat to the check sum
  % because we'll use the checksum as the cache file name
  [s,cache_name] = system(['md5sum ' tmpf ' | awk ''{printf "."$1".cache.mat"}''']);
  % clean up
  delete(tmpf);
  clear s tmpf variables;

  % If the checksum cache file exists then we've seen this input "state"
  % before, load in cached output "state"
  if(exist(cache_name,'file'))
    fprintf('Cache exists. Using cache...\n');
    % use cache
    load(cache_name);
  % Otherwise this is the first time we've seen this input "state", so we
  % execute the function as usual and save the output "state" to the cache 
  else
    fprintf('First time. Creating cache...\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Your function code goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = A+B;

    % get list of variables present in this scope at finish of function code,
    % this is the output "state"
    variables = who;
    % save output "state" to file, using md5 checksum cache file name
    save(cache_name,'-regexp',sprintf('^%s$|',variables{:}));
  end
end

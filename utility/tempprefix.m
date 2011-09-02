function T = tempprefix(D);
  % TEMPPREFIX returns a unique name prefix, starting with the given directory
  %   suitable for use as a prefix for temporary files.
  %
  % T = tempprefix(D)
  % T = tempprefix()  same as above but uses builtin tempdir as D
  %
  % Inputs:
  %   D  name of directory to find prefix
  % Outputs:
  %   T  unique temp prefix
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: tempname
  % 

  while(true)
    if(exist('D','var'))
      T = tempname(D);
    else
      T = tempname();
    end
    if(isempty(dir([T '*'])))
      break;
    end
  end
end

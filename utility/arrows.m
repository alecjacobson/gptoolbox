function ar = arrows(X,U,varargin)
  % ar = arrows(X,U,varargin)
  %
  % Example:
  %   % Drop in replacement for qvr(X,U,0) ...
  %   % Use qvr (or anything) to the axis.
  %   qvr(X,U,0,'Color','none');
  %   arrows(X,U,'HeadStyle','vback1','LineWidth',1,'HeadWidth',4,'HeadLength',4);
  ar = arrayfun(@(i)  ...
    setfield(annotation('arrow','position',[X(i,:) U(i,:)],varargin{:}),'Parent',gca), ...
    (1:size(X,1))');
end

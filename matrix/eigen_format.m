function varargout = eigen_format(A,name)
  format = ['  ' repmat('%0.17g,',1,size(A,2)) '\n'];
  if nargin > 1
      if islogical(A)
        t = '<bool,Eigen::Dynamic,Eigen::Dynamic>';
      else
        if size(A,1) == 1
          t= 'RowVectorX';
        elseif size(A,2) == 1
          t = 'VectorX';
        else
          t = 'MatrixX';
        end
        if all(A(:) == round(A(:)))
          t = [t 'i'];
        else
          t = [t 'd'];
        end
      end

    s = {sprintf('Eigen::%s %s(%d,%d);',t,name,size(A,1),size(A,2))};
    s{end+1} = [name '<< '];
  else
    s = {};
  end
  s{end+1} = sprintf(format,A');
  s{end}(end-1:end) = sprintf(';\n');
  s{end}(end) = [];

  if nargout >= 1
    varargout{1} = strjoin(s,sprintf('\n'));
  else
    fprintf('%s',strjoin(s,sprintf('\n')));
  end

end

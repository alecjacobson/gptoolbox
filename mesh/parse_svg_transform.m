function T = parse_svg_transform( t_str )
  
  pattern = '(scale|translate|rotate|skewX|skewY)\((-?(\d+\.\d*|\.\d+|\d+)(\s*[ ,]\s*)?){1,3}\)|(matrix)\((-?(\d+\.\d*|\.\d+|\d+)(\s*[ ,]\s*)?){6}\)';
  [startIdx, endIdx, tokens] = regexp(t_str, pattern, 'start', 'end', 'tokens');
  T = eye(3,3);
  for ti = 1:numel(tokens)
    token = tokens{ti};
    switch token{1}
    case 'matrix'
      [Ti,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      assert(foundT == 6);
      Ti = [reshape(Ti,2,3);0 0 1];
    case 'translate'
      [s,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      switch foundT
      case 2
        Ti = [1 0 s(1); 0 1 s(2); 0 0 1];
      case 1
        Ti = [1 0 s(1); 0 1 0; 0 0 1];
      otherwise
        error('scale transform must have 1 or 2 parameters');
      end
    case 'skewX'
      [s,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      a = s(1)/180*pi;
      Ti = [1 tan(a) 0;0 1 0;0 0 1];
    case 'skewY'
      [s,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      a = s(1)/180*pi;
      Ti = [1 0 0;tan(a) 1 0;0 0 1];
    case 'scale'
      [s,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      switch foundT
      case 2
        Ti = [s(1) 0 0; 0 s(2) 0; 0 0 1];
      case 1
        Ti = [s(1) 0 0; 0 s(1) 0; 0 0 1];
      otherwise
        error('scale transform must have 1 or 2 parameters');
      end
    case 'rotate'
      [s,foundT] = sscanf(strrep(token{2},',',' '),'%g');
      a = s(1)/180*pi;
      Ti = [[cos(a) -sin(a) 0;sin(a) cos(a) 0];0 0 1];
      if foundT == 3
        X = [1 0 s(2);0 1 s(3);0 0 1];
        Ti = X * Ti * inv(X);
      end
    end
    T = T * Ti;
  end
  T = T(1:2,:);
end

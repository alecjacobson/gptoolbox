function D = cubic_derivative_eval(C,t,k)
  t = reshape(t,numel(t),1);
  switch k
  case 0
    D = cubic_eval(C,t);
  case 1
    D = 3.*(1-t).^2.*(C(2,:)-C(1,:)) + 6*(1-t).*t.*(C(3,:)-C(2,:)) + 3.*t.^2.*(C(4,:)-C(3,:));
  case 2
    D = 6.*(1-t).*(C(3,:) - 2*C(2,:) + C(1,:)) + 6.*t.*(C(4,:) - 2*C(3,:) + C(2,:));
  case 3
    D = -6.*C(1,:) + 18.*C(2,:) - 18.*C(3,:) + 6.*C(4,:);
  otherwise
    D = zeros(numel(t),size(C,2));
  end
end

function [param,mosek_exists] = default_mosek_param()
  warning('deprecated. please call default_quadprog_param()');
  [param,mosek_exists] = default_quadprog_param();
end

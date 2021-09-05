function [step_size,X,obj,step,iter,max_iter] = line_search(obj_fun,proj_fun,X0,step0,max_step_size)
  % LINE_SEARCH Given an objective function (obj_fun) to minimize and
  % constraints to project onto (proj_fun) and an intial guess (X0) determine a
  % feasible step_size to move along a desired direcxtion (step0)
  %
  % [step_size,X,obj,step] = line_search(obj_fun,proj_fun,X0,step0,max_step_size)
  %
  % Inputs:
  %   obj_fun  objective function taking X as input
  %   proj_fun constraints projection function taking X as input
  %   X0  #X by 1 list of initial values
  %   step0  #X by 1 desired step direction
  %   max_step_size  maximum step size to consider
  % Outputs:
  %   step_size  scalar safe step size
  %   X  #X by 1 new values (X = proj_fun(X0+step_size*step0))
  %   obj  scalar objective value at X
  %   step  #X by 1 step conducted before projection (step_size*step0)
  %
  
  step_size = max_step_size;
  obj0 = obj_fun(proj_fun(X0));
  max_iter = 20;
  for iter = 1:max_iter
    step = step_size*step0;
    X = proj_fun(X0+step);
    obj  = obj_fun(X);
    if obj < obj0
      return;
    else
      step_size = step_size*0.5;
    end
  end
  step_size = 0;
  X = X0;
  obj = obj0;
  step = zeros(size(step));
end

function [x,fval] = pso(fcn,lb,ub,varargin)
  % [x,fval] = pso(fcn,lb,ub,varargin)
  nvars = numel(lb);
  lb = reshape(lb,1,[]);
  ub = reshape(ub,1,[]);
  max_iters = 10;
  population = 200;
  position = lb + rand(population,nvars).*(ub-lb);
  randn = @(m) (2*rand(population,m)-1);
  velocity = 0.1*(ub-lb).*randn(nvars);
  
  fcn_all = @(position) ...
    cell2mat(arrayfun( ...
      @(i) fcn(position(i,:)), ...
      1:population, ...
      'UniformOutput',0));

  best_f = fcn_all(position);
  best_position = position;
  [fval,xi] = min(best_f);
  x = best_position(xi,:);

  %Matryoshka
  %omega = 0.98;
  %phi_p = 0.01;
  %phi_g = 0.01;
  %Cylinder fit
  omega = 0.7;
  phi_p = 0.15;
  phi_g = 0.15;
  for iter = 1:max_iters
    velocity = ...
      omega * velocity + ...
      phi_p * randn(1) .* (best_position - position) + ...
      phi_g * randn(1) .* (x             - position);

    reflection = false;
    if reflection
      new_position = position+velocity;
      velocity(new_position<lb & velocity<0) = -velocity(new_position<lb & velocity<0);
      velocity(new_position>ub & velocity>0) = -velocity(new_position>ub & velocity>0);
      position = new_position;
    else
      position = min(max(lb,position),ub);
    end
    

    %switch nvars
    %case 2
    %  clf;
    %  hold on;
    %  scatter(position(:,1),position(:,2),'.');
    %  %scatter(best_position(:,1),best_position(:,2),'.');
    %  %scatter(x(:,1),x(:,2),'o');
    %  quiver(position(:,1),position(:,2),velocity(:,1),velocity(:,2));
    %  hold off;
    %  axis equal;
    %  axis([lb(1) ub(1) lb(2) ub(2)]);
    %  drawnow;
    %case 3
    %  clf;
    %  hold on;
    %  scatter3(position(:,1),position(:,2),position(:,3)'.');
    %  %scatter(best_position(:,1),best_position(:,2),'.');
    %  %scatter(x(:,1),x(:,2),'o');
    %  quiver3(position(:,1),position(:,2),position(:,3),velocity(:,1),velocity(:,2),velocity(:,3));
    %  hold off;
    %  axis equal;
    %  axis([lb(1) ub(1) lb(2) ub(2) lb(3) ub(3)]);
    %  view(15,15);
    %  camproj('persp');
    %  drawnow;
    %end

    f = fcn_all(position);
    improved = f<best_f;
    best_position(improved) = position(improved);
    [fval,xi] = min(best_f);
    x = best_position(xi,:);
  end
  
end

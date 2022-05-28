function [U,t,t1] = rest_on_ground(V,F,varargin)
  % [U,t] = rest_on_ground(V,F)
  % 
  % Given a solid mesh (with a meaningful centroid) rest the object on the
  % ground (z=0). This effectively simulates dropping the object on the ground
  % in a universe with no momentum, or imagine the object coming to rest on the
  % floor of a viscous ocean.
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into rows of V
  %   Optional:
  %     'Centroid' followed by center of mass
  % Outputs:
  %   U  #U by 3 list of new mesh vertex positions
  %   t  3-long list of contact points
  %

  % default values
  cen = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Centroid'}, ...
    {'cen'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  warning('Should reduce V,F to convex hull...');

  %save('bad-rest.mat','V','F')
  if isempty(cen)
    cen = centroid(V,F);
  end
  U = [V;cen];
  t = [];
  stable = false;
  cross2 = @(a,b,c) ...
    [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
     a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];

  iter = 0;
  max_iters = 3*size(F,1);
  t1 = [];
  while ~stable
    iter = iter+1;
    if iter > max_iters
      warning('Failed to converge in %d iterations',max_iters);
      break;
    end
    switch(numel(t))
    case 0
      [z,t] = min(U(:,3));
      U(:,3) = U(:,3)-z;
      assert(all(isreal(U(:))));
    case 1
      cen = U(end,:);
      arm = normalizerow(cen-U(t,:));
      up = [0 0 1];
      ax = normalizerow(cross(arm,up));
      bi = cross(up,ax);
      W = (U-U(t,:))*[bi;up]'; 
      A = atan2(W(:,2),W(:,1));
      Acen = A(end);
      if Acen>pi/2
        A = pi-A;
        ax = -ax;
      end
      A([t end]) = inf;
      A(A<0) = A(A<0)+2*pi;
      [z,s] = min(A);
      R2 = axisangle2matrix(ax,z);

    %clf;
    %hold on;
    %tsh = tsurf(F,U,'CData',A,fphong,falpha(0.1,0),fsoft);
    %ssh =     sct(U(end,:),'o','MarkerFaceColor','r');
    %ssh =     sct(U(t,:),'o','MarkerFaceColor','g');
    %ssh =     sct(U(s,:),'o','MarkerFaceColor','k');
    %%ssh = sct(U(g,:),'o','MarkerFaceColor','r');
    %qvr(U(t,:),100*bi,'Color','r','LineWidth',3);
    %qvr(U(t,:),100*up,'Color','g','LineWidth',3);
    %qvr(U(t,:),100*ax,'Color','b','LineWidth',3);
    %qvr(U(t,:),100*arm,'Color','c','LineWidth',3);
    %%hold off;
    %%l = light('Position',[0 1e-10 1],'Style','infinite');
    %%add_shadow(tsh,l,'Fade','none');
    %%sh = add_shadow([],l);
    %%sh{1}.MarkerFaceColor = ssh.MarkerFaceColor;
    %axis equal;
    %view(0,0);
    %hold off;
      

      U = (U-U(t,:))*R2+U(t,:);
      t = [t s];

      assert(all(isreal(U(:))));

    case 2
      cen = U(end,:);
      s = t(2);
      ax = normalizerow(U(s,:)-U(t(1),:));
      W = (U-U(s,:))*[normalizerow(cross(up,ax));up]';
      A = atan2(W(:,2),W(:,1));
      Acen = A(end);
      if Acen>pi/2
        A = pi-A;
        ax = -ax;
      end
      A([t end]) = inf;
      % ignore colinear vertices
      D = normrow(cross2(ax,U-U(t(1),:),2));

      A(D<1e-10) = inf;
      [z,g] = min(A);
      R3 = axisangle2matrix(ax,z);
      if all(t == [1779 8746])
    %clf;
    %hold on;
    %tsh = tsurf(F,U,'FaceColor',blue,falpha(0.1,0),fsoft);
    %ssh =     sct(U(end,:),'o','MarkerFaceColor','r');
    %ssh =     sct(U(t,:),'o','MarkerFaceColor','b');
    %ssh =     sct(U(g,:),'o','MarkerFaceColor','g');
    %%ssh = sct(U(g,:),'o','MarkerFaceColor','r');
    %%qvr(U(t,:),100*arm,'Color','r','LineWidth',3);
    %%qvr(U(t,:),100*up,'Color','g','LineWidth',3);
    %%qvr(U(t,:),100*ax,'Color','b','LineWidth',3);
    %%hold off;
    %%l = light('Position',[0 1e-10 1],'Style','infinite');
    %%add_shadow(tsh,l,'Fade','none');
    %%sh = add_shadow([],l);
    %%sh{1}.MarkerFaceColor = ssh.MarkerFaceColor;
    %axis equal;

      end

      U = (U-U(t(1),:))*R3+U(t(1),:);
      t = [t g];
      assert(all(isreal(U(:))));
    case 3
      if isempty(t1)
        t1 = sort(t);
        if (cross(U(t1(2),:)-U(t1(1),:),U(t1(3),:)-U(t1(1),:))*[0;0;1])>0
          t1 = t1([1 3 2]);
        end
      end
      cen = U(end,:);
      % check if stable
      H = convhull(U(t,1:2));
      stable = inpolygon(cen(1),cen(2),U(t(H),1),U(t(H),2));
      if ~stable
        % which points to remove
        E = t([H(1:end-1) H(2:end)]);
        debug = false;

        if debug
          clf;
          hold on;
          tsh = tsurf(F,U,'FaceColor',blue,falpha(0.1,0),fsoft);
          ssh = sct(cen,'o','MarkerFaceColor','r');
          ssh = sct(U(t,:),'o','MarkerFaceColor','b');
          hold off;
          view(2);
        end

        [sqrD,I,C] = point_mesh_squared_distance(cen(1:2),U(:,1:2),E);
        B = normrow(C-U(E(I,2),1:2))/normrow(U(E(I,1),1:2)-U(E(I,2),1:2));
        B = [B 1-B];
        tI = E(I,:);

        if debug
          hold on;
          plot_edges(U,tI,'LineWidth',3);
          txt(U(t,:),num2str(t'))
          hold off;
        end

        if B(1)<1e-8
          t = tI(2);
        elseif B(2)<1e-8
          t = tI(1);
        else
          t = tI;
        end
        if debug
          hold on;
          ssh = sct(U(t,:),'o','MarkerFaceColor','g');
          hold off;
          pause
        end
      end
      assert(all(isreal(U(:))));
    end
    
    
    
    tsh = [];
    %clf;
    %hold on;
    %tsh = tsurf(F,U,'FaceColor',blue,falpha(1.0,0),fsoft);
    %ssh = sct(U(end,:),'o','MarkerFaceColor','r');
    %ssh = sct(U(t,:),'o','MarkerFaceColor','g');
    %%ssh = sct(U(g,:),'o','MarkerFaceColor','r');
    %%qvr(U(t,:),100*arm,'Color','r','LineWidth',3);
    %%qvr(U(t,:),100*up,'Color','g','LineWidth',3);
    %%qvr(U(t,:),100*ax,'Color','b','LineWidth',3);
    %%hold off;
    %%view(2);
    %view(3);
    %l = light('Position',[0 1e-10 1],'Style','infinite');
    %add_shadow(tsh,l,'Fade','none');
    %%sh = add_shadow([],l);
    %%sh{1}.MarkerFaceColor = ssh.MarkerFaceColor;
    %axis equal;
    %camproj('persp');
    %drawnow;
    if stable
      t = sort(t);
      if (cross(U(t(2),:)-U(t(1),:),U(t(3),:)-U(t(1),:))*[0;0;1])>0
        t = t([1 3 2]);
      end
    end
  end
  % Strip off centroid
  U = U(1:end-1,:);

end

function [S,Sall] =  curve_smooth(S,E,varargin)
% Conformal mean curvature flow for curves
%
% [S,Sall] =  curve_smooth(S,E)
% [S,Sall] =  curve_smooth(S,E,'ParameterName',ParameterValue,...);
%
% Input:
%   S  #S by dim list of curve vertices
%   E  #E by 2 list of edge indices into S
%   Optional
%     'Lambda' followed by flow time step {0.1}
%     'MaxIters'  followed by maximum number of iters {2}
%     'FixBoundary'  whether to fix boundary {false}
% Output:
%   S  #S by dim list of final curve vertices
%   Sall #S by dim by iters list of all flow positions
%

% default values
lambda_flow = 0.1;
iters = 2;
fix_boundary = false;
% Map of parameter names to variable names
params_to_variables = containers.Map( ...
  {'Lambda','MaxIters','FixBoundary'}, ...
  {'lambda_flow','iters','fix_boundary'});
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

Sall(:,:,1) = S;
A = adjacency_matrix(E);
L = A-diag(sum(A,2));
if fix_boundary
  b = find(sum(A,1)==1);
  L(b,:) = 0;
end
%figure;
for ii=2:iters
 
    %plot([S(E(:,1),1) S(E(:,2),1)]',[S(E(:,1),2) S(E(:,2),2)]','-')
    %drawnow;
    
    A = adjacency_edge_cost_matrix(Sall(:,:,ii-1),E);
    M = diag(sum(A,2)/2);
    M = M./max(diag(M));
    % cMCF
    Sall(:,:,ii) = (M-lambda_flow*L)\(M*Sall(:,:,ii-1));
    %% rescale
    %c = mean(Sall(:,:,ii));
    %bbd = max(Sall(:,:,ii))-min(Sall(:,:,ii));
    %Sall(:,:,ii) = bsxfun(@plus,c,bsxfun(@minus,Sall(:,:,ii),c)* ...
    %  sqrt(sum((max(S)-min(S)).^2,2))/sqrt(sum(bbd.^2,2)));
    S = Sall(:,:,ii);
    
end

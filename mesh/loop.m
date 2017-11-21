function [VV,FF,SS,J] = loop(V,F,iter)
% LOOP perform loop subdivision. After n iterations of loop subivision, the
% resulting mesh will have
%   4^n |F| faces
%   2^(n-1)*(2*|E| + 3*|F|*(2^n-1)) edges
%   |V| + (2^n-1)|E| + (1+2^(n-1)*(2^n-3))*|F| vertices
%
% [VV,FF] = loop(V,F)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of face indices
% Outpus:
%   VV #VV by 3 new vertex positions
%   FF #FF by 3 list of face indices
%   SS #VV by #V matrix computing VV = SS *V 
%   J  #FF list of indices into F
%
% Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
%

if (~exist('iter','var'))
    iter = 1;
end
VV = V;
SS = speye(size(V,1));
FF = F;
J = (1:size(F,1))';

for i=1:iter
    % number of original vertices
    n = size(VV,1);
    % number of original faces
    nf = size(FF,1);

    % outline edges
    O = outline(FF);

    % indices of original vertices
    original = 1:n;
    % extract unique vertex indices on outline
    [boundary,~,~] = unique(O(:));
    % original interior vertices
    interior = original(~ismember(original,boundary));
    % number of original
    no = numel(original);
    % number of interior original
    ni = numel(interior);
    % number of boundary original
    nb = numel(boundary);

    % get adjacency matrix
    A = adjacency_matrix(FF);
    % number of edges
    ne = sum(A(:))/2;

    % valence values
    val = sum(A,2);
    beta = (1./val) .* (5/8 - ( 3/8 + (1/4)*(cos ( 2*pi./val ) )).^2);
    % subdivision matrix for interior (even) vertices
    Seven = sparse(n,n);

    %% as if all vertices are interior, then splice out interior
    %Sint = bsxfun(@times,A,beta) + diag(sparse(1-val.*beta));
    %S(interior,:) = Sint(interior,:);
    % construct interior coefficients directly
    Seven(interior,:) = ...
        bsxfun(@times,A(interior,:),beta(interior)) + ...
        sparse(1:ni,original(interior),1-val(interior).*beta(interior),ni,no);

    % construct boundary coefficients
    Sboundary = ...
        sparse([O(:,1);O(:,2)],[O(:,2);O(:,1)],1/8,n,n) + ...
        sparse([O(:,1);O(:,2)],[O(:,1);O(:,2)],3/8,n,n);
    Seven(boundary,:) = Sboundary(boundary,:);

    % odd vertex subdivision matrix

    % add a new vertex for each edge

    % get "face edges" aka half-edges
    FE = sort([FF(:,1) FF(:,2);FF(:,2) FF(:,3);FF(:,3) FF(:,1)],2);
    % get unique edge list, FE2E tells original face edges (with duplicates, EE)
    % where to find matching unique edge (E)
    [E,I,FE2E] = unique(FE,'rows');

    % #E by #V matrix, flaps(i,j) = 1 only if vertex j is a "flap" vertex of
    % edge i, i.e. vertex j *shares a face* both E(i,1) and E(i,2)
    FEflaps = [FF(:,3);FF(:,1);FF(:,2)];
    flaps = sparse(FE2E,FEflaps,1);


    % new vertices on boundary
    onboundary = sparse(sum(flaps,2) == 1);

    % don't consider flaps for boundary vertices
    flaps(onboundary,:) = 0;

    % THIS IS WRONG:
    %% flag for extra ordinary vertices
    %extraordinary = sum(A,2) > 7;
    %% boundary vertices are by default extraordinary
    %extraordinary(E(:,1)) = extraordinary(E(:,1)) | onboundary;
    %extraordinary(E(:,2)) = extraordinary(E(:,2)) | onboundary;
    %% both edge endpoints are extraordinary, for this case will ignore that
    %% they're extra ordinary
    %both = extraordinary(E(:,1)) & extraordinary(E(:,2));
    %neither = (~extraordinary(E(:,1))) & (~extraordinary(E(:,2)));
    %neither = neither | both;
    %Sodd = ...
    %    sparse( ...
    %    [1:size(E,1) 1:size(E,1)]', ...
    %    [E(:,1);E(:,2)], ...
    %    ([onboundary;onboundary]).*1/2 + ...
    %    (~[onboundary;onboundary]).*( ...
    %    [neither;neither] .* 3/8 + ...
    %    (~[neither;neither]) .* ( ...
    %    [extraordinary(E(:,1));extraordinary(E(:,2))].*1/2 + ...
    %    (~[extraordinary(E(:,1));extraordinary(E(:,2))]).*1/4)), ...
    %    ne, ...
    %    n) + ...
    %    flaps*1/8;

    % Without considering extraordinary vertices
    Sodd = ...
      sparse( ...
        [1:size(E,1) 1:size(E,1)]', ...
        [E(:,1);E(:,2)], ...
        ([onboundary;onboundary]).*1/2 + (~[onboundary;onboundary]).*3/8, ...
        ne, ...
        n) + ...
      flaps*1/8;

    % there are only one or two flap vertices per edge
    assert(all(sum(flaps,2) <= 2))
    % rows should sum to one
    assert(all(sum(Sodd,2) == 1))

    % indices of new points as if we really added a new point for each half edge
    i3 = n     + (1:nf)';
    i1 = n+nf  + (1:nf)';
    i2 = n+2*nf+ (1:nf)';
    % new face indices, 4 new faces for each original face. As if we simply
    % ignored the duplicates in m and had appended m to V
    FEF = [ FF(:,1) i3 i2 ; FF(:,2) i1 i3 ; FF(:,3) i2 i1 ; i1 i2 i3];
    J = [J;J;J;J];
    % reindex map from duplicate midpoint indices to unique midpoint indices
    FE2E = [(1:n)';FE2E+n];
    % reindex faces
    FF = FE2E(FEF);

    S = [Seven;Sodd];
    VV = S*VV;
    SS = S*SS;

    %trisurf(F,V(:,1),V(:,2),V(:,3),'FaceAlpha',0.1,'FaceColor','r','EdgeColor',[0.3 0 0]);
    %hold on;
    %trisurf(F,VVeven(:,1),VVeven(:,2),VVeven(:,3),'FaceAlpha',0.1,'FaceColor','b','EdgeColor',[0 0 0.3]);
    %trisurf(FF,VV(:,1),VV(:,2),VV(:,3),'FaceAlpha',0.1,'FaceColor','g','EdgeColor',[0 0.3 0.0]);
    %plot3(VVodd(:,1),VVodd(:,2),VVodd(:,3),'y.');
    %hold off;
    %view(2);

    %tsurf(F,V);
    %text(V(:,1),V(:,2),V(:,3),num2str(S(694,:)'),'BackgroundColor',[.8 .8 .8]);
end

end

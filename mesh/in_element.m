function [I,B1,B2,B3] = in_element(V,F,P,varargin)
  % IN_ELEMENT test for each p in P whether it lies inside each f in F defined
  % over V.
  % 
  % I = in_element(V,F,P)
  % [I,B1,B2,B3] = in_element(V,F,P,'ParameterName',ParameterValue,...)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim+1 list of element indices
  %   P  #P by dim list of query positions
  %   Optional:
  %     'Method' followed by one of the following {'knn'}:
  %       'brute-force' no acceleration O(#P * #F)
  %       'edge-walk' walk along edges ~O(#P * sqrt(#F)) Starting with a random
  %         barycenter step along the edges that intersect with the ray toward
  %         the query point. If a boundary is reached then search along all
  %         boundary edges and jump to farthest hit. **dim=2 only**
  %       'knn' use knnsearch to find closest element barycenters
  %       'spatial-hash' spatial hash on regular grid ~O(#P * sqrt(#F)) **dim=2
  %         only**
  %     'First'  only keep first match {false}
  %     'Quiet' suppress warnings {false}
  %     'Epsilon' epsilon used for determining inclusion {eps}
  % Outputs:
  %   I  #P by #F matrix of bools
  %   B1  #P by #F list of barycentric coordinates
  %   B2  #P by #F list of barycentric coordinates
  %   B3  #P by #F list of barycentric coordinates 
  %
  % Example:
  %   P = bsxfun(@plus,min(V),bsxfun(@times,rand(100,2),max(V)-min(V)));
  %   [I,B1,B2,B3] = in_element(V,F,P);
  %   % Only keep first
  %   [mI,J] = max(I,[],2);
  %   I = sparse(1:size(I,1),J,mI,size(I,1),size(I,2));
  %   % Mask barycentric coordinates
  %   B1 = B1.*I;
  %   B2 = B2.*I;
  %   B3 = B3.*I;
  %   Q = B1*V(F(:,1),:) + B2*V(F(:,2),:) + B3*V(F(:,3),:);
  %   Q = Q(any(I,2),:);
  %   tsurf(F,V);
  %   hold on;
  %   plot(Q(:,1),Q(:,2),'or','LineWidth',6);
  %   plot(P(:,1),P(:,2),'ob','LineWidth',2);
  %   hold off;
  %


  function [I] = in_element_brute_force(V,F,P)
    dim = size(V,2);
    assert(dim+1 == size(F,2));
  
    % number of elements 
    m = size(F,1);
    % number of query points 
    np = size(P,1);
    
    switch dim 
    case 3
      T = F;
      % tet face ares
      vol = abs(volume(V,T));
      allF = [ ...
        T(:,2) T(:,4) T(:,3); ...
        T(:,1) T(:,3) T(:,4); ...
        T(:,1) T(:,4) T(:,2); ...
        T(:,1) T(:,2) T(:,3); ...
        ];
      % List of tets for each face f of each tet t for each point p
      TP = cat(2, ...
        repmat(allF,[1 1 np]), ...
        permute(repmat(size(V,1)+(1:np),m*4,1),[1 3 2]));
      TP = reshape(permute(TP,[1 3 2]),m*4*np,dim+1);
      Pvol = abs(volume([V;P],TP));
      % Pvol(t,f,p) --> volume of tet t, face f with point p
      Pvol = reshape(Pvol,[size(T,1) dim+1 np]);
      % sumvol(p,t) --> sum of volumes of tets made with faces of t and p
      sumvol = permute(sum(Pvol,2),[3 1 2]);
      I = sparse(abs(bsxfun(@minus,sumvol,vol')) < sqrt(epsilon));
    case 2
      % triangle side lengths
      l = [ ...
        sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
        sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
        sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
        ];
  
      B = zeros([np m dim+1]);
      for ii = 1:(dim+1)
        jj = mod(ii+1,dim+1)+1;
        kk = mod(ii,dim+1)+1;
        ljj = pdist2(P,V(F(:,jj),:));
        lkk = pdist2(P,V(F(:,kk),:));
        
        % semiperimeters
        s = bsxfun(@plus,l(:,ii)',ljj + lkk)*0.5;
        % Heron's formula for area
        B(:,:,ii) = 2*sqrt(s.*(bsxfun(@minus,s,l(:,ii)').*(s-ljj).*(s-lkk)));
      end
      % sum of barycentric coordinates
      sumA = sum(B,3);
      % area of element
      dblA = doublearea(V,F);
      %% check whether sum is more than true are
      %I = ~bsxfun(@gt,sumA,dblA');
      I = sparse((bsxfun(@minus,sumA,dblA')) < sqrt(epsilon));
    end
    %B1 = sparse(B(:,:,1));
    %B2 = sparse(B(:,:,2));
    %B3 = sparse(B(:,:,3));
  end

  function I = in_element_hash_helper(V,F,P)
    assert(size(F,2) == 3, ...
      'F must contain triangles for Method=''spatial-hash''');
    num_bins = ceil(sqrt(size(F,1)));
    bin_x = ceil(sqrt(num_bins));
    bin_y = ceil(num_bins/bin_x); 
    num_bins = bin_x*bin_y;
    % spatial hash
    function VH = hash(V,MN,MX,bin_x,bin_y)
      [~,X] = histc(V(:,1),linspace(MN(1),MX(1),bin_x));
      [~,Y] = histc(V(:,2),linspace(MN(2),MX(2),bin_y));
      VH = sub2ind([bin_x bin_y],X,Y);
    end
    %% http://stackoverflow.com/a/5929567/148668
    %primes = [ 40960001, 59969537 45212177];
    %hash = @(X) mod( ...
    %  bitxor( int32(V(:,1)*primes(1)), int32(V(:,2)*primes(2))), ...
    %  num_bins);
    MN = min([V;P]);
    MX = max([V;P]);
    VH = hash(V,MN,MX,bin_x,bin_y);
    PH = hash(P,MN,MX,bin_x,bin_y);
    % This is wrong for triangles that span more hash cells than their vertices:
    % any time a triangle lands on the corner....
    FH = sparse(repmat(1:size(F,1),1,size(F,2))',VH(F(:)),1,size(F,1),num_bins)~=0;
    [~,FHI,FHJ] = find(FH);
    [FHJX,FHJY] = ind2sub([bin_x,bin_y],FHJ);
    % This assumes that #P >> #F
    I = sparse(size(P,1),size(F,1));
    for h = 1:num_bins
      %Vh = V(VH==h,:);
      IFh = FH(:,h);
      if any(IFh)
        IPh = PH==h;
        Ph = P(IPh,:);
        if any(IPh)
          Fh = F(IFh,:);
          Ih = in_element_brute_force(V,Fh,Ph);
          I(IPh,IFh) = Ih;
        end
      end
    end
  end

  % default values
  method = 'knn';
  first = false;
  quiet = false;
  epsilon = eps;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Method','First','Quiet','Epsilon'}, ...
    {'method','first','quiet','epsilon'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  switch method
  case 'brute-force'
    I = in_element_brute_force(V,F,P);
  case 'spatial-hash'
    % Try 45?? spatial grid, too (reduce corner cases)
    switch size(V,2)
    case 2
      R = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
    case 3
      R = [cos(pi/4) -sin(pi/4) 0 ;sin(pi/4) cos(pi/4) 0; 0 0 1];
    end
    I = in_element_hash_helper(V,F,P) | in_element_hash_helper(V*R,F,P*R);

    % Find any obviously incorrect values: not inside but winding number says
    % inside. Could still missing something if mesh overlaps itself.
    NI = ~any(I,2);
    switch size(F,2)
    case 3
      O = outline(F);
    case 4
      O = boundary_faces(F);
    end
    WI = abs(winding_number(V,O,P(NI,:))/(2*pi))>0.5;
    WI = sparse(find(NI),1,WI,size(P,1),1)~=0;
    % redo any that currently are not in any but winding number says are inside
    RI = NI & WI;
    I(RI,:) = in_element_brute_force(V,F,P(RI,:));

  case 'edge-walk'

    assert(size(F,2) == 3,'F must contain triangles for Method=''edge-walk''');
    % List of all "half"-edges: 3*#F by 2
    allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
    % Sort each row
    sortallE = sort(allE,2);
    % IC(i) tells us where to find sortallE(i,:) in uE:
    % so that sortallE(i,:) = uE(IC(i),:)
    [uE,~,IC] = unique(sortallE,'rows');
    % uE2F(e,f) = i means face f's ith edge is unique edge e
    uE2F = sparse(IC(:),repmat(1:size(F,1),1,3)',reshape(repmat(1:3,size(F,1),1),[],1));
    % uE2F(e,f) = 1 means face f is adjacent to unique edge e
    uE2F1 = sparse(IC(:),repmat(1:size(F,1),1,3)',1);
    % Face-face Adjacency matrix
    A = uE2F1'*uE2F;
    % A(f,g) = i means face f's ith edge is shared with g
    A = A-diag(diag(A));

    I = sparse(size(P,1),size(F,1));
    B1 = sparse(size(I,1),size(I,2));
    B2 = sparse(size(I,1),size(I,2));
    B3 = sparse(size(I,1),size(I,2));

    [~,is_b] = on_boundary(F);
    EF = repmat(1:size(F,1),1,3)';
    EFI = reshape(repmat(1:3,size(F,1),1),[],1);
    O = allE(is_b(:),:);
    OF = EF(is_b(:));
    OFI = EFI(is_b(:));

    % initial closest vertex
    BC = barycenter(V,F);
    % Centroid
    f_init = snap_points(mean(V),BC);
    %% Random initial guess
    %f_init = ceil(rand(1,1)*size(F,1));
    for p = 1:size(P,1)
      % current closest face, barycenter
      f = f_init;
      q = BC(f_init,:);
      % incoming edge
      e_in = [];
      while true
        % current point
        % edges to test
        E = ... %mod(bsxfun(@plus,setdiff([1;2;3],e_in),-1+(1:2)),3)+1;
          [2 3;3 1;1 2];
        out = ...
          lineSegmentIntersect([q P(p,:)],[V(F(f,E(:,1)),:) V(F(f,E(:,2)),:)]);
        out.intAdjacencyMatrix(e_in) = false;
        e_out = find(out.intAdjacencyMatrix);
        if isempty(e_out)
          % no hits so we're in the element
          I(p,f) = 1;
          B = barycentric_coordinates( ...
            P(p,:),V(F(f,1),:),V(F(f,2),:),V(F(f,3),:));
          B1(p,f) = B(1); B2(p,f) = B(2); B3(p,f) = B(3);
        else
          f_prev = f;
          if numel(e_out) > 1
            % Take farthests
            [sd,si] = sort(out.intNormalizedDistance1To2(e_out),'descend');
            e_out = e_out(si(1));
            % TODO: recurse on all in f
          end
          f = find(A(:,f_prev)==e_out);
          % Todo if 
          if isempty(f)
            % no neighbors so we're shooting outside.
            out = ...
              lineSegmentIntersect([q P(p,:)],[V(O(:,1),:) V(O(:,2),:)]);
            hits = find(out.intAdjacencyMatrix);
            [sd,si] = sort(out.intNormalizedDistance1To2(hits),'descend');
            f = OF(hits(si(1)));
            e_in = OFI(hits(si(1)));
            if f==f_prev
              error('Point is outside');
            end
          else
            if numel(f) > 1
              f = f(1);
              if ~quiet
                warning('Ignoring non-manifold edge: might miss multiply inside.');
              end
              % TODO: recurse on all in f
            end
            e_in = A(f_prev,f);
          end
        end
        %tsurf(F,V);
        %hold on;
        %tsurf(F(f,:),V,'CData',-1);
        %tsurf(F(f_prev,:),V,'CData',1);
        %plot_edges([q;P(p,:)],[1 2],'r');
        %plot_edges(V,[F(f_prev,E(:,1));F(f_prev,E(:,2))]','b');
        %plot(P(p,1),P(p,2),'*');
        %plot_edges(V,[F(f_prev,E(e_out,1));F(f_prev,E(e_out,2))]', ...
        %  'y','LineWidth',2);
        %hold off;
        %drawnow;
        %input('');
      end
    end

    return
  case 'knn'
    % assumes samples are all inside exactly one element
    BC = barycenter(V,F);
    I = sparse(size(P,1),size(F,1));
    B1 = sparse(size(P,1),size(F,1));
    B2 = sparse(size(P,1),size(F,1));
    B3 = sparse(size(P,1),size(F,1));
    B4 = sparse(size(P,1),size(F,1));
    % indices of points we haven't found yet
    IP = 1:size(P,1);

    prev_k = 0;
    k = 2;
    while true
      K = knnsearch(BC,P(IP,:),'K',k);
      K = K(:,prev_k+1:end);
      for ki = 1:size(K,2)
        switch size(F,2)
        case 3
          B = abs(barycentric_coordinates( ...
            P(IP,:), ...
            V(F(K(:,ki),1),:), ...
            V(F(K(:,ki),2),:), ...
            V(F(K(:,ki),3),:)));
        case 4
          B = abs(barycentric_coordinates( ...
            P(IP,:), ...
            V(F(K(:,ki),1),:), ...
            V(F(K(:,ki),2),:), ...
            V(F(K(:,ki),3),:), ...
            V(F(K(:,ki),4),:)));
        end
        found = abs(sum(B,2)-1)<sqrt(epsilon);
        I(sub2ind(size(I),IP,K(:,ki)')) = found;
        % DOESN'T REALLY PAY OFF TO COMPUTE THESE HERE
        %B1(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,1);
        %B2(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,2);
        %B3(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,3);
        %if size(F,2) == 4
        %  B4(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,4);
        %end
        % Peel off found
        IP = IP(~found);
        if isempty(IP)
          break;
        end
        K = K(~found,:);
      end

      %for p = 1:numel(IP)
      %  Kp = K(p,:);
      %  U = V(F(Kp,:),:);
      %  FKp = reshape((1:size(U,1))',size(U,1)/size(F,2),size(F,2));
      %  I(IP(p),Kp) = in_element_brute_force(U,FKp,P(IP(p),:));
      %end
      %IP = find(~any(I,2));

      if isempty(IP)
        break;
      end
      prev_k = k;
      k = min(prev_k*2,size(BC,1));
      if k == size(BC,1)
        if ~quiet
          warning('Some points not found');
        end
        break;
      end
    end

    %return;

  end

  if first
    % Only keep first
    [mI,J] = max(I,[],2);
    I = sparse(1:size(I,1),J,mI,size(I,1),size(I,2));
  end

  % Compute barycentric coordinates
  if nargout > 1
    B1 = sparse(size(I,1),size(I,2));
    B2 = sparse(size(I,1),size(I,2));
    B3 = sparse(size(I,1),size(I,2));
    B4 = sparse(size(I,1),size(I,2));
    work_I = I;
    while true
      % Peel off a layer
      [mIf,Jf] = max(work_I,[],2);
      If = find(mIf);
      Jf = Jf(If);
      mIf = mIf(If);
      if isempty(If)
        break
      end
      work_I = work_I - sparse(If,Jf,mIf,size(work_I,1),size(work_I,2));;
      Pf = P(If,:);
      Ff = F(Jf,:);
      switch size(Ff,2)
      case 4
        Bf = barycentric_coordinates(Pf, ...
          V(Ff(:,1),:),V(Ff(:,2),:),V(Ff(:,3),:),V(Ff(:,4),:));
      case 3
        Bf = barycentric_coordinates(Pf,V(Ff(:,1),:),V(Ff(:,2),:),V(Ff(:,3),:));
      end
      B1(sub2ind(size(B1),If,Jf)) = Bf(:,1);
      B2(sub2ind(size(B2),If,Jf)) = Bf(:,2);
      B3(sub2ind(size(B3),If,Jf)) = Bf(:,3);
      if size(Bf,2) >= 4
        B4(sub2ind(size(B3),If,Jf)) = Bf(:,4);
      end
    end
  end

end



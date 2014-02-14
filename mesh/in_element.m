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
  %     'BruteForce' followed by whether to use brute force computation without
  %       acceleration.
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
      ljj = all_pairs_distances(P,V(F(:,jj),:));
      lkk = all_pairs_distances(P,V(F(:,kk),:));
      
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
    I = sparse((bsxfun(@minus,sumA,dblA')) < sqrt(eps));
    %B1 = sparse(B(:,:,1));
    %B2 = sparse(B(:,:,2));
    %B3 = sparse(B(:,:,3));
  end

  function I = in_element_hash_helper(V,F,P)
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
    FH = sparse(repmat(1:size(F,1),1,3)',VH(F(:)),1,size(F,1),num_bins)~=0;
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
  brute_force = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'BruteForce'}, {'brute_force'});
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

  if brute_force
    I = in_element_brute_force(V,F,P);
  else
    % Try 45° spatial grid, too (reduce corner cases)
    R = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
    I = in_element_hash_helper(V,F,P) | in_element_hash_helper(V*R,F,P*R);

    % Find any obviously incorrect values: not inside but winding number says
    % inside. Could still missing something if mesh overlaps itself.
    NI = ~any(I,2);
    O = outline(F);
    WI = abs(winding_number(V,O,P(NI,:))/(2*pi))>0.5;
    WI = sparse(find(NI),1,WI,size(P,1),1)~=0;
    % redo any that currently are not in any but winding number says are inside
    RI = NI & WI;
    I(RI,:) = in_element_brute_force(V,F,P(RI,:));

    B1 = sparse(size(I,1),size(I,2));
    B2 = sparse(size(I,1),size(I,2));
    B3 = sparse(size(I,1),size(I,2));

  end

  work_I = I;
  while true
    % Peal off a layer
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
    Bf = barycentric_coordinates(Pf, V(Ff(:,1),:),V(Ff(:,2),:),V(Ff(:,3),:));
    B1(sub2ind(size(B1),If,Jf)) = Bf(:,1);
    B2(sub2ind(size(B2),If,Jf)) = Bf(:,2);
    B3(sub2ind(size(B3),If,Jf)) = Bf(:,3);
  end

end



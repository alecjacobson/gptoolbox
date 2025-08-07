function U = jiggle_duplicates(V,F,epsilon,sigma)
  function flag = selfintersect_free(V,F)
    [~,~,IF] = selfintersect(V,F,'DetectOnly',true,'FirstOnly',true);
    flag = isempty(IF);
  end
  assert(selfintersect_free(V,F));
  %assert(epsilon > sigma);

  function I = dups(U,epsilon)
    if epsilon > 0
      [uU,K,J] = unique(round(U/epsilon),'rows');
    else
      [uU,K,J] = unique(U,'rows');
    end
    C = accumarray(J,1);
    I = K(C > 1);
  end

  U = V;
  max_iter = 100;
  A = adjacency_matrix(F);
  for iter = 1:max_iter
    assert(selfintersect_free(U,F));
    I = dups(U,epsilon);
    if isempty(I)
      return;
    end
    W = U;
    W(I,:) = W(I,:) + sigma*randn(length(I),size(W,2));

    if selfintersect_free(W,F)
      U = W;
    else
      % OK so there were interesections
      [~,~,IF] = selfintersect(W,F,'DetectOnly',true);
      % do all involve a duplicate?
      involve_I = any(ismember(F(IF(:,1),:),I),2) | any(ismember(F(IF(:,2),:),I),2);
      if all(involve_I)
        % then I claim we can safely move the entries in I that aren't involed
        % in IF
        safe = ~ismember(I,unique(F(IF(:),:)));
        U(I(safe),:) = W(I(safe),:);
      end
    end
  end
  warning('jiggle_duplicates: max_iter reached');
  U = [];

end

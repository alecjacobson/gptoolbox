function [PP,CC,II,J] = spline_split(P,C,I,T)
  % [PP,CC,II,J] = spline_split(P,C,I,T)
  %

  IT = [I T];

  if size(IT,1) == 0
    [PP,CC,II,J] = deal(P,C,[],[]);
    return;
  end

  [IT,K] = sortrows(IT);
  I = IT(:,1);
  T = IT(:,2);


  CC = nan(size(C,1)+numel(I),size(C,2));
  PP = [P;nan(numel(I)*3,size(P,2))];
  % cursor into PP
  pc = size(P,1);

  c = 0;
  cc = 0;
  j = 0;
  II = nan(size(C,1)+numel(I),1);
  J = nan(numel(I),1);
  for ii = 1:numel(I)
    i = I(ii);
    if i > j
      untouched = (j+1):(i-1);
      %fprintf('   untouched: '); fprintf('%d ',untouched); fprintf('\n');
      CC(cc+(1:numel(untouched)),:) = C(untouched,:);
      II(cc+(1:numel(untouched))) = untouched;
      cc = cc + numel(untouched);
      % first time we're seeing a split from i
      j = i;
      cc = cc+1;
      CC(cc,:) = C(j,:);
      II(cc) = j;
      chewed = 0;
    end

    t = (T(ii)-chewed)/(1-chewed);
    chewed = T(ii);

    % CC(cc,1)---CC(cc,2)---CC(cc,3)---CC(cc,4)
    [c1,c2] = cubic_split(PP(CC(cc,:),:),t);
    % CC(cc,1)=c1(1)----CC(cc,2)=c1(2)----CC(cc,3)=c1(3)----CC(cc,4)=c1(4)=c2(1)
    % CC(cc+1,1)=c2(1)----CC(cc+1,2)=c2(2)----CC(cc+1,3)=c2(3)----CC(cc+1,4)=c2(4)
    new_p = pc+(1:3);
    PP(new_p,:) = [c1(3,:); c2(1,:); c2(2,:)];
    pc = pc + 3;
    % modify tangents
    PP(CC(cc,2),:) = c1(2,:);
    PP(CC(cc,3),:) = c2(3,:);
    CC(cc+1,:) = [new_p(2) new_p(3) CC(cc,3) CC(cc,4)];
    II(cc+1) = j;
    CC(cc,3:4) = [new_p(1) new_p(2)];
    II(cc) = j;
    J(ii) = cc;
    cc = cc+1;
  end
  % flush remaining
  if j < size(C,1)
    untouched = (j+1):size(C,1);
    %fprintf('   untouched: '); fprintf('%d ',untouched); fprintf('\n');
    CC(cc+(1:numel(untouched)),:) = C(untouched,:);
    II(cc+(1:numel(untouched))) = untouched;
    cc = cc + numel(untouched);
  end
  J(K) = J;




end

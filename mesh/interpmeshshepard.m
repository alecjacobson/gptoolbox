function [ d ] = interpmeshshepard( V,F,ind,val)

%% Compute biharmonic embedding
%B = biharmonic_embedding_yaron(V,F);

%% Compute weights (1/d with d biharmonic distances)
w = zeros(size(V,1),size(ind,2));

for i=1:size(ind,2)
    %w(:,i) = normrow(repmat(B(ind(i),:),size(V,1),1) - B);
    w(:,i) = geodesicdistance(V,F,[ind(i)]);
end

w = ones(size(V,1),size(ind,2))./w;

%% Shepard's interpolation
d = zeros(size(V,1),1);
for i=1:size(V,1)
    d(i) = (w(i,:)*val')./sum(w(i,:));
end

d(ind) = val;

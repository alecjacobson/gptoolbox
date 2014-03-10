function [id,d] = df_query( DF, p)
% DF_QUERY Use the spatial index build with df_build to provide approximate NN
% and distance from the NN.
% 
% [id,d] = df_query(DF,p)
%
% Input:
%   DF  output of df_build
%   p  query point
% Output:
%   id  index of the point closer to p (approximate)
%   d  distance from the point closer to p (approximate)

p = p-DF.MIN;
p = p./DF.S;
p = round(p)+1;

p(p<1) = 1;
p(p>DF.C') = DF.C(p>DF.C');

id = DF.N(p(1),p(2),p(3));
d  = DF.D(p(1),p(2),p(3));

end


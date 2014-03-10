function ET = et(V,F,E)
%  ET Compute the relation edge to triangle of a manifold mesh.
%
%  ET = et(V,F,E)
%
%  Input:
%   V  #V x 3  matrix of vertex coordinates
%   F  #F x 3  matrix of indices of vertices at triangle corners
%   E  #E x 2  matrix of indices of vertices at edge endpoints
%  Output:
%   ET #E x 2  map between an edge to its two indicent faces
%              (-1 in column 2 if the edge is on the border)
%
%  Example:
%  E = edges(V,F);
%  ET = et(V,F,E);

VT = vt(V,F);

ET = zeros(size(E,1),2);

for i=1:size(E,1)
    CFi1 = find(VT(:,E(i,1)))';
    CFi2 = find(VT(:,E(i,2)))';
    
    temp = intersect(CFi1,CFi2);
    if size(temp,2) == 1
        temp = [temp -1];
    end
    
    ET(i,:) = temp;
end

end


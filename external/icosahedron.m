function [vMat,fMat] = icosahedron()
% ICOSAHEDRON generates a triangulation corresponding to Icosahedron's
% surface
%
% Input:
%   none
%
% Output:
%   vMat: double[12,3] - vertices of the Icosahedron
%   fMat: double[20,3] - indices matrix of face definitions
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
IND_VEC=transpose(0:4);
Z_VEC=0.5*ones(5,1);
pi = 4 * atan(1.0);
tau = (realsqrt(5.0) + 1)/2;
r = tau - 0.5;
vMat([1,12],:)=[0.0, 0.0, 1.0;0.0, 0.0, -1.0];
%
alphaVec=-pi/5 + IND_VEC * pi/2.5;
vMat(2+IND_VEC,:)=[ cos(alphaVec)/r, sin(alphaVec)/r, Z_VEC/r];
%
alphaVec = IND_VEC * pi/2.5;
vMat(7+IND_VEC,:)=[cos(alphaVec)/r, sin(alphaVec)/r, -Z_VEC/r];
vMat(12,:)=[ 0.0, 0.0, -1.0];
fMat=[1 2 3;1 3 4;1 4 5;1 5 6;1 6 2;2 7 3;3 8 4;4 9 5;5 10 6;6 11 2;...
    7 8 3;8 9 4;9 10 5;10 11 6;11 7 2;7 12 8;...
    8 12 9;9 12 10;10 12 11;11 12 7];

function [V,F] = cat_meshes(varargin)
% CAT_MESHES concatenate many meshes
%
% [V,F] = cat_meshes(V1,F1,V2,F2, ...
%
% [V,F] = cat_meshes(fvStruct1, ...)
%
% Inputs:
%   V1  #V by dim list of vertex positions of mesh 1
%   F1  #F by simplex-size list of simplex indices of mesh 1
%   V2  #V by dim list of vertex positions of mesh 2
%   F2  #F by simplex-size list of simplex indices of mesh 2
%   ...
%   
%   fvStruct a struct (array) with the fields vertices and faces, e.g.:
%   fvStruct1(1).vertices = V1; fvStruct1(1).faces = F1;
%   fvStruct1(2).vertices = V2; fvStruct1(2).faces = F2;
%   fvStruct2(1).vertices = V3; fvStruct2(1).faces = F3;
%   ...
%   
% Outputs:
%   V  #V1+#V2+... by dim list of vertex positions
%   F  #F1+#F2+... by simplex-size list of simplex indices
%

assert(~isempty(varargin))

if isstruct(varargin{1})
    VF_fields = {'vertices','faces'};
    
    errorStructFields=['If the first input argument is a struct '...
        'with the fields vertices and faces the additonal ' ...
        'arguments must have the same format'];
    % Check, if all input arguments are structs
    assert(all(cellfun(@isstruct, varargin)), errorStructFields)
    % Check, if all structs contain the two fields vertices and faces
    assert(all(cellfun(@(x) all(ismember(fieldnames(x), ...
        VF_fields)), varargin)), errorStructFields)
    
    if length(varargin)==1
        errorArgAndStructLength = ['If the input is only one struct ' ...
            'it has to contain more than one mesh.'];
        assert(length(varargin{1})>1, ...
            errorArgAndStructLength)
    end
    
    % Order of the fields: vertices, faces
    varargin = cellfun(@(x) orderfields(x, VF_fields),varargin, 'UniformOutput',0);
    
    % Convert the structs into one cell array
    varargin = ...
        cellfun(@struct2cell, varargin, 'UniformOutput', false);
    varargin = cellfun(@squeeze, varargin, 'UniformOutput',0);
    varargin = reshape([varargin{:}],[],1)';
end

NoA = length(varargin);
assert(mod(NoA,2)==0);

cellfun(@(x) validateattributes(x, {'numeric'},...
    {'size',[NaN,3],'finite'}), varargin(1:2:end))
cellfun(@(x) validateattributes(x, {'numeric'},...
    {'integer'}), varargin(2:2:end))

V=[];
F=[];

for m = 1:NoA/2
    Vm = varargin{2*m-1};
    Fm = varargin{2*m};
    F = [F; Fm+size(V,1)];
    V = [V; Vm];
end

end


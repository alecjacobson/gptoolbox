function write_ply(vertex,face,filename, mode)

% write_ply - write data to PLY file.
%
%   write_ply(vertex,face,filename, mode);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   IMPORTANT: works only for triangular meshes.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<4
    mode = 'ascii';
end

if size(vertex,2)~=3
    vertex = vertex';
end
if size(vertex,2)~=3
    error('vertex does not have correct format.');
end


if size(face,2)~=3
    face = face';
end
if size(face,2)~=3
    error('face does not have correct format.');
end

% make a cube
clear Data;
Data.vertex.x = vertex(:,1);
Data.vertex.y = vertex(:,2);
Data.vertex.z = vertex(:,3);
Data.face.vertex_indices = {};
for i=1:size(face,1)
    Data.face.vertex_indices{end+1} = face(i,:)-1;
end
plywrite(Data,filename,mode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plywrite(Elements,Path,Format,Str)

%PLYWRITE  Write 3D data as a PLY file.
%   PLYWRITE(DATA,FILENAME) writes the structure DATA as a binary 
%   PLY file.  Every field of DATA is interpreted as an element
%   and every subfield as an element property.  Each subfield of
%   property data must either be an array or a cell array of 
%   arrays.  All property data in an element must have the same
%   length.
%
%   A common PLY data structure has the following fields:
%      DATA.vertex.x = x coordinates, [Nx1] real array
%      DATA.vertex.y = y coordinates, [Nx1] real array
%      DATA.vertex.z = z coordinates, [Nx1] real array
%
%      DATA.face.vertex_indices = vertex index lists, 
%         an {Mx1} cell array where each cell holds a one-
%         dimesional array (of any length) of vertex indices.
%   Some other common data fields:
%      DATA.vertex.nx = x coordinate of normal, [Nx1] real array
%      DATA.vertex.ny = y coordinate of normal, [Nx1] real array
%      DATA.vertex.nz = z coordinate of normal, [Nx1] real array
%
%      DATA.edge.vertex1 = index to a vertex, [Px1] integer array
%      DATA.edge.vertex2 = second vertex index, [Px1] integer array
%   Many other fields and properties can be added.  The PLY format 
%   is not limited to the naming in the examples above -- they are
%   only the conventional naming.
%
%   PLYWRITE(DATA,FILENAME,FORMAT) write the PLY with a specified 
%   data format, where FORMAT is
%      'ascii'                  ASCII text data
%      'binary_little_endian'   binary data, little endian
%      'binary_big_endian'      binary data, big endian (default)
%
%   PLYWRITE(DATA,FILENAME,FORMAT,'double') or
%   PLYWRITE(DATA,FILENAME,'double') write floating-point data as
%   double precision rather than in the default single precision.
%
%   Example:
%   % make a cube
%   clear Data;
%   Data.vertex.x = [0;0;0;0;1;1;1;1];
%   Data.vertex.y = [0;0;1;1;0;0;1;1];
%   Data.vertex.z = [0;1;1;0;0;1;1;0];
%   Data.face.vertex_indices = {[0,1,2,3],[7,6,5,4], ...
%         [0,4,5,1],[1,5,6,2],[2,6,7,3],[3,7,4,0]};
%   plywrite(Data,'cube.ply','ascii');
%
%   See also: PLYREAD

% Pascal Getreuer 2004

if nargin < 4
   Str = '';
   
   if nargin < 3
      Format = 'binary_big_endian';
   elseif strcmpi(Format,'double')
      Str = 'double';
      Format = 'binary_big_endian';
   end
end

[fid,Msg] = fopen(Path,'wt');

if fid == -1, error(Msg); end

PlyTypeNames = {'char','uchar','short','ushort','int','uint','float','double', ...
   'char8','uchar8','short16','ushort16','int32','uint32','float32','double64'};
FWriteTypeNames = {'schar','uchar','int16','uint16','int32','uint32','single','double'};
MatlabTypeNames = {'int8','uint8','int16','uint16','int32','uint32','single','double'};
PrintfTypeChar = {'%d','%u','%d','%u','%d','%u','%-.6f','%-.14e'};
IntegerDataMin = [-128,0,-2^15,-2^31,0];
IntegerDataMax = [127,255,2^16-1,2^31-1,2^32-1];

%%% write PLY header %%%
fprintf(fid,'ply\nformat %s 1.0\ncomment created by MATLAB plywrite\n',Format);
ElementNames = fieldnames(Elements);
NumElements = length(ElementNames);
Data = cell(NumElements,1);

for i = 1:NumElements
   eval(['tmp=isa(Elements.',ElementNames{i},',''struct'');']);
   
   if tmp
      eval(['PropertyNames{i}=fieldnames(Elements.',ElementNames{i},');']);
   else
      PropertyNames{i} = [];
   end
   
   if ~isempty(PropertyNames{i})
   	eval(['Data{i}{1}=Elements.',ElementNames{i},'.',PropertyNames{i}{1},';']);
      ElementCount(i) = prod(size(Data{i}{1}));
      Type{i} = zeros(length(PropertyNames{i}),1);
   else
      ElementCount(i) = 0;
   end
   
   fprintf(fid,'element %s %u\n',ElementNames{i},ElementCount(i));
   
   for j = 1:length(PropertyNames{i})
      eval(['Data{i}{j}=Elements.',ElementNames{i},'.',PropertyNames{i}{j},';']);
      
      if ElementCount(i) ~= prod(size(Data{i}{j}))
      	fclose(fid);
         error('All property data in an element must have the same length.');
      end
      
      if iscell(Data{i}{j})
         Type{i}(j) = 9;
         Data{i}{j} = Data{i}{j}{1};
      end
      
      for k = 1:length(MatlabTypeNames)
      	if isa(Data{i}{j},MatlabTypeNames{k})
         	Type{i}(j) = Type{i}(j) + k;
	         break;
         end
      end
      
      if ~rem(Type{i}(j),9)
         fclose(fid);
         error('Unsupported data structure.');
      end
      
      % try to convert float data to integer data
      if Type{i}(j) <= 8 			% array data
         if any(strcmp({'single','double'},MatlabTypeNames{Type{i}(j)}))
            if ~any(floor(Data{i}{j}) ~= Data{i}{j})		% data is integer
               MinValue = min(min(Data{i}{j}));
               MaxValue = max(max(Data{i}{j}));
               
               % choose smallest possible integer data format
               tmp = max(min(find(MinValue >= IntegerDataMin)),min(find(MaxValue <= IntegerDataMax)));
               
               if ~isempty(tmp)
                  Type{i}(j) = tmp;
               end
            end
         end
      else								% cell array data
         eval(['Data{i}{j}=Elements.',ElementNames{i},'.',PropertyNames{i}{j},';']);
         tmp = 1;
         
         for k = 1:prod(size(Data{i}{j}))
            tmp = tmp & all(floor(Data{i}{j}{k}) == Data{i}{j}{k});
	      end
         
         if tmp		% data is integer
	         MinValue = inf;
   	      MaxValue = -inf;
         
      	   for k = 1:prod(size(Data{i}{j}))
         	   MinValue = min(MinValue,min(Data{i}{j}{k}));
            	MaxValue = max(MaxValue,max(Data{i}{j}{k}));
	         end
            
            % choose smallest possible integer data format
            tmp = max(min(find(MinValue >= IntegerDataMin)),min(find(MaxValue <= IntegerDataMax)));
            
            if ~isempty(tmp)
               Type{i}(j) = tmp + 9;
            end
         end
      end
      
      % convert double to single if specified
      if rem(Type{i}(j),9) == 8 & ~strcmpi(Str,'double')
      	Type{i}(j) = Type{i}(j) - 1;
      end
      
      if Type{i}(j) <= 8
         fprintf(fid,'property %s %s\n',PlyTypeNames{Type{i}(j)},PropertyNames{i}{j});
      else
         fprintf(fid,'property list uchar %s %s\n',PlyTypeNames{Type{i}(j)-9},PropertyNames{i}{j});
      end
   end
end

fprintf(fid,'end_header\n');

switch Format
case 'ascii'
   Format = 0;
case 'binary_little_endian'
   fclose(fid);
   fid = fopen(Path,'a','ieee-le');
   Format = 1;
case 'binary_big_endian'
   fclose(fid);
   fid = fopen(Path,'a','ieee-be');
   Format = 2;
end

for i = 1:NumElements
   if ~isempty(PropertyNames{i})
   	if ~Format										% write ASCII data
      	for k = 1:ElementCount(i)
         	for j = 1:length(PropertyNames{i})
            	if Type{i}(j) <= 8
               	fprintf(fid,[PrintfTypeChar{Type{i}(j)},' '],Data{i}{j}(k));
            	else
               	fprintf(fid,'%u%s ',length(Data{i}{j}{k}),sprintf([' ',PrintfTypeChar{Type{i}(j)-9}],Data{i}{j}{k}));
					end
            end
            
         	fprintf(fid,'\n');
      	end
   	else												% write binary data
      	if all(Type{i} <= 8) & all(Type{i} == Type{i}(1))
         	% property data without list types (fast)
         	tmp = zeros(length(PropertyNames{i}),ElementCount(i));
         
         	for j = 1:length(PropertyNames{i})
            	tmp(j,:) = Data{i}{j}(:)';
         	end
         
         	fwrite(fid,tmp,FWriteTypeNames{Type{i}(j)});
     		elseif all(Type{i} > 8)
      		% only list types
         	Type{i} = Type{i} - 9;
            
         	if length(PropertyNames{i}) == 1
         		% only one list property
            	tmp = FWriteTypeNames{Type{i}(1)};
            
            	for k = 1:ElementCount(i)
               	fwrite(fid,length(Data{i}{1}{k}),'uchar');
               	fwrite(fid,Data{i}{1}{k},tmp);
            	end
         	else
         		% multiple list properties
	         	for k = 1:ElementCount(i)
   					for j = 1:length(PropertyNames{i})
      					fwrite(fid,length(Data{i}{j}{k}),'uchar');
                  	fwrite(fid,Data{i}{j}{k},FWriteTypeNames{Type{i}(j)});
						end
            	end
         	end
      	else
      		% mixed type
		      for k = 1:ElementCount(i)
   		      for j = 1:length(PropertyNames{i})
      		      if Type{i}(j) <= 8
         		      fwrite(fid,Data{i}{j}(k),FWriteTypeNames{Type{i}(j)});
	         	   else
   	         	   fwrite(fid,length(Data{i}{j}{k}),'uchar');
                     fwrite(fid,Data{i}{j}{k},FWriteTypeNames{Type{i}(j)-9});
                  end
					end
	         end
         end
      end
   end
end

fclose(fid);
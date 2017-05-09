function [V,F,UV,TF,N,NF] = readOBJ(filename,varargin)
  % READOBJ reads an OBJ file with vertex/face information
  %
  % [V,F,UV,TF,N,NF] = readOBJ(filename)
  % [V,F,UV,TF,N,NF] = readOBJ(filename,'ParameterName',ParameterValue,...)
  %
  % Input:
  %  filename  path to .obj file
  %  Optional:
  %    'Quads' whether to output face information in X by 4 matrices (faces
  %      with degree larger than 4 are still triangulated). A trailing zero
  %      will mean a triangle was read.
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #UV by 2 list of texture coordinates
  %  TF  #F by 3 list of triangle texture coordinates
  %  N  #N by 3 list of normals
  %  NF  #F by 3 list of triangle corner normal indices into N
  %
  % Examples:
  %   % read a quad/triangle mesh and display it
  %   [V,F] = readOBJ('quads.obj','Quads',true);
  %   % Turn triangles into degenerate quads 
  %   DF = (F==0).*F(:,[4 1 2 3])+(F~=0).*F;
  %   trisurf(DF,V(:,1),V(:,2),V(:,3));
  %
  %
  % See also: load_mesh, readOBJfast, readOFF

  % default values
  quads = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Quads'}, ...
    {'quads'});
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

  numv = 0;
  numf = 0;
  numuv = 0;
  numtf = 0;
  numn = 0;
  numnf = 0;

  % simplex size
  if quads
    ss = 4;
  else
    ss = 3;
  end
  % Amortized array allocation
  V = zeros(10000,3);
  F = zeros(10000,ss);
  UV = zeros(10000,3);
  TF = zeros(10000,ss);
  N = zeros(10000,3);
  NF = zeros(10000,ss);

  triangulated = false;
  all_ss = true;
  fp = fopen( filename, 'r' );
  type = fscanf( fp, '%s', 1 );
  count = 0;
  while strcmp( type, '' ) == 0
      line = fgets(fp);
      if strcmp( type, 'v' ) == 1
          v = sscanf( line, '%lf %lf %lf' );
          numv = numv+1;
          if(numv>size(V,1))
            V = cat(1,V,zeros(10000,3));
          end
          V(numv,:) = [v(1:3)'];
      elseif strcmp( type, 'vt')
          v = sscanf( line, '%f %f %f' );
          numuv = numuv+1;
          if size(UV,2)>2 && length(v) == 2
              UV = UV(:,1:2);
          end
          if(numuv>size(UV,1))
            UV = cat(1,UV,zeros(10000,length(v)));
          end
          UV(numuv,:) = [v'];
      elseif strcmp( type, 'vn')
          n = sscanf( line, '%f %f %f' );
          numn = numn+1;
          if(numn>size(N,1))
            N = cat(1,N,zeros(10000,3));
          end
          N(numn,:) = [n'];
      elseif strcmp( type, 'f' ) == 1
          [t, count] = sscanf(line,'%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');
          if (count>2)
              tf = t(2:3:end);
              nf = t(3:3:end);
              t = t(1:3:end);
          else
              [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
              if (count>2)
                  tf = t(2:2:end);
                  t = t(1:2:end);
                  nf = -ones(numel(t),1);
              else
                [t, count] = sscanf(line, '%d//%d %d//%d %d//%d %d//%d %d//%d');
                if (count>2)
                    nf = t(2:2:end);
                    t = t(1:2:end);
                    tf = -ones(numel(t),1);
                else
                    [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
                    if (count>2)
                      tf = -ones(numel(t),1);
                      nf = -ones(numel(t),1);
                    else
                      [t, count] = sscanf( line, '%d// %d// %d// %d// %d// %d// %d// %d// %d// %d// %d//\n' );
                      tf = -ones(numel(t),1);
                      nf = -ones(numel(t),1);
                    end
                end
              end
          end
          t = t + (t<0).*   (numv+1);
          tf = tf + (tf<0).*(numuv+1);
          nf = nf + (nf<0).*(numn+1);

          assert(numel(t) == numel(tf));
          assert(numel(t) == numel(nf));
          if numel(t) > ss
            if ~triangulated
              warning('Trivially triangulating high degree facets');
            end
            triangulated = true;
          end
          j = 2;
          i = 1;
          %Vt = V(t,:);
          %[~,A] = affine_fit(Vt);
          %VtA = Vt*A;
          %VtA0 = Vt*A;
          %[~,alpha] = curvature(VtA);
          %flip = -sign(sum(alpha));
          %E = [1:size(VtA,1);2:size(VtA,1) 1]';
          %[dV,dF] = triangle(VtA,E,[]);
          %if size(dF,1)>2
          %  tsurf(dF,dV);
          %  hold on;
          %  plot(VtA0([1:end 1],1),VtA0([1:end 1],2),'LineWidth',3);
          %  hold off
          %  pause
          %end
          while true
            if numel(t) > ss
              corners = [1 2 3];

              %plot(VtA0([1:end 1],1),VtA0([1:end 1],2));
              %hold on;
              %plot(VtA([1:3],1),VtA([1:3],2),'LineWidth',3);
              %hold off;
              %expand_axis(2);
              %pause;

              %[~,alpha] = curvature(VtA,[1 2;2 3]);
              %alpha = flip * alpha(2);
              %others = VtA(setdiff(1:end,corners),:);
              %these = VtA(corners,:);
              %w = inpolygon(others(:,1),others(:,2),these(:,1),these(:,2));
              %alpha
              %if alpha>=0 && ~any(w)
              %  % lazy
              %  t = t([2:end 1]);
              %  VtA = VtA([2:end 1],:);
              %  continue;
              %end
            else
              if all_ss && numel(t)<ss
                warning('Small degree facet found');
                all_ss = false;
              end
              corners = 1:numel(t);
            end
            numf = numf+1;
            if(numf>size(F,1))
              F = cat(1,F,zeros(10000,ss));
            end
            F(numf,1:numel(corners)) = [t(corners)'];
            numtf = numtf+1;
            if(numtf>size(TF,1))
              TF = cat(1,TF,zeros(10000,ss));
            end
            TF(numtf,1:numel(corners)) = [tf(corners)'];
            numnf = numnf+1;
            if(numnf>size(NF,1))
              NF = cat(1,NF,zeros(10000,ss));
            end
            NF(numnf,1:numel(corners)) = [nf(corners)'];
            if numel(t) <= ss
              break;
            end
            t = t([1 3:end]);
            tf = tf([1 3:end]);
            nf = nf([1 3:end]);
            %VtA = VtA([1 3:end],:);
            if numel(t) < 3
              break;
            end
          end
      elseif strcmp( type, '#' ) == 1
          %fscanf( fp, '%s\n', 1 );
          % ignore line
      end
      type = fscanf( fp, '%s', 1 );
  end
  fclose( fp );

  %try
  %    F = cell2mat(F);
  %end
  V = V(1:numv,:);
  F = F(1:numf,:);
  UV = UV(1:numuv,:);
  TF = TF(1:numtf,:);
  N = N(1:numn,:);
  NF = NF(1:numnf,:);

  %% transform into array if all faces have the same number of vertices

  if (size(UV,1)>0) UV = UV; end

end

classdef SVG < handle
  % members
  properties(GetAccess=public, SetAccess=public)
    P
    C
    I
    F
    S
    W
    D
    viewBox
    draw_tol = 1
  end
  properties(GetAccess=public, SetAccess=private)
    draw_data = struct('V',[],'E',[],'I',[],'K',[],'T',[],'F',[]);
  end

  methods(Access=public)
    function this = SVG(varargin)
      % if the input is a string
      if nargin == 1 && ischar(varargin{1})
        filename = varargin{1};
        [this.P,this.C,this.I,this.F,this.S,this.W,this.D,this.viewBox] = ...
          readSVG_cubics(filename);
        bbd = normrow(max(this.P)-min(this.P));
        this.draw_tol = bbd*0.001;
      else
        error('Not supported');
      end
    end

    function draw(this,draw_tol)
      if nargin < 2
        draw_tol = this.draw_tol;
      end
      if this.draw_tol ~= draw_tol || isempty(this.draw_data.V)
        % fprintf('building cached data...\n');
        this.draw_tol = draw_tol;
        S = sparse(1:numel(this.I),this.I,1,numel(this.I),max(this.I));
        this.draw_data.V = {};
        this.draw_data.F = {};
        this.draw_data.E = {};
        for i = 1:max(this.I)
          Ii = find(S(:,i));
          [Vi,Ei] = spline_to_poly(this.P,this.C(Ii,:),this.draw_tol);
          EIi = repmat(i,size(Ei,1),1);
          if all(isnan(this.F(i,:))) || size(Ei,1) < 3 
            Fi = [];
          else
            [Vi,Fi] = triangulate_interior(Vi,Ei);
            Ei = outline(Fi);
          end
          if all(isnan(this.S(i,:))) || isnan(this.W(i))
            Ei = [];
          end
          FIi = repmat(i,size(Fi,1),1);

          this.draw_data.V{i} = Vi;
          this.draw_data.E{i} = Ei;
          this.draw_data.F{i} = Fi;
        end
      else
        %fprintf('using cached data...\n');
      end

      was_hold = ishold;
      if ~was_hold
        clf;
      end
      hold on;
      for i = 1:max(this.I)
        tsurf(this.draw_data.F{i},this.draw_data.V{i},'FaceVertexCData',this.F(i,:),'FaceColor','flat','EdgeColor','none');
        if ~isempty(this.draw_data.E{i})
          %tsurf(this.draw_data.E{i},this.draw_data.V{i},'FaceColor','none','EdgeColor',this.S(i,:),'LineWidth',this.W(i));
          % This should be precomputed.
          plot_taper(this.draw_data.V{i},this.draw_data.E{i},this.W(i)/2,'FaceColor',this.S(i,:),'EdgeColor','none');
        end
      end

      hold off;
      if ~was_hold
        hold off;
      end

    end

  end

end

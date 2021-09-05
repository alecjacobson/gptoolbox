function [PR,CR] = flatten_splines(P,C,I,F,S,SW,D)
  % FLATTEN_SPLINES  Given a set of splines (P,C) defining a non-white vector
  % graphics image. "flatten" the curves into a single shape. Input white-filled
  % shapes are treated as a negative space.
  %
  %
  % Inputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %   I  #C list of indices into total number of svg objects
  %   F  #max(I) list of fill colors, NaN means none
  %   S  #max(I) list of stroke colors, NaN means none
  %   W  #max(I) list of widths colors, NaN means none
  %   D  #max(I) list of display, 0 means 'none', 1 otherwise
  % Outputs:
  %   PR  #PR by 2 list of control point locations
  %   CR  #CR by 4 list of indices into PR of cubics
  %

  
  % Example:
  %   [P,C,I,F,S,SW,D] = readSVG_cubics('~/Dropbox/models/rotor.svg');
  %   [PR,CR] = flatten_splines(P,C,I,F,S,SW,D);
  %   plot_spline(PR,CR);
  %   set(gca,'YDir','reverse')
  %   axis equal;

  is_filled = ~any(isnan(F),2);
  is_stroked = (~any(isnan(S),2)) & SW>0;
  [~,~,UJ] = unique(P,'rows');
  U = sparse(UJ,1:size(P,1),1);
  
  tol = 1e-7;
  
  PR = P;
  CR = [];
  IR = [];
  for ii = 1:max(I)
    if ~D(ii)
      continue;
    end
    if ~is_filled(ii)
      % skip non-fills
      continue;
    end
    Cii = C(I==ii,:);
    Pii = P;
    %if ~(Cii(1) == Cii(end) || all(P(Cii(1),:) == P(Cii(end),:)))
    if any(U*sparse(Cii(:,[1 4]),1,repmat([1 -1],size(Cii,1),1),size(P,1),1))
      % close it up!
      Cii = [Cii;Cii(end) size(Pii,1)+[1 2] Cii(1)];
      Pii = [Pii;interp1([0 1],Pii(Cii([1 end]),:),[1/3 2/3])];
    end
    A = spline_area(Pii,Cii);
    % if negative area, flip 
    if A<0
      Cii = fliplr(Cii);
    end
    A = spline_area(Pii,Cii);
  
    if isempty(CR) && ~all( F(ii,:) == 255)
      % first fill just keep
      CR = Cii;
      IR = repmat(ii,size(Cii,1),1);
    else
      PA = PR;
      CA = CR;
      IA = IR;
      PB = Pii;
      CB = Cii;
      IB = repmat(ii,size(Cii,1),1);
      if all( F(ii,:) == 255) 
        operation = 'minus';
      else
        operation = 'union';
      end
      [PAB,CAB,DAB,JAB] = trim_with_spline(PA,CA,PB,CB,tol);
      IAB = IA(JAB);
      [PBA,CBA,DBA,JBA] = trim_with_spline(PB,CB,PA,CA,tol);
      IBA = IB(JBA);
      switch operation
      case 'union'
        % Keep only those that are outside
        CAB = CAB(~DAB,:);
        IAB = IAB(~DAB,:);
        % Keep only those that are outside
        CBA = CBA(~DBA,:);
        IBA = IBA(~DBA,:);
      case 'minus'
        % Keep only those that are outside
        CAB = CAB(~DAB,:);
        IAB = IAB(~DAB,:);
        % Keep those that are inside yet flipped
        CBA = fliplr(CBA(DBA,:));
        IBA = IBA(DBA,:);
      end
      PBA_new = PBA(size(P,1)+1:end,:);
      CBA(CBA>size(P,1)) = (CBA(CBA>size(P,1))-size(P,1))+size(PAB,1);
      PR = [PAB;PBA_new];
      CR = [CAB;CBA];
      IR = [IA;IB];
    end
  end
  
end

function [FF] = triangulate_poly_pyramid(V,E,F)
  % TRIANGULATE_POLY_PYRAMID Triangulate a PLC described by (E,F) as output by
  % the SVR or Pyramid program
  %
  % [FF] = triangulate_poly_pyramid(E,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2list of segment indices, indexing V
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets,  each facet is a polygon
  %       **NOTE: facets index E *not* V, contrary to typical (V,F) meshes and
  %       contrary to writePOLY_tetgen prototype
  %     OR
  %   F  #F by constant-degree list of facets
  % Outputs:
  %   FF  #FF by 3 list of triangulated facets indexing V
  % 


  FF = [];
  if isstruct(F)
    for f = 1:numel(F.facets)
      %FV = outline_loop(E(F.facets{f},:),'Unoriented',true);
      %% Be sure we have a row vector
      %FV = FV(:)';

      % FACETS AREN"T IN ORDER THIS WON"T WORK
      FV = zeros(numel(F.facets{f}),1);
      FV = E(F.facets{f}(1),:);
      for e = 2:numel(F.facets{f})-1
        if FV(e) == E(F.facets{f}(e),1)
          FV(e+1) = E(F.facets{f}(e),2);
        elseif FV(e) == E(F.facets{f}(e),2)
          FV(e+1) = E(F.facets{f}(e),1);
        elseif e == 2
          % Shit, we might have gotten the order of the first edge wrong
          FV = fliplr(FV);
          if FV(e) == E(F.facets{f}(e),1)
            FV(e+1) = E(F.facets{f}(e),2);
          elseif FV(e) == E(F.facets{f}(e),2)
            FV(e+1) = E(F.facets{f}(e),1);
          else
            assert(false)
          end
        else
          assert(false);
        end
      end
      e = numel(F.facets{f});
      assert( (FV(e) == E(F.facets{f}(e),1)) || FV(e) == E(F.facets{f}(e),2));

      [I,IA,FVF] = unique(FV);
      %[I,IA,EVF] = unique(E(F.facets{f},:));
      %EVF = reshape(EVF,size(E(F.facets{f},:)));
      VF = V(I,:);
      [C,S] = princomp(VF);
      % Should be 2D
      assert(max(abs(S(:,3)))<1e-6);
      TS = delaunay(S(:,1),S(:,2));
      BC = barycenter(S(:,1:2),TS);
      %O = EVF;
      O = [FVF(1:end);FVF([2:end 1])]';
      w = abs(winding_number(S(:,1:2),O,BC)/(2*pi));

      %set(tsurf(TS,[S(:,1:2)],1),'CData',w);
      %colorbar
      %hold on;
      %plot([S(O(:,1),1) S(O(:,2),1)]',[S(O(:,1),2) S(O(:,2),2)]','-ok','LineWidth',4)
      %hold off;
      %axis equal;
      %input('');

      TS = TS(w>0.5,:);
      FF = [FF;I(TS)];
    end
  else
    error
  end

end

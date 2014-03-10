%http://ltfat.sourceforge.net/

  function [poly] = mask2poly(mask, indTime, indFreq)
    % Convert region mask to region of interest (ROI) polygon
    % !!! I'm not shure that this algorithm works in every case

    assert(~isempty(mask))
    maskRef = mask;
    
    mask = zeros(size(mask, 1)+2, size(mask, 2)+2);
    mask(2:end-1,2:end-1) = maskRef;

    [ind1, ind2] = find(mask);

    % position of segments on the edge for first dimension
    seg1 = sparse(size(mask,1)+1, size(mask,2));
    % position of segments on the edge for second dimension
    seg2 = sparse(size(mask,1), size(mask,2)+1) ;

    for n = 1:length(ind1)
      if ~mask(ind1(n)-1, ind2(n))
        seg1(ind1(n), ind2(n)) = 1;
      end
      if ~mask(ind1(n)+1, ind2(n))
        seg1(ind1(n)+1, ind2(n)) = 1;
      end
      if ~mask(ind1(n), ind2(n)-1)
        seg2(ind1(n), ind2(n)) = 1;
      end
      if ~mask(ind1(n), ind2(n)+1)
        seg2(ind1(n), ind2(n)+1) = 1;
      end
    end

    ind = 0;

    while nnz(seg1)>1
      curve = chainSeg;
      if ~isempty(curve)
        ind = ind+1;
        poly(ind).x = curve(:,2);
        poly(ind).y = curve(:,1);
        poly(ind).hole = true;
      end
    end

    poly(1).hole = false;

    function curve = chainSeg()
    % chaining of segments

      % choose one start segment
      [ind1, ind2] = find(seg1);
      ind1 = ind1(1);
      ind2 = ind2(1);
      curve = [ind1, ind2; ind1, ind2+1];

      seg1(ind1, ind2) = 0;

      % pos is a variable to remember in which direction we are
      % progressing
      % [1 1] if we're going in a growing index number
      % [-1 0] otherwise
      pos = [1, 1];

      % precise if the last added segment comes from seg1 or seg2
      last1 = true;

      while true
        if last1
          % the last segment was from dimension 1
          if seg1(ind1, ind2+pos(1))
            ind2 = ind2+pos(1);
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2+pos(2)];
          elseif seg2(ind1-1, ind2+pos(2))
            ind1 = ind1-1;
            ind2 = ind2+pos(2);
            pos = [-1, 0];
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1, ind2];
            last1 = false;
          elseif seg2(ind1, ind2+pos(2))
            ind2 = ind2+pos(2);
            pos = [1, 1];
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1+1, ind2];
            last1 = false;
          else
            break;
          end
        else
          % the last segment was from dimension 2
          if seg2(ind1+pos(1), ind2)
            ind1 = ind1+pos(1);
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1+pos(2), ind2];
          elseif seg1(ind1+pos(2), ind2-1)
            ind1 = ind1+pos(2);
            ind2 = ind2-1;
            pos = [-1, 0];
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2];
            last1 = true;
          elseif seg1(ind1+pos(2), ind2)
            ind1 = ind1+pos(2);
            pos = [1, 1];
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2+1];
            last1 = true;
          else
            break;
          end
        end
      end

      curve = curve-1.5;

      if size(curve, 1)==2
        % the first segment couldn't be linked to some other segment,
        % there is no polygon
        curve = [];
      end

    end

  end
% Copyright (c) 2011, Nikolay S.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

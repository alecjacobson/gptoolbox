function [ps,ix] = dpsimplify(p,tol)

% Recursive Douglas-Peucker Polyline Simplification, Simplify
%
% [ps,ix] = dpsimplify(p,tol)
%
% dpsimplify uses the recursive Douglas-Peucker line simplification 
% algorithm to reduce the number of vertices in a piecewise linear curve 
% according to a specified tolerance. The algorithm is also know as
% Iterative Endpoint Fit. It works also for polylines and polygons
% in higher dimensions.
%
% In case of nans (missing vertex coordinates) dpsimplify assumes that 
% nans separate polylines. As such, dpsimplify treats each line
% separately.
%
% For additional information on the algorithm follow this link
% http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
%
% Input arguments
%
%     p     polyline n*d matrix with n vertices in d 
%           dimensions.
%     tol   tolerance (maximal euclidean distance allowed 
%           between the new line and a vertex)
%
% Output arguments
%
%     ps    simplified line
%     ix    linear index of the vertices retained in p (ps = p(ix))
%
% Examples
%
% 1. Simplify line 
%
%     tol    = 1;
%     x      = 1:0.1:8*pi;
%     y      = sin(x) + randn(size(x))*0.1;
%     p      = [x' y'];
%     ps     = dpsimplify(p,tol);
%
%     plot(p(:,1),p(:,2),'k')
%     hold on
%     plot(ps(:,1),ps(:,2),'r','LineWidth',2);
%     legend('original polyline','simplified')
%
% 2. Reduce polyline so that only knickpoints remain by 
%    choosing a very low tolerance
%
%     p = [(1:10)' [1 2 3 2 4 6 7 8 5 2]'];
%     p2 = dpsimplify(p,eps);
%     plot(p(:,1),p(:,2),'k+--')
%     hold on
%     plot(p2(:,1),p2(:,2),'ro','MarkerSize',10);
%     legend('original line','knickpoints')
%
% 3. Simplify a 3d-curve
% 
%     x = sin(1:0.01:20)'; 
%     y = cos(1:0.01:20)'; 
%     z = x.*y.*(1:0.01:20)';
%     ps = dpsimplify([x y z],0.1);
%     plot3(x,y,z);
%     hold on
%     plot3(ps(:,1),ps(:,2),ps(:,3),'k*-');
%
%
%
% Author: Wolfgang Schwanghart, 13. July, 2010.
% w.schwanghart[at]unibas.ch

% Alec: Perhaps this was written when MATLAB had fewer features. There are many
% bizarre uses of cells/anonymous functions. Probably the performance could be
% increased a lot.
%


if nargin == 0
    help dpsimplify
    return
end

error(nargchk(2, 2, nargin))

% error checking
if ~isscalar(tol) || tol<0;
    error('tol must be a positive scalar')
end


% nr of dimensions
nrvertices    = size(p,1); 
dims    = size(p,2);

% Alec: this is insane. If a and b are both zero (lying on an axis) then this
% will be NaN.
%% anonymous function for starting point and end point comparision
%% using a relative tolerance test
%compare = @(a,b) abs(a-b)/max(abs(a),abs(b)) <= eps;
% Alec: Handle a=0 & b=0 case.
compare = @(a,b) abs(a-b)/max(max(abs(a),abs(b)),eps) <= eps;

% what happens, when there are NaNs?
% NaNs divide polylines.
Inan      = any(isnan(p),2);
% any NaN at all?
Inanp     = any(Inan);

% if there is only one vertex
if nrvertices == 1 || isempty(p);
    ps = p;
    ix = 1;

% if there are two 
elseif nrvertices == 2 && ~Inanp;
    % when the line has no vertices (except end and start point of the
    % line) check if the distance between both is less than the tolerance.
    % If so, return the center.
    if dims == 2;
        d    = hypot(p(1,1)-p(2,1),p(1,2)-p(2,2));
    else
        d    = sqrt(sum((p(1,:)-p(2,:)).^2));
    end
    
    if d <= tol;
        ps = sum(p,1)/2;
        ix = 1;
    else
        ps = p;
        ix = [1;2];
    end
    
elseif Inanp;
    
    % case: there are nans in the p array
    % --> find start and end indices of contiguous non-nan data
    Inan = ~Inan;
    sIX = strfind(Inan',[0 1])' + 1; 
    eIX = strfind(Inan',[1 0])'; 
 
    if Inan(end)==true;
        eIX = [eIX;nrvertices];
    end
    
    if Inan(1);
        sIX = [1;sIX];
    end
    
    % calculate length of non-nan components
    lIX = eIX-sIX+1;   
    % put each component into a single cell
    c   = mat2cell(p(Inan,:),lIX,dims);
    
    % now call dpsimplify again inside cellfun. 
    if nargout == 2;
        [ps,ix]   = cellfun(@(x) dpsimplify(x,tol),c,'uniformoutput',false);
        ix        = cellfun(@(x,six) x+six-1,ix,num2cell(sIX),'uniformoutput',false);
    else
        ps   = cellfun(@(x) dpsimplify(x,tol),c,'uniformoutput',false);
    end
    
    % write the data from a cell array back to a matrix
    ps = cellfun(@(x) [x;nan(1,dims)],ps,'uniformoutput',false);    
    ps = cell2mat(ps);
    ps(end,:) = [];
    
    % ix wanted? write ix to a matrix, too.
    if nargout == 2;
        ix = cell2mat(ix);
    end
    
       
else
    

% if there are no nans than start the recursive algorithm
ixe     = size(p,1);
ixs     = 1;

% logical vector for the vertices to be retained
% Alec: This is being updated by simplifyrec as a side-effect
I   = true(ixe,1);

% call recursive function
p   = simplifyrec(p,tol,ixs,ixe);
ps  = p(I,:);

% if desired return the index of retained vertices
if nargout == 2;
    ix  = find(I);
end

end

% _________________________________________________________
function p  = simplifyrec(p,tol,ixs,ixe)
    
    % check if startpoint and endpoint are the same 
    % better comparison needed which included a tolerance eps
    % Alec: why use cells for this? Are ixs and ixe ever more than length = 1?
    
    c1 = num2cell(p(ixs,:));
    c2 = num2cell(p(ixe,:));   
    
    % same start and endpoint with tolerance
    sameSE = all(cell2mat(cellfun(compare,c1(:),c2(:),'UniformOutput',false)));

    
    if sameSE; 
        % calculate the shortest distance of all vertices between ixs and
        % ixe to ixs only
        if dims == 2;
            d    = hypot(p(ixs,1)-p(ixs+1:ixe-1,1),p(ixs,2)-p(ixs+1:ixe-1,2));
        else
            d    = sqrt(sum(bsxfun(@minus,p(ixs,:),p(ixs+1:ixe-1,:)).^2,2));
        end
    else    
        % calculate shortest distance of all points to the line from ixs to ixe
        % subtract starting point from other locations
        pt = bsxfun(@minus,p(ixs+1:ixe,:),p(ixs,:));

        % end point
        a = pt(end,:)';

        beta = (a' * pt')./(a'*a);
        b    = pt-bsxfun(@times,beta,a)';
        if dims == 2;
            % if line in 2D use the numerical more robust hypot function
            d    = hypot(b(:,1),b(:,2));
        else
            d    = sqrt(sum(b.^2,2));
        end
    end
    
    % identify maximum distance and get the linear index of its location
    [dmax,ixc] = max(d);
    ixc  = ixs + ixc; 
    
    % if the maximum distance is smaller than the tolerance remove vertices
    % between ixs and ixe
    if dmax <= tol;
        if ixs ~= ixe-1;
            I(ixs+1:ixe-1) = false;
        end
    % if not, call simplifyrec for the segments between ixs and ixc (ixc
    % and ixe)
    else   
        p   = simplifyrec(p,tol,ixs,ixc);
        p   = simplifyrec(p,tol,ixc,ixe);

    end

end
end

% Copyright (c) 2009, Wolfgang Schwanghart
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

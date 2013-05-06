function [M,IND] = combn(V,N)
% COMBN - all combinations of elements
%   M = COMBN(V,N) returns all combinations of N elements of the elements in
%   vector V. M has the size (length(V).^N)-by-N.
%
%   [M,I] = COMBN(V,N) also returns the index matrix I so that M = V(I).
%
%   V can be an array of numbers, cells or strings.
%
%   Example:
%     M = COMBN([0 1],3) returns the 8-by-3 matrix:
%       0     0     0
%       0     0     1
%       0     1     0
%       0     1     1
%       ...
%       1     1     1
%
%   All elements in V are regarded as unique, so M = COMBN([2 2],3) returns 
%   a 8-by-3 matrix with all elements equal to 2.
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N. For larger
%   values of n and N, one could loop over the output of COMBNSUB
%   retrieving one or more rows of the output at a single time.
% 
%   See also PERMS, NCHOOSEK
%        and COMBNSUB, ALLCOMB, and PERMPOS on the File Exchange

% tested in Matlab R13, R14, 2010b
% version 4.3 (apr 2013)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
% operation on floating points, according to the IEEE 754 standard, could return 
% erroneous values. His excellent solution was to add (1/2) to the values
% of A. 
% 3.2 (may 2007) changed help and error messages slightly
% 4.0 (may 2008) again a faster implementation, based on ALLCOMB, suggested on the
%     newsgroup comp.soft-sys.matlab on May 7th 2008 by "Helper". It was
%     pointed out that COMBN(V,N) equals ALLCOMB(V,V,V...) (V repeated N
%     times), ALLCMOB being faster. Actually version 4 is an improvement
%     over version 1 ... 
% 4.1 (jan 2010) removed call to FLIPLR, using refered indexing N:-1:1
%     (is faster, suggestion of Jan Simon, jan 2010), removed REPMAT, and
%     let NDGRID handle this
% 4.2 (apr 2011) corrrectly return a column vector for N = 1 (error pointed
%      out by Wilson).
% 4.3 (apr 2013) make a reference to COMBNSUB

error(nargchk(2,2,nargin)) ;

if isempty(V) || N == 0,
    M = [] ;
    IND = [] ;
elseif fix(N) ~= N || N < 1 || numel(N) ~= 1 ;
    error('combn:negativeN','Second argument should be a positive integer') ;
elseif N==1,
    % return column vectors
    M = V(:) ; 
    IND = (1:numel(V)).' ;
else
    % speed depends on the number of output arguments
    if nargout<2,
        M = local_allcomb(V,N) ;
    else
        % indices requested
        IND = local_allcomb(1:numel(V),N) ;
        M = V(IND) ;
    end
end

% LOCAL FUNCTIONS

function Y = local_allcomb(X,N)
% See ALLCOMB, available on the File Exchange
if N>1
    % create a list of all possible combinations of N elements
    [Y{N:-1:1}] = ndgrid(X) ;
    % concatenate into one matrix, reshape into 2D and flip columns
    Y = reshape(cat(N+1,Y{:}),[],N) ;
else
    % no combinations have to be made
    Y = X(:) ;
end

% =========================================================================
% Previous algorithms


% Version 3.2
%     % COMBN is very fast using a single matrix multiplication, without any
%       explicit for-loops. 
%     nV = numel(V) ;
%     % use a math trick
%     A = [0:nV^N-1]+(1/2) ;
%     B = [nV.^(1-N:0)] ;
%     IND = rem(floor((A(:) * B(:)')),nV) + 1 ;
%     M = V(IND) ;       

% Version 2.0 
%     for i = N:-1:1
%         X = repmat(1:nV,nV^(N-i),nV^(i-1));
%         IND(:,i) = X(:);
%     end
%     M = V(IND) ;

% Version 1.0
%     nV = numel(V) ;
%     % don waste space, if only one output is requested
%     [IND{1:N}] = ndgrid(1:nV) ;
%     IND = fliplr(reshape(cat(ndims(IND{1}),IND{:}),[],N)) ;
%     M = V(IND) ;


% Combinations using for-loops
% can be implemented in C or VB
% nv = length(V) ;
% C = zeros(nv^N,N) ; % declaration
% for ii=1:N,     
%     cc = 1 ;
%     for jj=1:(nv^(ii-1)),
%         for kk=1:nv,
%             for mm=1:(nv^(N-ii)),
%                 C(cc,ii) = V(kk) ;
%                 cc = cc + 1 ;
%             end
%         end
%     end
% end  


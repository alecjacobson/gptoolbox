function [c] = slerp(a,b,t)
  % normalize a and b and keep 0's
  norma = sqrt(sum(a.^2,2));
  na = a;
  na(norma>0,:) = a(norma>0,:)./repmat(norma(norma>0),1,size(a,2));

  normb = sqrt(sum(b.^2,2));
  nb = b;
  nb(normb>0,:) = b(normb>0,:)./repmat(normb(normb>0),1,size(b,2));

  %angle = acos( dot(a,b) );
  angle = acos( sum(na.*nb,2) );
  % hard case
  %if(max(pi==angle))
  %  error('SLERP: angle between vectors cannot be exactly PI...');
  %else
    % easy degenerate case
    c = na;
    %c(angle==0) = na(angle==0);
    % errrrrr
    c(angle~=0 & angle~=pi,:) = ...
      repmat( ...
        sin((1.0-t(angle~=0 & angle~=pi)).*angle(angle~=0 & angle~=pi))./ ...
        sin(angle(angle~=0 & angle~=pi)),1,size(na,2)).* ...
      na(angle~=0 & angle~=pi,:) + ...
      repmat((sin(t(angle~=0 & angle~=pi).*angle(angle~=0 & angle~=pi))./ ...
      sin(angle(angle~=0 & angle~=pi))),1,size(nb,2)).*...
      nb(angle~=0 & angle~=pi,:); 
  %end

  % lerp magnitudes
  c = c.*repmat((1.0-t).*sqrt(sum(a.^2,2)) + ...
    t.*sqrt(sum(b.^2,2)),1,size(a,2));

end

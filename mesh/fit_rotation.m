function R = fit_rotation(S)
  % FIT_ROTATION find rotation for a given covariance matrix
  %
  % R = fit_rotation(S)
  %
  % Inputs:
  %   S  n by n covariance matrix
  % Outputs:
  %   R  n by n rotation matrix closest to S
  %

  % svd 
  [su,ss,sv]=svd(S);
  R = sv*su';
  % if reflection then flip last column
  if( det(R) < 0 )
    su(:,end) = -su(:,end);
    R = sv*su';
  end
  % should definitely be rotation now
  %assert( det(R) >= 0 );
end

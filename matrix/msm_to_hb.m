function msm_to_hb ( output_filename, A, rhs, title, key, type, ifmt, job ) 

%*****************************************************************************80
%
%% MSM_TO_HB writes a MATLAB Sparse Matrix to a Harwell Boeing Sparse Matrix file.
%
%  Usage: 
%
%    For full control:
%
%      msm_to_hb ( output_filename, A, rhs, 'title', 'key', 'type', ifmt, job ) 
%
%    or, to use defaults:
%
%      msm_to_hb ( output_filename, A )
%
%  Modified:
%
%    28 April 2004
%
%  Author:
%
%    Xiaoye Li, UC Berkeley. 
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the name of the file to which the information
%    should be written.
%
%    Input, sparse matrix A, the NROW by NCOL matrix, stored in MATLAB sparse 
%    matrix format, which is to be written to the file.
%
%    Input, real RHS(NRHS,NROW), the right-hand side array, accessed only if ( 2 < JOB).
%
%    Input, string TITLE, a title of up to 72 characters.
%    TITLE defaults to 'Title'.
%
%    Input, string KEY, a key for the matrix, of up to 8 characters.
%    KEY defaults to 'Key'.
%
%    Input, string TYPE, the HB type for the matrix, of exactly 3 characters.
%    TYPE defaults to 'RUA' (real, unsymmetric, assembled).
%
%    Input, integer IFMT, specifies the output format of the numerical values.
%    * IFMT < 100 chooses the format Dxx.yy, in which yy is precisely IFMT 
%      and xx is IFMT+7.
%    * 100 < IFMT chooses the format Fxx.yy, in which the length of the mantissa 
%      yy is the integer mod(ifmt,100) and the length of the integer part is IFMT/100.
%    For example:
%    * IFMT =   4 means  E11.4   [-]x.xxxxE+ee    
%    * IFMT = 104 means  F7.4    [-]x.xxxx
%    IFMT defaults to 8.
%
%    Input, integer JOB, indicates what is to be written out.
%    * 1, write structure only, the index arrays JA and IA.
%    * 2, write structure and matrix, A, JA, IA
%    * 3, write structure, matrix, and one right hand side: A, JA, IA, RHS.
%    * K+2, write structure, matrix and K successive right-hand sides.
%    JOB defaults to 2.
%
  [ nrow, ncol ] = size ( A );
  nnzeros = nnz ( A );
  n = ncol;
%
%  Coordinate form --> compressed column format: ja,ia,a
%
  k = 0;
  ja(1) = k + 1;
  for j = 1 : n 
    [ rows, temp, vals ] = find ( A(:,j) );
    sz = size ( rows ); 
    for k2 = 1 : sz
      k = k + 1;
      ia(k) = rows(k2);
      a(k) = vals(k2);
    end
    ja(j+1) = k + 1;
  end

  fid = fopen ( output_filename, 'wt+' );

  if ( fid < 0 ); 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MSM_TO_HB - Fatal error!\n' );
    fprintf ( 1, '  Cannot open the output file %s.\n', output_filename );
    error(['Can''t open file "' output_filename '" for writing.']); 
  end;
%
%  Default format
%
  if ( nargin < 4 )
    title = 'Title';
  end

  if ( nargin < 5 )
    key = 'Key';
  end

  if ( nargin < 6 )
    type = 'RUA';
  end

  if ( nargin < 7 )
    ifmt = 8;
  end

  if ( nargin < 8 )
    job = 2;
  end
%
%  Determine the FORTRAN format for the column pointer array.
%
  len = ceil ( log10 ( 0.1 + nnzeros + 1 ) ) + 1;
  nperline = min ( floor ( 80 / len ), ncol+1 );
  ptr_len = len;
  ptr_nperline = nperline;
  ptrcrd = floor ( ncol / nperline ) + 1;
  s1 = int2str ( len );
  s2 = int2str ( nperline );
  ptrfmt = [ '(' s1 'I' s2 ')' ];
%
%  Determine the FORTRAN format for the row index array.
%
  nperline = min ( floor ( 80 / len ), nnzeros );
  ind_len = len;
  ind_nperline = nperline;
  indcrd = floor ( ( nnzeros-1 ) / nperline ) + 1;
  s1 = int2str ( nperline );
  s2 = int2str ( len );
  indfmt = [ '(' s1 'I' s2 ')' ];
%
%  Determine the FORTRAN format for array and RHS values.
%
  valcrd = 0;
  rhscrd = 0;
  c_valfmt = [];

  if ( 1 < job )
    if ( 100 <= ifmt )
      ihead = floor ( ifmt / 100 );
      ifmt = ifmt - 100 * ihead;
      len = ihead + ifmt + 2;
      nperline = floor ( 80 / len );
      c_len = len;
      for i = 1 : nperline
        c_valfmt = [ c_valfmt '%' int2str ( c_len ) '.' int2str ( ifmt ) 'f' ];
      end
      valfmt = [ int2str ( nperline ) 'F' int2str ( len ) '.' int2str ( ifmt ) ];
    else
      len = ifmt + 7;
      nperline = floor ( 80 / len );
      c_len = len;
      s1 = int2str ( c_len );
      s2 = int2str ( ifmt );
      for i = 1 : nperline
        c_valfmt = [ c_valfmt '%' s1 '.' s2 'E' ];
      end
      s1 = int2str ( nperline );
      s2 = int2str ( len );
      s3 = int2str ( ifmt );
      valfmt = [ s1 'E' s2 '.' s3 ];
    end
    valcrd = floor ( ( nnzeros - 1 ) / nperline ) + 1;
    valfmt = [ '(' valfmt ')' ];
    c_valfmt = [ c_valfmt '\n' ];
  end

  nrhs = job - 2;
  if ( 1 <= nrhs )
    rhscrd = floor ( ( nrhs * nrow - 1 ) / nperline ) + 1;
  end

  totcrd = ptrcrd + indcrd + valcrd + rhscrd;
%
%  Write the header.
%
%  Line 1.
%
  t = title; 
  m = size ( t, 2 );

  for i = m+1 : 72
    t = [ t ' ' ];
  end
  fprintf ( fid, '%72s', t );
  t = key; 
  m = size ( t, 2 );
  for i = m+1 : 8
    t = [ t ' ' ]; 
  end
  fprintf ( fid, '%8s\n', t );
%
%  Line 2
%
  fprintf ( fid, '%14d%14d%14d%14d%14d\n', totcrd, ptrcrd, indcrd, valcrd, rhscrd );
%
%  Line 3
%
  t = type; 
  m = size ( t, 2 );
  for i = m+1 : 14 
    t = [ t ' ' ]; 
  end
  fprintf ( fid, '%14s', t );
  fprintf ( fid, '%14i%14i%14i%14i\n', nrow, ncol, nnzeros, nrhs );
%
%  Line 4
%
  t = ptrfmt; 
  m = size ( t, 2 );
  for i = m+1 : 16
    t = [ t ' ' ];
  end
  fprintf ( fid, '%16s', t );
  t = indfmt; 
  m = size ( t, 2 );
  for i = m+1 : 16
    t = [ t ' ' ];
  end
  fprintf ( fid, '%16s', t );
  t = valfmt;
  m = size ( t, 2 );
  for i = m+1 : 20
    t = [ t ' ' ];
  end
  fprintf ( fid, '%20s', t );
  fprintf ( fid, '%20s\n', t );
%
%  Column pointers.
%
  t = [];
  s1 = int2str ( ptr_len );
  for j = 1 : ptr_nperline
    t = [ t '%' s1 'd' ];
  end
  t = [ t '\n' ];
  fprintf ( fid, t, ja(1:ncol+1) ); 
  if ( ptr_nperline < ncol + 1 & ...
    ( ncol + 1 ) / ptr_nperline ~= floor ( ( ncol+1 ) / ptr_nperline ) )
    fprintf ( fid, '\n' );
  end
%
%  Row indices.
%
  t = [];
  s1 = int2str ( ind_len );
  for j = 1 : ind_nperline
    t = [ t '%' s1 'd' ]; 
  end
  t = [ t '\n' ]; 
  fprintf ( fid, t, ia(1:nnzeros) ); 
  if ( ind_nperline < nnzeros & ...
    nnzeros / ind_nperline ~= floor ( nnzeros / ind_nperline ) )
    fprintf ( fid, '\n' );
  end
%
%  Numerical values of nonzero elements of the matrix.
%
  if ( 2 <= job )
    if ( job == 2 )
    	fprintf ( fid, c_valfmt, a(1:nnzeros) );
    else
    	fprintf ( fid, c_valfmt, a(1:nnzeros) );
	if ( nperline < nnzeros & ...
        nnzeros/nperline ~= floor ( nnzeros / nperline ) )
	    fprintf ( fid, '\n' );
	end
	fprintf ( fid, c_valfmt, rhs(1:nrhs*nrow) );
    end
  end

  fprintf ( fid, '\n' );

  fclose ( fid );

  return
end
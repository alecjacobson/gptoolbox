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

% 
% 		   GNU LESSER GENERAL PUBLIC LICENSE
%                        Version 3, 29 June 2007
% 
%  Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
%  Everyone is permitted to copy and distribute verbatim copies
%  of this license document, but changing it is not allowed.
% 
% 
%   This version of the GNU Lesser General Public License incorporates
% the terms and conditions of version 3 of the GNU General Public
% License, supplemented by the additional permissions listed below.
% 
%   0. Additional Definitions.
% 
%   As used herein, "this License" refers to version 3 of the GNU Lesser
% General Public License, and the "GNU GPL" refers to version 3 of the GNU
% General Public License.
% 
%   "The Library" refers to a covered work governed by this License,
% other than an Application or a Combined Work as defined below.
% 
%   An "Application" is any work that makes use of an interface provided
% by the Library, but which is not otherwise based on the Library.
% Defining a subclass of a class defined by the Library is deemed a mode
% of using an interface provided by the Library.
% 
%   A "Combined Work" is a work produced by combining or linking an
% Application with the Library.  The particular version of the Library
% with which the Combined Work was made is also called the "Linked
% Version".
% 
%   The "Minimal Corresponding Source" for a Combined Work means the
% Corresponding Source for the Combined Work, excluding any source code
% for portions of the Combined Work that, considered in isolation, are
% based on the Application, and not on the Linked Version.
% 
%   The "Corresponding Application Code" for a Combined Work means the
% object code and/or source code for the Application, including any data
% and utility programs needed for reproducing the Combined Work from the
% Application, but excluding the System Libraries of the Combined Work.
% 
%   1. Exception to Section 3 of the GNU GPL.
% 
%   You may convey a covered work under sections 3 and 4 of this License
% without being bound by section 3 of the GNU GPL.
% 
%   2. Conveying Modified Versions.
% 
%   If you modify a copy of the Library, and, in your modifications, a
% facility refers to a function or data to be supplied by an Application
% that uses the facility (other than as an argument passed when the
% facility is invoked), then you may convey a copy of the modified
% version:
% 
%    a) under this License, provided that you make a good faith effort to
%    ensure that, in the event an Application does not supply the
%    function or data, the facility still operates, and performs
%    whatever part of its purpose remains meaningful, or
% 
%    b) under the GNU GPL, with none of the additional permissions of
%    this License applicable to that copy.
% 
%   3. Object Code Incorporating Material from Library Header Files.
% 
%   The object code form of an Application may incorporate material from
% a header file that is part of the Library.  You may convey such object
% code under terms of your choice, provided that, if the incorporated
% material is not limited to numerical parameters, data structure
% layouts and accessors, or small macros, inline functions and templates
% (ten or fewer lines in length), you do both of the following:
% 
%    a) Give prominent notice with each copy of the object code that the
%    Library is used in it and that the Library and its use are
%    covered by this License.
% 
%    b) Accompany the object code with a copy of the GNU GPL and this license
%    document.
% 
%   4. Combined Works.
% 
%   You may convey a Combined Work under terms of your choice that,
% taken together, effectively do not restrict modification of the
% portions of the Library contained in the Combined Work and reverse
% engineering for debugging such modifications, if you also do each of
% the following:
% 
%    a) Give prominent notice with each copy of the Combined Work that
%    the Library is used in it and that the Library and its use are
%    covered by this License.
% 
%    b) Accompany the Combined Work with a copy of the GNU GPL and this license
%    document.
% 
%    c) For a Combined Work that displays copyright notices during
%    execution, include the copyright notice for the Library among
%    these notices, as well as a reference directing the user to the
%    copies of the GNU GPL and this license document.
% 
%    d) Do one of the following:
% 
%        0) Convey the Minimal Corresponding Source under the terms of this
%        License, and the Corresponding Application Code in a form
%        suitable for, and under terms that permit, the user to
%        recombine or relink the Application with a modified version of
%        the Linked Version to produce a modified Combined Work, in the
%        manner specified by section 6 of the GNU GPL for conveying
%        Corresponding Source.
% 
%        1) Use a suitable shared library mechanism for linking with the
%        Library.  A suitable mechanism is one that (a) uses at run time
%        a copy of the Library already present on the user's computer
%        system, and (b) will operate properly with a modified version
%        of the Library that is interface-compatible with the Linked
%        Version.
% 
%    e) Provide Installation Information, but only if you would otherwise
%    be required to provide such information under section 6 of the
%    GNU GPL, and only to the extent that such information is
%    necessary to install and execute a modified version of the
%    Combined Work produced by recombining or relinking the
%    Application with a modified version of the Linked Version. (If
%    you use option 4d0, the Installation Information must accompany
%    the Minimal Corresponding Source and Corresponding Application
%    Code. If you use option 4d1, you must provide the Installation
%    Information in the manner specified by section 6 of the GNU GPL
%    for conveying Corresponding Source.)
% 
%   5. Combined Libraries.
% 
%   You may place library facilities that are a work based on the
% Library side by side in a single library together with other library
% facilities that are not Applications and are not covered by this
% License, and convey such a combined library under terms of your
% choice, if you do both of the following:
% 
%    a) Accompany the combined library with a copy of the same work based
%    on the Library, uncombined with any other library facilities,
%    conveyed under the terms of this License.
% 
%    b) Give prominent notice with the combined library that part of it
%    is a work based on the Library, and explaining where to find the
%    accompanying uncombined form of the same work.
% 
%   6. Revised Versions of the GNU Lesser General Public License.
% 
%   The Free Software Foundation may publish revised and/or new versions
% of the GNU Lesser General Public License from time to time. Such new
% versions will be similar in spirit to the present version, but may
% differ in detail to address new problems or concerns.
% 
%   Each version is given a distinguishing version number. If the
% Library as you received it specifies that a certain numbered version
% of the GNU Lesser General Public License "or any later version"
% applies to it, you have the option of following the terms and
% conditions either of that published version or of any later version
% published by the Free Software Foundation. If the Library as you
% received it does not specify a version number of the GNU Lesser
% General Public License, you may choose any version of the GNU Lesser
% General Public License ever published by the Free Software Foundation.
% 
%   If the Library as you received it specifies that a proxy can decide
% whether future versions of the GNU Lesser General Public License shall
% apply, that proxy's public statement of acceptance of any version is
% permanent authorization for you to choose that version for the
% Library.

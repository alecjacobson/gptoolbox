/* 
RES = upConv(IM, FILT, EDGES, STEP, START, STOP, RES);
  >>> See upConv.m for documentation <<<
  This is a MATLAB interface to the internal_expand function. 
  EPS, 7/96.
*/

#define V4_COMPAT
#include <matrix.h>  /* MATLAB matrices */
#include <mex.h>

#include "convolve.h"

#define notDblMtx(it) (!mxIsNumeric(it) || !mxIsDouble(it) || mxIsSparse(it) || mxIsComplex(it))

void mexFunction(int nlhs,	     /* Num return vals on lhs */
		 mxArray *plhs[],    /* Matrices on lhs      */
		 int nrhs,	     /* Num args on rhs    */
		 const mxArray *prhs[]     /* Matrices on rhs */
		 )
  {
  double *image,*filt, *temp, *result, *orig_filt;
  int x_fdim, y_fdim, x_idim, y_idim;
  int orig_x = 0, orig_y, x, y;
  int x_rdim, y_rdim;
  int x_start = 1;
  int x_step = 1;
  int y_start = 1;
  int y_step = 1;
  int x_stop, y_stop;
  const mxArray *arg0,*arg1,*arg3,*arg4,*arg6;
  double *mxMat;
  char edges[15] = "reflect1";

  if (nrhs<2) mexErrMsgTxt("requres at least 2 args.");

  /* ARG 1: IMAGE  */
  arg0 = prhs[0];
  if notDblMtx(arg0) mexErrMsgTxt("IMAGE arg must be a non-sparse double float matrix.");
  image = mxGetPr(arg0);
  x_idim = (int) mxGetM(arg0); /* X is inner index! */
  y_idim = (int) mxGetN(arg0);

  /* ARG 2: FILTER */
  arg1 = prhs[1];
  if notDblMtx(arg1) mexErrMsgTxt("FILTER arg must be non-sparse double float matrix.");  filt = mxGetPr(arg1);
  x_fdim = (int) mxGetM(arg1); 
  y_fdim = (int) mxGetN(arg1);

  /* ARG 3 (optional): EDGES */
  if (nrhs>2) 
	  {
  	  if (!mxIsChar(prhs[2]))
  	  	mexErrMsgTxt("EDGES arg must be a string.");
  	  mxGetString(prhs[2],edges,15);
	  }

  /* ARG 4 (optional): STEP */
  if (nrhs>3)
      {
      arg3 = prhs[3];
      if notDblMtx(arg3) mexErrMsgTxt("STEP arg must be double float matrix.");
      if (mxGetM(arg3) * mxGetN(arg3) != 2)
    	 mexErrMsgTxt("STEP arg must contain two elements.");
      mxMat = mxGetPr(arg3);
      x_step = (int) mxMat[0];
      y_step = (int) mxMat[1];
      if ((x_step<1) || (y_step<1))
         mexErrMsgTxt("STEP values must be greater than zero.");
      }

  /* ARG 5 (optional): START */
  if (nrhs>4)
      {
      arg4 = prhs[4];
      if notDblMtx(arg4) mexErrMsgTxt("START arg must be double float matrix.");
      if (mxGetM(arg4) * mxGetN(arg4) != 2)
	mexErrMsgTxt("START arg must contain two elements.");
      mxMat = mxGetPr(arg4);
      x_start = (int) mxMat[0];
      y_start = (int) mxMat[1];
      if ((x_start<1) || (y_start<1))
         mexErrMsgTxt("START values must be greater than zero.");
      }
  x_start--;  /* convert to standard C indexes */
  y_start--;

  /* ARG 6 (optional): STOP */
  if (nrhs>5)
      {
      if notDblMtx(prhs[5]) mexErrMsgTxt("STOP arg must be double float matrix.");
      if (mxGetM(prhs[5]) * mxGetN(prhs[5]) != 2)
    	 mexErrMsgTxt("STOP arg must contain two elements.");
      mxMat = mxGetPr(prhs[5]);
      x_stop = (int) mxMat[0];
      y_stop = (int) mxMat[1];
      if ((x_stop<x_start) || (y_stop<y_start))
         mexErrMsgTxt("STOP values must be greater than START values.");
      }
  else
      { /* default: make res dims a multiple of STEP size */
      x_stop = x_step * ((x_start/x_step) + x_idim);
      y_stop = y_step * ((y_start/y_step) + y_idim);
      }

  /* ARG 6 (optional): RESULT image */
  if (nrhs>6)
      {
      arg6 = prhs[6];
      if notDblMtx(arg6) mexErrMsgTxt("RES arg must be double float matrix.");

      /* 7/10/97: Returning one of the args causes problems with MATLAB's memory 
	 manager, so we don't return anything if the result image is passed */
      /*  plhs[0] = arg;  */
      result = mxGetPr(arg6);
      x_rdim =  (int) mxGetM(arg6); /* X is inner index! */
      y_rdim = (int) mxGetN(arg6);
      if  ((x_stop>x_rdim) || (y_stop>y_rdim))
	mexErrMsgTxt("STOP values must within image dimensions.");
      }
  else
      {
      x_rdim = x_stop; 
      y_rdim = y_stop;
      /*  x_rdim = x_step * ((x_stop+x_step-1)/x_step);
          y_rdim = y_step * ((y_stop+y_step-1)/y_step);  */

      plhs[0] = (mxArray *) mxCreateDoubleMatrix(x_rdim,y_rdim,mxREAL);
      if (plhs[0] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
      result = mxGetPr(plhs[0]);
      }
	  
  if ( (((x_stop-x_start+x_step-1) / x_step) != x_idim) ||
       (((y_stop-y_start+y_step-1) / y_step) != y_idim) )
    {
      mexPrintf("Im dims: [%d %d]\n",x_idim,y_idim);
      mexPrintf("Start:   [%d %d]\n",x_start,y_start);
      mexPrintf("Step:    [%d %d]\n",x_step,y_step);
      mexPrintf("Stop:    [%d %d]\n",x_stop,y_stop);
      mexPrintf("Res dims: [%d %d]\n",x_rdim,y_rdim);
      mexErrMsgTxt("Image sizes and upsampling args are incompatible!");
    }

  /* upConv has a bug for even-length kernels when using the 
     reflect1, extend, or repeat edge-handlers */
  if ((!strcmp(edges,"reflect1") || !strcmp(edges,"extend") || !strcmp(edges,"repeat"))
      &&
      ((x_fdim%2 == 0) || (y_fdim%2 == 0)))
      {
      orig_filt = filt;
      orig_x = x_fdim; 
      orig_y = y_fdim;
      x_fdim = 2*(orig_x/2)+1;
      y_fdim = 2*(orig_y/2)+1;
      filt = mxCalloc(x_fdim*y_fdim, sizeof(double));
      if (filt == NULL)
	mexErrMsgTxt("Cannot allocate necessary temporary space");
      for (y=0; y<orig_y; y++)
	for (x=0; x<orig_x; x++)
	    filt[y*x_fdim + x] = orig_filt[y*orig_x + x];
      }

  if ((x_fdim > x_rdim) || (y_fdim > y_rdim))
    {
    mexPrintf("Filter: [%d %d], ",x_fdim,y_fdim);
    mexPrintf("Result: [%d %d]\n",x_rdim,y_rdim);
    mexErrMsgTxt("FILTER dimensions larger than RESULT dimensions.");
    }
 
  temp = mxCalloc(x_fdim*y_fdim, sizeof(double));
  if (temp == NULL)
    mexErrMsgTxt("Cannot allocate necessary temporary space");

  /*
  printf("(%d, %d), (%d, %d), (%d, %d), (%d, %d), (%d, %d), %s\n",
	 x_idim,y_idim,x_fdim,y_fdim,x_rdim,y_rdim,
	 x_start,x_step,y_start,y_step,edges);
	 */

  if (strcmp(edges,"circular") == 0)
	internal_wrap_expand(image, filt, x_fdim, y_fdim,
			     x_start, x_step, x_stop, y_start, y_step, y_stop,
			     result, x_rdim, y_rdim);
  else internal_expand(image, filt, temp, x_fdim, y_fdim,
		       x_start, x_step, x_stop, y_start, y_step, y_stop,
		       result, x_rdim, y_rdim, edges);

  if (orig_x) mxFree((char *) filt);
  mxFree((char *) temp);

  return;
  }      




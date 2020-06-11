#include <mex.h>
#include "paste.h"
#include <iostream>

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  unsigned char * IM;
  size_t h,w,c;
  if(!paste(IM,h,w,c))
  {
    mexErrMsgTxt("Clipboard doesn't contain image.");
  }
  unsigned char * I = IM;
  size_t ceff = c;
  unsigned char * A = NULL;
  bool has_alpha = false;
  if(c==4)
  {
    I = new unsigned char[w*h*ceff];
    A = new unsigned char[w*h*1];
    ceff = 3;
    has_alpha = true;
    for(size_t y = 0; y < h ; y++)
    {
      for(size_t x = 0; x < w; x++)
      {
        I[y+h*(x+0*w)] = IM[y+h*(x+0*w)];
        I[y+h*(x+1*w)] = IM[y+h*(x+1*w)];
        I[y+h*(x+2*w)] = IM[y+h*(x+2*w)];
        A[y+h*(x+0*w)] = IM[y+h*(x+3*w)];
      }
    }
  }

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt("Too many output parameters.");
    }
    case 2:
    {
      if(has_alpha)
      {
        size_t dims[] = {h,w};
        plhs[1] = mxCreateNumericArray(2,dims,mxUINT8_CLASS,mxREAL);
        unsigned char * Ap = (unsigned char *)mxGetData(plhs[1]);
        std::copy(A,A+w*h,Ap);
      }else
      {
        plhs[1] = mxCreateDoubleMatrix( 0, 0, mxREAL );
      }
      // Fall through
    }
    case 1:
    {
      size_t dims[] = {h,w,ceff};
      plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
      unsigned char * Ip = (unsigned char *)mxGetData(plhs[0]);
      std::copy(I,I+w*h*ceff,Ip);
      // Fall through
    }
    case 0: break;
  }
  delete[] I;
  delete[] IM;
  delete[] A;
}

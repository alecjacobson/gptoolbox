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
  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt("Too many output parameters.");
    }
    case 1:
    {
      size_t dims[] = {h,w,c};
      plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
      unsigned char * IMp = (unsigned char *)mxGetData(plhs[0]);
      std::copy(IM,IM+w*h*c,IMp);
      // Fall through
    }
    case 0: break;
  }
}

//***************************************************************************
//
// Matlab C routine file:  skelgrad.cpp
//
// Written 8/04 by N. Howe
//
// Input:
//   img:  binary silhouette image
//
// Output:
//   skg:  skeleton gradient transform
//   skr:  skeleton radius
//
//***************************************************************************

#include "mex.h"

//***************************************************************************

#define SQR(x) (x)*(x)
#define MIN(x,y) (((x) < (y)) ? (x):(y))
#define MAX(x,y) (((x) > (y)) ? (x):(y))
#define ABS(x) (((x) < 0) ? (-(x)):(x))
#define errCheck(a,b) if (!(a)) mexErrMsgTxt((b));
#define cutBounds(i) (((i)>0) ? (((i)<ncut) ? cut[(i)]:mxGetInf()):-mxGetInf())
#define MOD(x,n) (((x)%(n)<0) ? ((x)%(n)+(n)):((x)%(n)))

//***************************************************************************

enum direction {North, South, East, West, None};

const direction dircode[16] = {
  None, West, North, West, East, None, North, West,
  South, South, None, South, East, East, North, None
};

//************************************************************
//
//  quicksort() recursively calls itself to sort the array.
//  It splits the array around a randomly chosen pivot,
//  and calls itself on each piece.
//
//  R/W arr:  The array to sort.
//  R/O n:    The length of the array to sort.
//

void quicksort(int *arr, int n) {

  if (n > 8) {
    int pivot;
    int pivid;
    int lo = 1;    // everything to the left of lo is smaller than pivot
    int hi = n-1;  // everything to the right of hi is larger than the pivot
    int tmp;

    pivid = rand()%n;
    mxAssert((pivid >= 0),"Randomness too small");
    mxAssert((pivid < n),"Randomness too big");
    tmp = arr[0];
    arr[0] = arr[pivid];
    arr[pivid] = tmp;
    pivot = arr[0];  // store value for convenience
    while (hi != lo-1) {  // keep going until everything has been categorized.
      if (arr[lo] < pivot) {
        // smaller than pivot, so stays here
        lo++;
      } else {
        // larger than pivot, so moves to the other end
        tmp = arr[hi];
        arr[hi] = arr[lo];
        arr[lo] = tmp;
        hi--;
      }
    }
    tmp = arr[hi];
    arr[hi] = arr[0];
    arr[0] = tmp;
    quicksort(arr,hi);          // sort smaller side of array
    quicksort(arr+lo,n-lo);  // sort bigger side of array
  } else {
    // use insertion sort
    int i, j;
    int tmp;

    for (i = 0; i < n; i++) {
      for (j = i; (j > 0)&&(arr[j] < arr[j-1]); j--) {
        tmp = arr[j];
        arr[j] = arr[j-1];
        arr[j-1] = tmp;
      }
    }
  }
}
// end of quicksort()

//****************************************************************************
//
// joint_neighborhood() returns a byte value representing the immediate 
// neighborhood of the specified joint in bit form, clockwise from NW.
//
// R/O arr:  binary pixel array
// R/O i,j:  coordinates of point
// R/O nrow, ncol:  dimensions of binary image
// Returns:  bit representation of 8-neighbors
//

template <class T>
inline int joint_neighborhood(T *arr, int i, int j, int nrow, int ncol) {
  int p = i+j*nrow;
  int condition = 8*(i <= 0)+4*(j <= 0)+2*(i >= nrow)+(j >= ncol);

  switch (condition) {
  case 0:  // all points valid
    return (arr[p-nrow-1]?1:0)+(arr[p-1]?2:0)+(arr[p]?4:0)+(arr[p-nrow]?8:0);
  case 1:  // right side not valid
    return (arr[p-nrow-1]?1:0)+(arr[p-nrow]?8:0);
  case 2:  // bottom not valid
    return (arr[p-nrow-1]?1:0)+(arr[p-1]?2:0);
  case 3:  // bottom and right not valid
    return (arr[p-nrow-1]?1:0);
  case 4:  // left side not valid
    return (arr[p-1]?2:0)+(arr[p]?4:0);
  case 5:  // left and right sides not valid
    return 0;
  case 6:  // left and bottom sides not valid
    return (arr[p-1]?2:0);
  case 7:  // left, bottom, and right sides not valid
    return 0;
  case 8:  // top side not valid
    return (arr[p]?4:0)+(arr[p-nrow]?8:0);
  case 9:  // top and right not valid
    return (arr[p-nrow]?8:0);
  case 10:  // top and bottom not valid
  case 11:  // top, bottom and right not valid
    return 0;
  case 12:  // top and left not valid
    return (arr[p]?4:0);
  case 13:  // top, left and right sides not valid
  case 14:  // top, left and bottom sides not valid
  case 15:  // no sides valid
    return 0;
  default:
    mexErrMsgTxt("Impossible condition.\n");
    return -1;
  }
}

//***************************************************************************
//
// compute_skeleton_gradient() does the main computation.
// Written as a template to support multiple image types.
//
// R/O img:  the image
// R/O nrow, ncol:  image dimensions
// R/O nlhs:  number of arguments to return
// W/O plhs:  array of return arguments (0 = gradient, 1 = radius)
// 

template <class T> void
compute_skeleton_gradient(T *img, int nrow, int ncol, int nlhs, 
			  mxArray **plhs) {
  int i, j, ei, ej, inear;
  int ijunc, iedge, iseq, lastdir, mind, minjunc, pspan;
  int jnrow = nrow+1, jncol = ncol+1;
  int njunc, jhood, nedge, nnear;
  int *jx, *jy, *edgej, *seqj, *edgelen;
  int *dNE, *dNW, *dSE, *dSW, *nearj;
  bool *seenj;
  double *skg, *rad;

  // count junctions
  njunc = 0;
  for (j = 0; j < jncol; j++) {
    for (i = 0; i < jnrow; i++) {
      jhood = joint_neighborhood(img,i,j,nrow,ncol);
      if ((jhood != 0)&&(jhood != 15)) {
	njunc++;
      }
    }
  }
  //mexPrintf("Counted %d junctions.\n",njunc);
 
  // allocate scratch space
  jx = (int*)mxMalloc(njunc*sizeof(int));
  jy = (int*)mxMalloc(njunc*sizeof(int));
  seqj = (int*)mxMalloc(njunc*sizeof(int));
  edgej = (int*)mxMalloc(njunc*sizeof(int));
  seenj = (bool*)mxMalloc(jnrow*jncol*sizeof(bool));
  dNE = (int*)mxMalloc(njunc*sizeof(int));
  dNW = (int*)mxMalloc(njunc*sizeof(int));
  dSE = (int*)mxMalloc(njunc*sizeof(int));
  dSW = (int*)mxMalloc(njunc*sizeof(int));
  nearj = (int*)mxMalloc(njunc*sizeof(int));
  for (i = 0; i < jnrow*jncol; i++) {
    seenj[i] = false;
  }
  //mexPrintf("Space allocated\n");

  // register junctions
  ijunc = 0;
  nedge = 0;
  for (j = 0; j < jncol; j++) {
    for (i = 0; i < jnrow; i++) {
      jhood = joint_neighborhood(img,i,j,nrow,ncol);
      //mexPrintf("(%d,%d) neighborhood:  %d.\n",i,j,jhood);
      mxAssert(i+j*jnrow>=0&&i+j*jnrow<jnrow*jncol,"Out of bounds.");
      if ((jhood != 0)&&(jhood != 15)&&(jhood != 5)&&(jhood != 10)
	  &&!seenj[i+j*jnrow]) {
	// found new edge; traverse it
	iseq = 0;
	ei = i;
	ej = j;
	lastdir = North;
	//mexPrintf("Beginning traverse.\n");
	while (!seenj[ei+ej*jnrow]||(jhood==5)||(jhood==10)) {
	  //mexPrintf("Traversing at (%d,%d).\n",ei,ej);

	  if (!seenj[ei+ej*jnrow]) {
	    // register this junction
	    mxAssert((ijunc>=0)&&(ijunc<njunc),"Junction index error.");
	    jx[ijunc] = ej;
	    jy[ijunc] = ei;
	    edgej[ijunc] = nedge;
	    seqj[ijunc] = iseq;
	    iseq++;
	    ijunc++;
	    seenj[ei+ej*jnrow] = true;
	  }

	  // traverse clockwise
	  switch (dircode[jhood]) {
	  case North:
	    ei--;
	    lastdir = North;
	    break;
	  case South:
	    ei++;
	    lastdir = South;
	    break;
	  case East:
	    ej++;
	    lastdir = East;
	    break;
	  case West:
	    ej--;
	    lastdir = West;
	    break;
	  case None:
	    switch (lastdir) {
	    case East:  // go North
	      ei--;
	      lastdir = North;
	      break;
	    case West:  // go South
	      ei++;
	      lastdir = South;
	      break;
	    case South:  // go East
	      ej++;
	      lastdir = East;
	      break;
	    case North:  // go West
	      ej--;
	      lastdir = West;
	      break;
	    }
	    break;
	  }
	  mxAssert((ei>=0)&&(ej>=0)&&(ei<jnrow)&&(ej<jncol),"Traversed out.");
	  jhood = joint_neighborhood(img,ei,ej,nrow,ncol);
	  //mexPrintf("Traversal direction:  %d; new neighborhood:  %d.\n",
	  //    lastdir,jhood);
	}
	nedge++;
	}
    }
  }
  mxAssert(njunc == ijunc,"Junction miscount.");
  //mexPrintf("Junctions counted.\n");

  // count perimeter along each edge
  edgelen = (int*)mxMalloc(nedge*sizeof(int));
  for (iedge = 0; iedge < nedge; iedge++) {
    edgelen[iedge] = 0;
  }
  for (ijunc = 0; ijunc < njunc; ijunc++) {
    edgelen[edgej[ijunc]]++;
  }

  // create output
  plhs[0] = mxCreateNumericMatrix(nrow,ncol,mxDOUBLE_CLASS,mxREAL);
  skg = mxGetPr(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateNumericMatrix(nrow,ncol,mxDOUBLE_CLASS,mxREAL);
    rad = mxGetPr(plhs[1]);
  }
  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
      if (img[i+j*nrow]) {
	// compute distance to all junction points
	// keeping track of minimum
	mind = SQR(jnrow+jncol);
	minjunc = -1;
	for (ijunc = 0; ijunc < njunc; ijunc++) {
	  dNE[ijunc] = SQR(i-jy[ijunc])+SQR(j-jx[ijunc]);
	  dNW[ijunc] = SQR(i-jy[ijunc])+SQR(j+1-jx[ijunc]);
	  dSE[ijunc] = SQR(i+1-jy[ijunc])+SQR(j-jx[ijunc]);
	  dSW[ijunc] = SQR(i+1-jy[ijunc])+SQR(j+1-jx[ijunc]);
	  if (dNE[ijunc] < mind) {
	    mind = dNE[ijunc];
	    minjunc = ijunc;
	  }
	  if (dNW[ijunc] < mind) {
	    mind = dNW[ijunc];
	    minjunc = ijunc;
	  }
	  if (dSE[ijunc] < mind) {
	    mind = dSE[ijunc];
	    minjunc = ijunc;
	  }
	  if (dSW[ijunc] < mind) {
	    mind = dSW[ijunc];
	    minjunc = ijunc;
	  }
	}
	//mexPrintf("Point (%d,%d):  junction %d.\n",i,j,minjunc);
	mxAssert((minjunc >=0)&&(minjunc<njunc),"Bad minimum junction.");
	mxAssert((edgej[minjunc] >=0)&&(edgej[minjunc]<nedge),
		 "Bad minimum junction edge.");

	// store radius if desired
	if (nlhs > 1) {
	  rad[i+j*nrow] = mind;
	}

	// find all other junction points at minimal distance
	nnear = pspan = 0;
	for (ijunc = 0; ijunc < njunc; ijunc++) {
	  if ((dNE[ijunc] <= dNE[minjunc])||(dNW[ijunc] <= dNW[minjunc])||
	      (dSE[ijunc] <= dSE[minjunc])||(dSW[ijunc] <= dSW[minjunc])) {
	    // we have a candidate junction
	    if (edgej[ijunc] != edgej[minjunc]) {
	      pspan = -1;
	      break;
	    } else {
	      nearj[nnear] = seqj[ijunc];
	      nnear++;
	    }
	  }
	}

	if (pspan >= 0){
	  // compute perimeter span -- find largest gap and take remainder
	  quicksort(nearj,nnear);
	  //mexPrintf("Positions:  ");
	  //for (inear = 0; inear < nnear; inear++) {
	  //  mexPrintf("%d ",nearj[inear]);
	  //}
	  pspan = nearj[0]-nearj[nnear-1]+edgelen[edgej[minjunc]];
	  //mexPrintf("\nDifferences:  %d ",pspan);
	  for (inear = 1; inear < nnear; inear++) {
	    if (pspan < nearj[inear]-nearj[inear-1]) {
	      pspan = nearj[inear]-nearj[inear-1];
	    }
	    //mexPrintf("%d ",nearj[inear]-nearj[inear-1]);
	  }
	  //mexPrintf(" => Result:  %d (from %d; ep = %d).\n",
	  //    edgelen[edgej[minjunc]]-pspan,pspan,
	  //    edgelen[edgej[minjunc]]);
	  pspan = edgelen[edgej[minjunc]]-pspan;
	  skg[i+j*nrow] = pspan;
	} else {
	  skg[i+j*nrow] = mxGetInf();
	}

	//mexPrintf("Final span:  %g.\n",pspan);
      } else {
	skg[i+j*nrow] = 0;
	if (nlhs > 1) {
	  rad[i+j*nrow] = 0;
	}
      }
    }
  }

  // free space
  mxFree(jx);
  mxFree(jy);
  mxFree(seqj);
  mxFree(edgej);
  mxFree(seenj);
  mxFree(dNE);
  mxFree(dNW);
  mxFree(dSE);
  mxFree(dSW);
  mxFree(nearj);
}
// end of compute_skeleton_gradient()

//***************************************************************************
//
// Gateway driver to call the calculation from Matlab.
//
// This is the Matlab entry point.
// 

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int nrow, ncol;
  double *img;

  // check for proper number of arguments
  errCheck(nrhs == 1,"Exactly one input argument required.");
  errCheck(nlhs <= 2,"Too many output arguments.");

  // check format of arguments
  errCheck(mxIsUint8(prhs[0])||mxIsLogical(prhs[0])||mxIsDouble(prhs[0]),
	   "Input must be double binary image.");
  nrow = mxGetM(prhs[0]);
  ncol = mxGetN(prhs[0]);
  img = mxGetPr(prhs[0]);

  // process image
  if (mxIsDouble(prhs[0])) {
    compute_skeleton_gradient(img,nrow,ncol,nlhs,plhs);
  } else {
    compute_skeleton_gradient((unsigned char *)img,nrow,ncol,nlhs,plhs);
  }
}
// end of mexFunction()

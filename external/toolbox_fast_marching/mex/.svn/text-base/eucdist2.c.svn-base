/*
 * Copyright 1993-2003 The MathWorks, Inc.
 * $Revision: 1.5.4.3 $  $Date: 2004/02/01 21:44:32 $
 */

/*
 * EUCDIST2 MEX-file
 *
 * Distance transform, Euclidean version
 *
 * D = EUCDIST2(BW) computes the Euclidean distance transform on the input 
 * binary image BW.  Specifically, it computes the distance to the nearest
 * nonzero-valued pixel.  BW must be 2-D.
 *
 * [D,L] = EUCDIST2(BW) returns a linear index array L representing a
 * nearest-neighbor map.  L(r,c) is the linear index of the nonzero-valued
 * element of BW closest to (r,c).
 *
 * Input-output specs
 * ------------------
 * BW:    2-D logical matrix
 *        empty allowed
 *
 * D:     2-D real double matrix, same size as BW
 *        contains nonnegative values
 *
 * Algorithm notes
 * ---------------
 * Heinz Breu, Joseph Gil, David Kirkpatrick, and Michael Werman,
 * "Linear Time Euclidean Distance Transform Algorithms," IEEE
 * Transactions on Pattern Analysis and Machine Intelligence,
 * vol 17 no 5 May 1995, pp. 529-533.  Here we use the second
 * algorithm, "Transform Based on Partial Voronoi Diagram Construction,"
 * pp. 531-533.
 */

#include <math.h>
#include "mex.h"

typedef struct pixel_T 
{
    int x;
    int y;
}
pixel_T;

#define FAKE_COL_VAL -1

/*
 * From "Function Remove(u,v,w,r)", middle of second column, page 532.
 * Reversed the terms x and y compared to the equations in the paper.
 *
 * This function returns true if the Voronoi cell centered on
 * v does not intersect the line x=c.
 */
bool remove_candidate(pixel_T u, pixel_T v, pixel_T w, int c)
{
    /*
     * Perform calculations in double, since this operation is
     * sensitive to overflow when done in int32.  See g198383.
     * -sle, 2004/01/13
     */
    double term1;
    double term2;
    double ux = (double) u.x;
    double uy = (double) u.y;
    double vx = (double) v.x;
    double vy = (double) v.y;
    double wx = (double) w.x;
    double wy = (double) w.y;
    
    term1 = (wy - vy) * (vy*vy - uy*uy - 2*c*(vx - ux) +
                           vx*vx - ux*ux);

    term2 = (vy - uy) * (wy*wy - vy*vy - 2*c*(wx - vx) +
                           wx*wx - vx*vx);

    return term1 >= term2;
}

/*
 * Update the candidates array for column c.  Return the number
 * of candidates.  Also update the col_values array, which is
 * an array of maximal (or minimal) column values for pixels in
 * each row, up to and including the current column c.  bw is a pointer
 * to the binary image array; M is the number of rows of that array.
 */
void update_candidates(pixel_T *candidates, int *num_candidates,
                       int *col_values, int c, mxLogical *bw, int M)
{
    mxLogical *pr;
    int r;
    
    /*
     * Update the max_col_values array.
     */
    pr = bw + M*c;
    for (r = 0; r < M; r++)
    {
        if (pr[r] != 0)
        {
            col_values[r] = c;
        }
    }

    /*
     * Update the candidates array.
     */
    *num_candidates = 0;
    for (r = 0; r < M; r++)
    {
        if (col_values[r] != FAKE_COL_VAL)
        {
            candidates[*num_candidates].x = col_values[r];
            candidates[*num_candidates].y = r;
            (*num_candidates)++;
        }
    }
}

/*
 * This is the algorithm "Create L from Candidates," bottom of column 2,
 * page 532 to top of column 1, page 533.  The following changes have
 * been made to the printed version of the algorithm:
 *   -  Script "L" replaced by "ell"
 *   -  "r" replaced by "c"
 *   -  "l" replaced by "p" (because SLE never names variables "l"!)
 *   -  "c" replaced by "num_candidates"
 *   -  fixed typo in paper:  "RemoveL([k-1],w,r)" should have been
 *      "Remove(L[k-1],L[k],w,r)."
 *
 * Outputs are Voronoi cell array ell and num_cells.  Inputs are
 * candidates, num_candidates, and c, which is the column currently
 * being processed.
 */
void compute_ell(pixel_T *ell, int *num_cells, pixel_T *candidates,
                 int num_candidates, int c)
{
    int k;
    int p;
    pixel_T w;

    if (num_candidates == 0)
    {
        *num_cells = 0;
    }
    else if (num_candidates == 1)
    {
        *num_cells = 1;
        ell[0] = candidates[0];
    }
    else
    {
        ell[0] = candidates[0];
        ell[1] = candidates[1];
        k = 1;
        p = 2;
        while (p < num_candidates)
        {
            w = candidates[p];
            while ((k >= 1) && remove_candidate(ell[k-1], ell[k], w, c))
            {
                k--;
            }
            k++;
            p++;
            ell[k] = w;
        }
        *num_cells = k+1;
    }
}

/*
 * Scan column c.  For each location in column c, compute the squared
 * distance to the nearest Voronoi cell center.  If that squared distance
 * is smaller than what's already in the array Out at that location,
 * then replace it.
 */
void assign_distances(double *Out, int M, int c, pixel_T *ell, int num_cells,
                      double *nn)
{
    int r;
    int current_cell = 0;
    double *pr;
    double *pr_nn;
    double sq_dist;
    double temp_sq_dist;
    double dx;
    double dy;
    int k;
    bool do_labels = nn != NULL;

    if (num_cells > 0)
    {
        pr = Out + M*c;
        if (do_labels)
        {
            pr_nn = nn + M*c;
        }
        for (r = 0; r < M; r++)
        {
            /*
             * Compute squared distance to current cell.
             */
            dx = c - ell[current_cell].x;
            dy = r - ell[current_cell].y;
            sq_dist = dx*dx + dy*dy;
            
            /*
             * See if we need to update the current cell.
             */
            for (k = current_cell+1; k < num_cells; k++)
            {
                dx = c - ell[k].x;
                dy = r - ell[k].y;
                temp_sq_dist = dx*dx + dy*dy;
                if (temp_sq_dist < sq_dist)
                {
                    current_cell = k;
                    sq_dist = temp_sq_dist;
                }
                else
                {
                    break;
                }
            }
            
            if (sq_dist < pr[r])
            {
                pr[r] = sq_dist;
                if (do_labels)
                {
                    pr_nn[r] = M * ell[current_cell].x + 
                        ell[current_cell].y + 1.0;
                }
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxLogical *bw;
    mxArray *Out;
    mxArray *NN_map;
    int M;
    int N;
    int *col_values;
    pixel_T *candidates;
    pixel_T *ell;
    int num_candidates;
    int num_cells;
    double *pr;
    double *nn = NULL;
    int k;
    int c;
    double inf = mxGetInf();

    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("Images:eucdist2:wrongNumInputs",
                          "EUCDIST2 requires one input argument.");
    }
    bw = mxGetLogicals(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);

    Out = mxCreateDoubleMatrix(M, N, mxREAL);
    plhs[0] = Out;
    pr = (double *) mxGetData(Out);
    for (k = 0; k < M*N; k++)
    {
        pr[k] = inf;
    }

    if (nlhs > 1)
    {
        /*
         * Return the nearest neighbor map.
         */
        NN_map = mxCreateDoubleMatrix(M, N, mxREAL);
        plhs[1] = NN_map;
        nn = (double *) mxGetData(NN_map);
    }

    candidates = (pixel_T *) mxMalloc(M * sizeof(*candidates));
    ell = (pixel_T *) mxMalloc(M * sizeof(*ell));
    col_values = (int *) mxMalloc(M * sizeof(*col_values));

    /*
     * Initialize the col_values array
     */
    for (k = 0; k < M; k++)
    {
        col_values[k] = FAKE_COL_VAL;
    }
    
    /*
     * First pass
     */
    for (c = 0; c < N; c++)
    {
        update_candidates(candidates, &num_candidates,
                          col_values, c, bw, M);
        compute_ell(ell, &num_cells, candidates, num_candidates, c);
        assign_distances(mxGetPr(Out), M, c, ell, num_cells, nn);
    }

    /*
     * Reinitialize the col_values array
     */
    for (k = 0; k < M; k++)
    {
        col_values[k] = FAKE_COL_VAL;
    }
    
    /*
     * Second pass
     */
    for (c = N-1; c >= 0; c--)
    {
        update_candidates(candidates, &num_candidates,
                          col_values, c, bw, M);
        compute_ell(ell, &num_cells, candidates, num_candidates, c);
        assign_distances(pr, M, c, ell, num_cells, nn);
    }

    /*
     * Return distance instead of squared distance.
     */
    for (k = 0; k < M*N; k++)
    {
        pr[k] = sqrt(pr[k]);
    }
}

#define _USE_MATH_DEFINES
#include <queue>
#include <vector>
#include <cmath>
#include <mex.h>
#include <iostream>
#define EXPORT_FCNS
#include "helper.h"

#include "libsrncp.h"

using namespace std;
#include <tuple>

//mex -v libsrncp.cpp

void _main();

constexpr double sqr(double x)
{
	return x * x;
}

static int verbose = 0;

FILE *logfile;

// For Matlab MEX compilation
#ifdef __cplusplus
extern "C" {
#endif
    
/**
 * \brief Unwrap and return phase of an adjacent voxel (v2) if the reference voxel (v1) is already unwrapped
 * \param v1 reference voxel phase value
 * \param v2 adjacent voxel phase value
 */
double PhaseUnwrap2px(double v1, double v2)
// delta = v2 - v1
// for positive delta:
// -> calculate n = integer number of 2pi intervals in delta interval
// -> if delta - 2pi*n <= pi, correction = -2pi*n
// -> if delta - 2pi*n > pi, correction = -2*pi - 2pi*n = -2pi*(n+1)
// return v2 + correction
// for negative v2 - v1:
// do the same with delta = -(v2 - v1) and return v2 - correction
{
    if (v2 > v1) 
        return v2 + (floor(-(v2 - v1 + M_PI) / (2 * M_PI)) + 1) * (2 * M_PI);
    else
        return v2 - (floor(-(v1 - v2 + M_PI) / (2 * M_PI)) + 1) * (2 * M_PI);
};

/**
 * \brief This function takes an angle in the range [-3*pi, 3*pi] and wraps it to the range [-pi, pi].
 * \param angle - input angle
 */
const inline double gamma(const double angle)
// This function takes an angle in the range [-3*pi, 3*pi] and
// wraps it to the range [-pi, pi].
{
    if(angle > M_PI)
        return angle - M_PI*2;
    else if(angle < -M_PI)
        return angle + M_PI*2;
    else 
        return angle;
}

/**
* \brief Fast 3D phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path (3D-SRNCP)
* Arrays phase and mask should have column-major order, as in Matlab
* i.e. to the ith row, jth column, kth plane in array corresponds the linear index [i+j*h+k*h*w]
* \param h height (number of rows) of phase and mask
* \param w width (number of columns) of phase and mask
* \param d depth (number of 2d planes) of phase and mask
* \param phase pointer to the array that stores the wrapped phase image and that will store the unwrapped image after processing
* \param mask pointer to the array that pointer to an array that contains 1 for valid and 0 for invalid voxels
* \param mask_largest_unwrapped_group a flag variable, if set to 1, then the 'mask' array will be set to 1 for the voxels in the biggest group of contiguous unwrapped voxels and to 0 elsewhere
*/
EXPORTED_FUNCTION void Unwrap3D(int h, int w, int d, double* phase, int* mask, int mask_largest_unwrapped_group)
{
    int i, j, k, m, found_nonzero_el, imh = h, imw = w, imd = d; // save original dimensions for IM
    int minrow = 0, maxrow = 0, mincol = 0, maxcol = 0, minslc = 0, maxslc = 0;
                
    // Create a second pointer to the arrays in the function input
    // to be able to shift it
    double* IM_phase = phase;
    int* IM_mask = mask;
  
    // 'Crop' arrays 'IM_phase' and 'IM_mask'
    // Find the non-zero voxel (1) in the 'IM_mask' array that has minimal row, column and plane indexes
    // Find the non-zero voxel (2) in the 'IM_mask' array that has maximal row, column and plane indexes
	if (IM_mask != nullptr)
	{
        // If mask is specified, set phase of voxels outside it to 0
		for (i = 0; i < h; i++)
			for (j = 0; j < w; j++)
				for (k = 0; k < d; k++)
                    if (IM_mask[i + j * h + k * h * w] == 0) IM_phase[i + j * h + k * h * w] = 0;

		found_nonzero_el = 0; //to break a triple loop
		for (i = 0; (i < h) && (!found_nonzero_el); i++)
			for (j = 0; (j < w) && (!found_nonzero_el); j++)
				for (k = 0; (k < d) && (!found_nonzero_el); k++)
					if (IM_mask[i + j * h + k * h * w] != 0)
					{
						minrow = i;
						found_nonzero_el = 1;
					}
		if (!found_nonzero_el) //mask is empty            
			minrow, maxrow, mincol, maxcol, minslc, maxslc = 0;
		else
		{
			found_nonzero_el = 0;
			for (j = 0; (j < w) && (!found_nonzero_el); j++)
				for (i = 0; (i < h) && (!found_nonzero_el); i++)
					for (k = 0; (k < d) && (!found_nonzero_el); k++)
						if (IM_mask[i + j * h + k * h * w] != 0)
						{
							mincol = j;
							found_nonzero_el = 1;
						}
			found_nonzero_el = 0;
			for (k = 0; (k < d) && (!found_nonzero_el); k++)
				for (i = 0; (i < h) && (!found_nonzero_el); i++)
					for (j = 0; (j < w) && (!found_nonzero_el); j++)
						if (IM_mask[i + j * h + k * h * w] != 0)
						{
							minslc = k;
							found_nonzero_el = 1;
						}
			found_nonzero_el = 0;
			for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
				for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
					for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
						if (IM_mask[i + j * h + k * h * w] != 0)
						{
							maxrow = i;
							found_nonzero_el = 1;
						}
			found_nonzero_el = 0;
			for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
				for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
					for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
						if (IM_mask[i + j * h + k * h * w] != 0)
						{
							maxcol = j;
							found_nonzero_el = 1;
						}
			found_nonzero_el = 0;
			for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
				for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
					for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
						if (IM_mask[i + j * h + k * h * w] != 0)
						{
							maxslc = k;
							found_nonzero_el = 1;
						}

			if (minrow < 0) minrow = 0;
			if (mincol < 0) mincol = 0;
			if (minslc < 0) minslc = 0;
			if (maxrow > h - 1) maxrow = h - 1;
			if (maxcol > w - 1) maxcol = w - 1;
			if (maxslc > d - 1) maxslc = d - 1;
           
            // In the 'IM_phase' and 'IM_mask' arrays, shift the begin pointer to the 1st voxel
			IM_phase += minrow + mincol * imh + minslc * imw * imh;
			IM_mask += minrow + mincol * imh + minslc * imw * imh;

			h = maxrow - minrow + 1; // Cropping: define new height
			w = maxcol - mincol + 1; // Cropping: define new width
			d = maxslc - minslc + 1; // Cropping: define new depth
		};
	};

    // the voxel quality map
    vector <double> vqmap(h * w * d);
    int inmask;
    int ci,cj,ck;
    double val;
    if (verbose == 1) logfile = fopen("3dsrncp.log","w");

    //srand(1); // FOR DEBUGGING ONLY!

    // Calculate voxel quality (second differences) map 
    vector<int> sdflag(3*3*3);
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (k = 0; k < d; k++)
            {
                // set the quality by default to 10000 (means not calculated)                
                // maximal theoretical value for 3D case is (2pi)^2 * 13 = 513.219
                vqmap[i + j * h + k * h * w] = 10000;
                
                //do not calculate qualities of voxels for which not all 9 neighbours are inside array boundaries or in the mask
                inmask = 1;
                if ((i == 0) || (i == h-1) || (j == 0) || (j == w-1) || (k == 0) || (k == d-1)) inmask = 0;
                else if (IM_mask != nullptr)
                    for (ci = -1; ci < 2; ci++)
                        for (cj = -1; cj < 2; cj++)
                            for (ck = -1; ck < 2; ck++)
                                if (IM_mask[(i+ci)+(j+cj)*imh+(k+ck)*imh*imw] == 0) inmask = 0;
                if (inmask == 0) continue;

                // calculate and sum up horizontal, vertical, normal and 10 diagonal second differences
                val = 0;                
                fill(sdflag.begin(), sdflag.end(), 0);
                for (ci = -1; ci < 2; ci++)
                    for (cj = -1; cj < 2; cj++)
                        for (ck = -1; ck < 2; ck++)
                        {
                            if (sdflag[(1+ci) + (1+cj)*3 + (1+ck)*3*3] == 1) continue;
                            sdflag[(1-ci) + (1-cj)*3 + (1-ck)*3*3] = 1;
                            sdflag[(1+ci) + (1+cj)*3 + (1+ck)*3*3] = 1;
                            val += sqr(gamma(IM_phase[(i-ci)+(j-cj)*imh+(k-ck)*imh*imw]-IM_phase[i+j*imh+k*imh*imw]) -
                                    gamma(IM_phase[i+j*imh+k*imh*imw]-IM_phase[(i+ci)+(j+cj)*imh+(k+ck)*imh*imw]));
                        }
                
                vqmap[i + j * h + k * h * w] = val;
            }

    typedef tuple<int, int, int, int, double> Edge;    
    // the list of edges that we will sort later
    vector <Edge> EdgeList(h * w * d * 3);    

    // Reference tables used for grouping voxels together
    // the array that stores the index of the first voxel in the group to which belongs the indexed voxel
    vector <int> VoxelGroupFirst(h * w * d);    
    // the array that stores the index of the next voxel in the group to which belongs the indexed voxel
    vector <int> VoxelGroupNext(h * w * d);
    // the array that stores the index of the last voxel in the group to which belongs the indexed voxel
    // the indexed voxel should be first in its group, otherwise reference is not valid
    vector <int> VoxelGroupLast(h * w * d);
    // the array that stores the size of the group to which belongs the indexed voxel
    // the indexed voxel should be first in its group, otherwise reference is not valid
    vector <int> VoxelGroupSize(h * w * d);
    
    // Calculate edge qualities and initialize voxel group index matrix VoxelGroupInd to -1 (set voxel as non-unwrapped)
    int valid_ex, valid_ey, valid_ez; // flags to check if we can calculate x, y, z edge quality
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (k = 0; k < d; k++)
            {
                // each voxel initially belongs to its own group where it is alone
                // so the first voxel in the group is the indexed voxel (i,j,k)
                VoxelGroupFirst[i + j * h + k * h * w] = i + j * h + k * h * w;
                // .. and because it is alone in its group, it has no next voxel,
                VoxelGroupNext[i + j * h + k * h * w] = -1;
                // .. and the last voxel in the group is the same as the first,
                VoxelGroupLast[i + j * h + k * h * w] = i + j * h + k * h * w;
                // .. and the size of group is 1
                VoxelGroupSize[i + j * h + k * h * w] = 1;
                
                // are both voxels on x edge are 'valid' (inside valid matrix boundaries and, if given, IM_mask)?
                valid_ex = 1;
                // if the second voxel on the x edge (i+1,j,k) is outside matrix, x edge is 'invalid'
                if (i == h-1) valid_ex = 0;
                else if (IM_mask != nullptr)
                    // if either of the voxels in the edge are not in the mask IM_mask, x edge is also 'invalid'
                    if ((IM_mask[i+j*imh+k*imh*imw] == 0)||(IM_mask[(i+1)+j*imh+k*imh*imw] == 0)) valid_ex = 0;
                // same reasoning for y edge
                valid_ey = 1;
                if (j == w-1) valid_ey = 0;
                else if (IM_mask != nullptr)
                    if ((IM_mask[i+j*imh+k*imh*imw] == 0)||(IM_mask[i+(j+1)*imh+k*imh*imw] == 0)) valid_ey = 0;
                // same reasoning for z edge
                valid_ez = 1;
                if (k == d-1) valid_ez = 0;
                else if (IM_mask != nullptr)
                    if ((IM_mask[i+j*imh+k*imh*imw] == 0)||(IM_mask[i+j*imh+(k+1)*imh*imw] == 0)) valid_ez = 0;
                
                // if both voxels in the x edge are valid
                if (valid_ex == 1)
                {
                    val = vqmap[i + j * h + k * h * w] + vqmap[(i+1) + j * h + k * h * w];
                    // if both voxels are on the border and have the default quality of 10000
                    if ((vqmap[i + j * h + k * h * w] > 9999)&&(vqmap[(i+1) + j * h + k * h * w] > 9999))
                        // add some random value
                        val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                }
                // if one of voxels are not valid, then the edge should not be processed
                // set a high value to place it in the end of the list after sorting
                else val = 1e11;
                // store the value in the list
                EdgeList[i + j * h + k * h * w + 0 * h * w * d] = {i,j,k,0,val};

                // same for y edge
                if (valid_ey == 1)
                {
                    val = vqmap[i + j * h + k * h * w] + vqmap[i + (j+1) * h + k * h * w];
                    if ((vqmap[i + j * h + k * h * w] > 9999)&&(vqmap[i + (j+1) * h + k * h * w] > 9999))
                        val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                }
                else val = 1e11;
                EdgeList[i + j * h + k * h * w + 1 * h * w * d] = {i,j,k,1,val};
                
                // same for z edge
                if (valid_ez == 1)
                {
                    val = vqmap[i + j * h + k * h * w] + vqmap[i + j * h + (k+1) * h * w];
                    if ((vqmap[i + j * h + k * h * w] > 9999)&&(vqmap[i + j * h + (k+1) * h * w] > 9999))
                        val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                }
                else val = 1e11;
                EdgeList[i + j * h + k * h * w + 2 * h * w * d] = {i,j,k,2,val};
            }
    // we do not need voxel quality map anymore
    vqmap.clear();

    // sort the edge list in the order of descending quality (ascending value)
    sort(EdgeList.begin(), EdgeList.end(), [](const Edge &a, const Edge &b) {return (get<4>(a) < get<4>(b));});

    int t; // Edge type (0=x, 1=y, 2=z)
    int g_i, g_ci; // Group index for the first voxel in the edge (i,j,k) and second voxel (ci,cj,ck);
    int n,vi,vj,vk,vci,vcj,vck;
    
    int g_max = 0; // linearized index of the first voxel in the group with maximal number of voxels;
    
    if (verbose == 1) for (m = 0; m < h * w * d * 3; m++)
    {
        if (get<4>(EdgeList[m]) > 1e10) break;
        //i, j, k: coordinates of first voxel in the edge
        i = get<0>(EdgeList[m]);
        j = get<1>(EdgeList[m]);
        k = get<2>(EdgeList[m]);
        t = get<3>(EdgeList[m]);
        // ci, cj, ck: coordinates of second voxel in the edge
        ci = i; cj = j; ck = k;
        if (t == 0) ci++;
        else if (t == 1) cj++;
        else ck++;
        if (verbose == 1) fprintf(logfile,"edge %d: %d %d %d - %d %d %d val %.8e\n",m,i,j,k,ci,cj,ck,get<4>(EdgeList[m]));
    }
    
    for (m = 0; m < h * w * d * 3; m++)
    {
        if (get<4>(EdgeList[m]) > 1e10) break;
        //i, j, k: coordinates of first voxel in the edge
        i = get<0>(EdgeList[m]);
        j = get<1>(EdgeList[m]);
        k = get<2>(EdgeList[m]);
        t = get<3>(EdgeList[m]);
        // ci, cj, ck: coordinates of second voxel in the edge
        ci = i; cj = j; ck = k;
        if (t == 0) ci++;
        else if (t == 1) cj++;
        else ck++;
        //g_i, g_ci: linearized column-major order indexes of (i,j,k) and (ci,cj,ck)
        g_i = i + j * h + k * h * w;
        g_ci = ci + cj * h + ck * h * w;
        
        // decode indexes of the first voxel in the group of first voxel in the edge
        vi = (VoxelGroupFirst[g_i] % (h * w)) % h;
        vj = (int)(VoxelGroupFirst[g_i] % (h * w)) / h;
        vk = (int)VoxelGroupFirst[g_i] / (h * w);
        // decode indexes of the first voxel in the group of second voxel in the edge
        vci = (VoxelGroupFirst[g_ci] % (h * w)) % h;
        vcj = (int)(VoxelGroupFirst[g_ci] % (h * w)) / h;
        vck = (int)VoxelGroupFirst[g_ci] / (h * w);            
        if (verbose == 1) fprintf(logfile,"proc %d: %d %d %d (1st %d %d %d) - %d %d %d (1st %d %d %d) val %.8e\n",m,i,j,k,vi,vj,vk,ci,cj,ck,vci,vcj,vck,get<4>(EdgeList[m]));        
        
        // If both (i,j,k) and (ci, cj, ck) are already in the same group, nothing to do
        // This edge will not be used
        if (VoxelGroupFirst[g_i] == VoxelGroupFirst[g_ci]) continue;
        
        // if (ci, cj, ck) voxel have not been unwrapped before, 
        // then unwrap (ci, cj, ck) with respect to (i, j, k) and add it to its group;  
        if (VoxelGroupSize[VoxelGroupFirst[g_ci]] == 1)
        {
            // PhaseUnwrap2px(x,y) returns the unwrapped value of y with respect to x;
            // Unwrap (ci,cj,ck) with respect to (i,j,k)
            IM_phase[ci + cj * imh + ck * imh * imw] = PhaseUnwrap2px(IM_phase[i + j * imh + k * imh * imw], IM_phase[ci + cj * imh + ck * imh * imw]);
            if (VoxelGroupSize[VoxelGroupFirst[g_i]] == 1)
                {if (verbose == 1) fprintf(logfile,"ref. %d %d %d unwrap %d %d %d\n", i,j,k, ci,cj,ck);}
            else
                if (verbose == 1) fprintf(logfile,"ref. %d %d %d (group of %d el.) unwrap %d %d %d\n",i,j,k, VoxelGroupSize[VoxelGroupFirst[g_i]], ci,cj,ck);
            // Add voxel (ci,cj,ck) to the group of (i,j,k), update reference tables
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_i]]] = g_ci;
            VoxelGroupLast[VoxelGroupFirst[g_i]] = g_ci;
            VoxelGroupSize[VoxelGroupFirst[g_i]]++;
            VoxelGroupFirst[g_ci] = VoxelGroupFirst[g_i];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_i];
            continue;
        }

        // if (i, j, k) voxel have not been unwrapped before, 
        // then unwrap (i, j, k) with respect to (ci, cj, ck) and add it to its group;  
        if (VoxelGroupSize[VoxelGroupFirst[g_i]] == 1)
        {
            // Unwrap (i,j,k) with respect to (ci,cj,ck)
            IM_phase[i + j * imh + k * imh * imw] = PhaseUnwrap2px(IM_phase[ci + cj * imh + ck * imh * imw], IM_phase[i + j * imh + k * imh * imw]);
            if (verbose == 1) fprintf(logfile,"unwrap %d %d %d ref. %d %d %d (group of %d el.)\n", i,j,k, ci,cj,ck, VoxelGroupSize[VoxelGroupFirst[g_ci]]);
            // Add voxel (i,j,k) to the group of (ci,cj,ck), update reference tables
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_ci]]] = g_i;
            VoxelGroupLast[VoxelGroupFirst[g_ci]] = g_i;
            VoxelGroupSize[VoxelGroupFirst[g_ci]]++;
            VoxelGroupFirst[g_i] = VoxelGroupFirst[g_ci];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_ci]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_ci];
            continue;
        }        
        
        // if both (i, j, k) and (ci, cj, ck) have been unwrapped before and belong to different groups,
        // unwrap and add group of (ci, cj, ck) to the group of (i, j, k) if the group of (i, j, k) is bigger
        if (VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[VoxelGroupFirst[g_ci]])
        {
            if (verbose == 1) fprintf(logfile,"add to group of %d %d %d (%d el., other %d el.):\n",i,j,k,VoxelGroupSize[VoxelGroupFirst[g_i]],VoxelGroupSize[VoxelGroupFirst[g_ci]]);
            // when we unwrap the edge linking two groups, the value of the voxel in the second group is shifted by 'val'
            // we save this 'val' value to subsequently add it to all voxels in the second group (of voxel (ci, cj, ck)), 
            // thus not breaking the relative phase differences between voxels in that group
            val = PhaseUnwrap2px(IM_phase[i + j * imh + k * imh * imw], IM_phase[ci + cj * imh + ck * imh * imw]) - IM_phase[ci + cj * imh + ck * imh * imw];
            // Joining the group of (ci, cj, ck) with the group of (i, j, k)
            // Link the last voxel in the group of (i, j, k) to the first voxel in the group of (ci, cj, ck)
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_i]]] = VoxelGroupFirst[g_ci];
            // Update the last voxel in the group of (i, j, k): now it is the last voxel in the group of (ci, cj, ck)
            VoxelGroupLast[VoxelGroupFirst[g_i]] = VoxelGroupLast[VoxelGroupFirst[g_ci]];
            // Update the size of the group of voxel (i, j, k)  
            VoxelGroupSize[VoxelGroupFirst[g_i]] += VoxelGroupSize[VoxelGroupFirst[g_ci]];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_i];
            // Index of the first voxel in the group of (ci, cj, ck) 
            n = VoxelGroupFirst[g_ci];
            while (n != -1)
            {
                // Get ci,cj,ck coordinates of the current voxel in the group of (ci, cj, ck) 
                ci = (n % (h * w)) % h;
                cj = (int)(n % (h * w)) / h;
                ck = (int)n / (h * w);
                if (verbose == 1) fprintf(logfile," + %d %d %d:\n",ci,cj,ck);
                // Update the reference to the first voxel in the group
                // as both groups are merged to the group of (i, j, k) voxel
                VoxelGroupFirst[n] = VoxelGroupFirst[g_i];
                // Add the 'val' offset to the phase value of the current voxel in the group of (ci, cj, ck) 
                IM_phase[ci + cj * imh + ck * imh * imw] += val;
                // Move to next voxel in the group of (ci, cj, ck) 
                n = VoxelGroupNext[n];
            }
            continue;
        }
        // unwrap and add group of (i, j, k) to the group of (ci, cj, ck) if the group of (ci, cj, ck) is bigger
        if (verbose == 1) fprintf(logfile,"add to group of %d %d %d (%d el., other %d el.):\n",ci,cj,ck,VoxelGroupSize[VoxelGroupFirst[g_ci]],VoxelGroupSize[VoxelGroupFirst[g_i]]);
        // when we unwrap the edge linking two groups, the value of the voxel in the second group is shifted by 'val'
        // we save this 'val' value to subsequently add it to all voxels in the second group (of voxel (i, j, k)),
        // thus not breaking the relative phase differences between voxels in that group
        val = PhaseUnwrap2px(IM_phase[ci + cj * imh + ck * imh * imw], IM_phase[i + j * imh + k * imh * imw]) - IM_phase[i + j * imh + k * imh * imw];
        // Joining the group of (i, j, k) with the group of (ci, cj, ck)
        // Link the last voxel in the group of (ci, cj, ck) to the first voxel in the group of (i, j, k)
        VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_ci]]] = VoxelGroupFirst[g_i];
        // Update the last voxel in the group of (ci, cj, ck): now it is the last voxel in the group of (i, j, k)
        VoxelGroupLast[VoxelGroupFirst[g_ci]] = VoxelGroupLast[VoxelGroupFirst[g_i]];
        // Update the size of the group of voxel (ci, cj, ck)  
        VoxelGroupSize[VoxelGroupFirst[g_ci]] += VoxelGroupSize[VoxelGroupFirst[g_i]];
        // Update eventually the reference g_max to the largest group
        if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_ci]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_ci];
        // Index of the first voxel in the group of (i, j, k) 
        n = VoxelGroupFirst[g_i];
        while (n != -1)
        {
            // Get i,j,k coordinates of the current voxel in the group of (i, j, k) 
            i = (n % (h * w)) % h;
            j = (int)(n % (h * w)) / h;
            k = (int)n / (h * w);
            if (verbose == 1) fprintf(logfile," + %d %d %d:\n",i,j,k);
            // Update the reference to the first voxel in the group
            // as both groups are merged to the group of (ci, cj, ck) voxel
            VoxelGroupFirst[n] = VoxelGroupFirst[g_ci];
            // Add the 'val' offset to the phase value of the current voxel in the group of (i, j, k) 
            IM_phase[i + j * imh + k * imh * imw] += val;
            // Move to next voxel in the group of (i, j, k) 
            n = VoxelGroupNext[n];
        }
    }
    // If the flag is set to 1, mark in 'mask' as 1 all voxels in the biggest group of unwrapped voxels
    if ((mask_largest_unwrapped_group == 1)&&(IM_mask != nullptr))
    {
        // Reinitialize mask with 0
        fill(mask, mask + imh * imw * imd, 0);
        // Go to the first voxel in the biggest group
        n = VoxelGroupFirst[g_max];
        while (n != -1)
        {
           // Get i,j,k coordinates of the current voxel in the biggest group
            i = (n % (h * w)) % h;
            j = (int)(n % (h * w)) / h;
            k = (int)n / (h * w);
            // Set 'mask' to 1
            IM_mask[i + j * imh + k * imh * imw] = 1;
            // Move to next voxel in the biggest group
            n = VoxelGroupNext[n];
        }   
    }
    
    if (verbose == 1) fclose(logfile);
}

/**
* \brief Fast 4D phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path (4D-SRNCP)
* Arrays phase and mask should have column-major order, as in Matlab
* i.e. to the ith row, jth column, kth plane and pth slab in array corresponds the linear index [i+j*h+k*h*w+p*h*w*d]
* \param h height (number of rows) of phase and mask
* \param w width (number of columns) of phase and mask
* \param d depth (number of 2d planes) of phase and mask
* \param s spissitude (number of slabs) of phase and mask
* \param phase pointer to the array that stores the wrapped phase image and that will store the unwrapped image after processing
* \param mask pointer to the array that contains 1 for valid and 0 for invalid voxels
* \param mask_largest_unwrapped_group a flag variable, if set to 1, then the 'mask' array will be set to 1 for the voxels in the biggest group of contiguous unwrapped voxels and to 0 elsewhere
*/

EXPORTED_FUNCTION void Unwrap4D(int h, int w, int d, int s, double* phase, int* mask, int mask_largest_unwrapped_group)
{
    int i, j, k, p, m, found_nonzero_el, imh = h, imw = w, imd = d, ims = s; // save original dimensions for IM
    int minrow = 0, maxrow = 0, mincol = 0, maxcol = 0, minslc = 0, maxslc = 0, minvol = 0, maxvol = 0;
                
    // Create a second pointer to the arrays in the function input
    // to be able to shift it
    double* IM_phase = phase;
    int* IM_mask = mask;
  
    // 'Crop' arrays 'IM_phase' and 'IM_mask'
    // Find the non-zero voxel (1) in the 'IM_mask' array that has minimal row, column and plane indexes
    // Find the non-zero voxel (2) in the 'IM_mask' array that has maximal row, column and plane indexes
	if (IM_mask != nullptr)
	{
        // If mask is specified, set phase of voxels outside it to 0
		for (i = 0; i < h; i++)
			for (j = 0; j < w; j++)
				for (k = 0; k < d; k++)
					for (p = 0; p < s; p++)
						if (IM_mask[i + j * h + k * h * w + p * h * w * d] == 0)
                            IM_phase[i + j * h + k * h * w + p * h * w * d] = 0;

        found_nonzero_el = 0; //to break a quadruple loop
		for (i = 0; (i < h) && (!found_nonzero_el); i++)
			for (j = 0; (j < w) && (!found_nonzero_el); j++)
				for (k = 0; (k < d) && (!found_nonzero_el); k++)
					for (p = 0; (p < s) && (!found_nonzero_el); p++)
						if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
						{
							minrow = i;
							found_nonzero_el = 1;
						}
		if (!found_nonzero_el) //mask is empty            
			minrow, maxrow, mincol, maxcol, minslc, maxslc, minvol, maxvol = 0;
		else
		{
			found_nonzero_el = 0;
			for (j = 0; (j < w) && (!found_nonzero_el); j++)
				for (i = 0; (i < h) && (!found_nonzero_el); i++)
					for (k = 0; (k < d) && (!found_nonzero_el); k++)
						for (p = 0; (p < s) && (!found_nonzero_el); p++)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								mincol = j;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (k = 0; (k < d) && (!found_nonzero_el); k++)
				for (i = 0; (i < h) && (!found_nonzero_el); i++)
					for (j = 0; (j < w) && (!found_nonzero_el); j++)
						for (p = 0; (p < s) && (!found_nonzero_el); p++)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								minslc = k;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (p = 0; (p < s) && (!found_nonzero_el); p++)
				for (i = 0; (i < h) && (!found_nonzero_el); i++)
					for (j = 0; (j < w) && (!found_nonzero_el); j++)
						for (k = 0; (k < d) && (!found_nonzero_el); k++)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								minvol = p;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
				for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
					for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
						for (p = s - 1; (p > -1) && (!found_nonzero_el); p--)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								maxrow = i;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
				for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
					for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
						for (p = s - 1; (p > -1) && (!found_nonzero_el); p--)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								maxcol = j;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
				for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
					for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
						for (p = s - 1; (p > -1) && (!found_nonzero_el); p--)
							if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
							{
								maxslc = k;
								found_nonzero_el = 1;
							}
			found_nonzero_el = 0;
			for (p = s - 1; (p > -1) && (!found_nonzero_el); p--)
				for (i = h - 1; (i > -1) && (!found_nonzero_el); i--)
					for (j = w - 1; (j > -1) && (!found_nonzero_el); j--)
                        for (k = d - 1; (k > -1) && (!found_nonzero_el); k--)
                            if (IM_mask[i + j * h + k * h * w + p * h * w * d] != 0)
                            {
                                maxvol = p;
                                found_nonzero_el = 1;
                            }
            
            if (minrow < 0) minrow = 0;
            if (mincol < 0) mincol = 0;
            if (minslc < 0) minslc = 0;
            if (minvol < 0) minvol = 0;
            if (maxrow > h - 1) maxrow = h - 1;
            if (maxcol > w - 1) maxcol = w - 1;
            if (maxslc > d - 1) maxslc = d - 1;
            if (maxvol > s - 1) maxvol = s - 1;

            // In the 'IM_phase' and 'IM_mask' arrays, shift the begin pointer to the 1st voxel
			IM_phase += minrow + mincol * imh + minslc * imw * imh + minvol * imd * imw * imh;
			IM_mask += minrow + mincol * imh + minslc * imw * imh + minvol * imd * imw * imh;

			h = maxrow - minrow + 1; // Cropping: define new height
			w = maxcol - mincol + 1; // Cropping: define new width
			d = maxslc - minslc + 1; // Cropping: define new depth
            s = maxvol - minvol + 1; // Cropping: define new spissitude
		};
	};

    // the voxel quality map
    vector <double> vqmap(h * w * d * s);
    int inmask;
    int ci,cj,ck,cp;
    double val;
    if (verbose == 1) logfile = fopen("3dsrncp.log","w");

    //srand(1); // FOR DEBUGGING ONLY!

    // Calculate voxel quality (second differences) map 
    vector<int> sdflag(3*3*3*3);
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (k = 0; k < d; k++)
                for (p = 0; p < s; p++)
                {
                    // set the quality by default to 10000 (means not calculated)
                    // maximal theoretical value for 4D case is (2pi)^2 * 40 = 1579.14
                    vqmap[i + j * h + k * h * w + p * h * w * d] = 10000;
                    
                    //do not calculate qualities of voxels for which not all 9 neighbours are inside array boundaries or in the mask
                    inmask = 1;
                    if ((i == 0) || (i == h-1) || (j == 0) || (j == w-1) || (k == 0) || (k == d-1) || (p == 0) || (p == s-1)) inmask = 0;
                    else if (IM_mask != nullptr)
                        for (ci = -1; ci < 2; ci++)
                            for (cj = -1; cj < 2; cj++)
                                for (ck = -1; ck < 2; ck++)
                                    for (cp = -1; cp < 2; cp++)
                                        if (IM_mask[(i+ci)+(j+cj)*imh+(k+ck)*imh*imw+(p+cp)*imh*imw*imd] == 0) inmask = 0;
                    if (inmask == 0) continue;
                    
                    // calculate and sum up horizontal, vertical, normal and 10 diagonal second differences
                    val = 0;
                    fill(sdflag.begin(), sdflag.end(), 0);
                    for (ci = -1; ci < 2; ci++)
                        for (cj = -1; cj < 2; cj++)
                            for (ck = -1; ck < 2; ck++)
                                for (cp = -1; cp < 2; cp++)
                                {
                                    if (sdflag[(1+ci) + (1+cj)*3 + (1+ck)*3*3 + (1+cp)*3*3*3] == 1) continue;
                                    sdflag[(1-ci) + (1-cj)*3 + (1-ck)*3*3 + (1-cp)*3*3*3] = 1;
                                    sdflag[(1+ci) + (1+cj)*3 + (1+ck)*3*3 + (1+cp)*3*3*3] = 1;
                                    val += sqr(gamma(IM_phase[(i-ci)+(j-cj)*imh+(k-ck)*imh*imw+(p-cp)*imh*imw*imd]-IM_phase[i+j*imh+k*imh*imw+p*imh*imw*imd]) -
                                            gamma(IM_phase[i+j*imh+k*imh*imw+p*imh*imw*imd]-IM_phase[(i+ci)+(j+cj)*imh+(k+ck)*imh*imw+(p+cp)*imh*imw*imd]));
                                }
                    
                    vqmap[i + j * h + k * h * w + p * h * w * d] = val;
                }

    typedef tuple<int, int, int, int, int, double> Edge;    
    // the list of edges that we will sort later
    vector <Edge> EdgeList(h * w * d * s * 4);

    // Reference tables used for grouping voxels together
    // the array that stores the index of the first voxel in the group to which belongs the indexed voxel
    vector <int> VoxelGroupFirst(h * w * d * s); 
    // the array that stores the index of the next voxel in the group to which belongs the indexed voxel
    vector <int> VoxelGroupNext(h * w * d * s);
    // the array that stores the index of the last voxel in the group to which belongs the indexed voxel
    // the indexed voxel should be first in its group, otherwise reference is not valid
    vector <int> VoxelGroupLast(h * w * d * s);
    // the array that stores the size of the group to which belongs the indexed voxel
    // the indexed voxel should be first in its group, otherwise reference is not valid
    vector <int> VoxelGroupSize(h * w * d * s);
    
    // Calculate edge qualities and initialize voxel group index matrix VoxelGroupInd to -1 (set voxel as non-unwrapped)
    int valid_ex, valid_ey, valid_ez, valid_es; // flags to check if we can calculate x, y, z, s edge quality
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (k = 0; k < d; k++)
                for (p = 0; p < s; p++)
                {
                    // each voxel initially belongs to its own group where it is alone
                    // so the first voxel in the group is the indexed voxel (i,j,k)
                    VoxelGroupFirst[i + j * h + k * h * w + p * h * w * d] = i + j * h + k * h * w + p * h * w * d;
                    // .. and because it is alone in its group, it has no next voxel,
                    VoxelGroupNext[i + j * h + k * h * w + p * h * w * d] = -1;
                    // .. and the last voxel in the group is the same as the first,
                    VoxelGroupLast[i + j * h + k * h * w + p * h * w * d] = i + j * h + k * h * w + p * h * w * d;
                    // .. and the size of group is 1
                    VoxelGroupSize[i + j * h + k * h * w + p * h * w * d] = 1;
                    
                    // are both voxels on x edge are 'valid' (inside valid matrix boundaries and, if given, IM_mask)?
                    valid_ex = 1;
                    // if the second voxel on the x edge (i+1,j,k,p) is outside matrix, x edge is 'invalid'
                    if (i == h-1) valid_ex = 0;
                    else if (IM_mask != nullptr)
                        // if either of the voxels in the edge are not in the mask IM_mask, x edge is also 'invalid'
                        if ((IM_mask[i+j*imh+k*imh*imw+p*imh*imw*imd] == 0)||(IM_mask[(i+1)+j*imh+k*imh*imw+p*imh*imw*imd] == 0)) valid_ex = 0;
                    // same reasoning for y edge
                    valid_ey = 1;
                    if (j == w-1) valid_ey = 0;
                    else if (IM_mask != nullptr)
                        if ((IM_mask[i+j*imh+k*imh*imw+p*imh*imw*imd] == 0)||(IM_mask[i+(j+1)*imh+k*imh*imw+p*imh*imw*imd] == 0)) valid_ey = 0;
                    // same reasoning for z edge
                    valid_ez = 1;
                    if (k == d-1) valid_ez = 0;
                    else if (IM_mask != nullptr)
                        if ((IM_mask[i+j*imh+k*imh*imw+p*imh*imw*imd] == 0)||(IM_mask[i+j*imh+(k+1)*imh*imw+p*imh*imw*imd] == 0)) valid_ez = 0;
                    // same reasoning for s edge
                    valid_es = 1;
                    if (p == s-1) valid_es = 0;
                    else if (IM_mask != nullptr)
                        if ((IM_mask[i+j*imh+k*imh*imw+p*imh*imw*imd] == 0)||(IM_mask[i+j*imh+k*imh*imw+(p+1)*imh*imw*imd] == 0)) valid_es = 0;
                    
                    // if both voxels in the x edge are valid
                    if (valid_ex == 1)
                    {
                        val = vqmap[i + j * h + k * h * w + p * h * w * d] + vqmap[(i+1) + j * h + k * h * w + p * h * w * d];
                        // if both voxels are on the border and have the default quality of 10000
                        if ((vqmap[i + j * h + k * h * w + p * h * w * d] > 9999)&&(vqmap[(i+1) + j * h + k * h * w + p * h * w * d] > 9999))
                            // add some random value
                            val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                    }
                    // if one of voxels are not valid, then the edge should not be processed
                    // set a high value to place it in the end of the list after sorting
                    else val = 1e11;
                    // store the value in the list
                    EdgeList[i + j * h + k * h * w + p * h * w * d + 0 * h * w * d * s] = {i,j,k,p,0,val};
                    
                    // same for y edge
                    if (valid_ey == 1)
                    {
                        val = vqmap[i + j * h + k * h * w + p * h * w * d] + vqmap[i + (j+1) * h + k * h * w + p * h * w * d];
                        if ((vqmap[i + j * h + k * h * w + p * h * w * d] > 9999)&&(vqmap[i + (j+1) * h + k * h * w + p * h * w * d] > 9999))
                            val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                    }
                    else val = 1e11;
                    EdgeList[i + j * h + k * h * w + p * h * w * d + 1 * h * w * d * s] = {i,j,k,p,1,val};
                    
                    // same for z edge
                    if (valid_ez == 1)
                    {
                        val = vqmap[i + j * h + k * h * w + p * h * w * d] + vqmap[i + j * h + (k+1) * h * w + p * h * w * d];
                        if ((vqmap[i + j * h + k * h * w + p * h * w * d] > 9999)&&(vqmap[i + j * h + (k+1) * h * w + p * h * w * d] > 9999))
                            val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                    }
                    else val = 1e11;
                    EdgeList[i + j * h + k * h * w + p * h * w * d + 2 * h * w * d * s] = {i,j,k,p,2,val};
                    
                    // same for s edge
                    if (valid_es == 1)
                    {
                        val = vqmap[i + j * h + k * h * w + p * h * w * d] + vqmap[i + j * h + k * h * w + (p+1) * h * w * d];
                        if ((vqmap[i + j * h + k * h * w + p * h * w * d] > 9999)&&(vqmap[i + j * h + k * h * w + (p+1) * h * w * d] > 9999))
                            val += 1000*(double)rand()*(double)rand()/RAND_MAX/RAND_MAX;
                    }
                    else val = 1e11;
                    EdgeList[i + j * h + k * h * w + p * h * w * d + 3 * h * w * d * s] = {i,j,k,p,3,val};               }
    // we do not need voxel quality map anymore
    vqmap.clear();

    // sort the edge list in the order of descending quality (ascending value)
    sort(EdgeList.begin(), EdgeList.end(), [](const Edge &a, const Edge &b) {return (get<5>(a) < get<5>(b));});
    
    int t; // Edge type (0=x, 1=y, 2=z, 3=s)
    int g_i, g_ci; // Group index for the first voxel in the edge (i,j,k) and second voxel (ci,cj,ck);
    int n,vi,vj,vk,vp,vci,vcj,vck,vcp;
    
    int g_max = 0; // linearized index of the first voxel in the group with maximal number of voxels;
    
    if (verbose == 1) for (m = 0; m < h * w * d * s * 4; m++)
    {
        if (get<5>(EdgeList[m]) > 1e10) break;
        //i, j, k, p: coordinates of first voxel in the edge
        i = get<0>(EdgeList[m]);
        j = get<1>(EdgeList[m]);
        k = get<2>(EdgeList[m]);
        p = get<3>(EdgeList[m]);
        t = get<4>(EdgeList[m]);
        // ci, cj, ck: coordinates of second voxel in the edge
        ci = i; cj = j; ck = k; cp = p;
        if (t == 0) ci++;
        else if (t == 1) cj++;
        else if (t == 2) ck++;
        else cp++;
        if (verbose == 1) fprintf(logfile,"edge %d: %d %d %d %d - %d %d %d %d val %.8e\n",m,i,j,k,p,ci,cj,ck,cp,get<5>(EdgeList[m]));
    }
    //if (verbose == 1) fclose(logfile); return;
    for (m = 0; m < h * w * d * s * 4; m++)
    {
        if (get<5>(EdgeList[m]) > 1e10) break;
        //i, j, k: coordinates of first voxel in the edge
        i = get<0>(EdgeList[m]);
        j = get<1>(EdgeList[m]);
        k = get<2>(EdgeList[m]);
        p = get<3>(EdgeList[m]);
        t = get<4>(EdgeList[m]);
        // ci, cj, ck: coordinates of second voxel in the edge
        ci = i; cj = j; ck = k; cp = p;
        if (t == 0) ci++;
        else if (t == 1) cj++;
        else if (t == 2) ck++;
        else cp++;
        //g_i, g_ci: linearized column-major order indexes of (i,j,k,p) and (ci,cj,ck,cp)
        g_i = i + j * h + k * h * w + p * h * w * d;
        g_ci = ci + cj * h + ck * h * w + cp * h * w * d;
        
        // decode indexes of the first voxel in the group of first voxel in the edge        
        vi = ((VoxelGroupFirst[g_i] % (h * w * d)) % (h * w)) % h;
        vj = (int)((VoxelGroupFirst[g_i] % (h * w * d)) % (h * w)) / h;
        vk = (int)(VoxelGroupFirst[g_i] % (h * w * d)) / (h * w);
        vp = (int)VoxelGroupFirst[g_i] / (h * w * d);
        // decode indexes of the first voxel in the group of second voxel in the edge
        vci = ((VoxelGroupFirst[g_ci] % (h * w * d)) % (h * w)) % h;
        vcj = (int)((VoxelGroupFirst[g_ci] % (h * w * d)) % (h * w)) / h;
        vck = (int)(VoxelGroupFirst[g_ci] % (h * w * d)) / (h * w);
        vcp = (int)VoxelGroupFirst[g_ci] / (h * w * d);
        if (verbose == 1) fprintf(logfile,"proc %d: %d %d %d %d (1st %d %d %d %d) - %d %d %d %d (1st %d %d %d %d) val %.8e\n",m,i,j,k,p,vi,vj,vk,vp,ci,cj,ck,cp,vci,vcj,vck,vcp,get<5>(EdgeList[m]));        
        
        // If both (i,j,k,p) and (ci, cj, ck, cp) are already in the same group, nothing to do
        // This edge will not be used
        if (VoxelGroupFirst[g_i] == VoxelGroupFirst[g_ci]) continue;
        
        // if (ci, cj, ck, cp) voxel have not been unwrapped before, 
        // then unwrap (ci, cj, ck, cp) with respect to (i, j, k, p) and add it to its group;  
        if (VoxelGroupSize[VoxelGroupFirst[g_ci]] == 1)
        {
            // PhaseUnwrap2px(x,y) returns the unwrapped value of y with respect to x;
            // Unwrap (ci,cj,ck,cp) with respect to (i,j,k,p)
            IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd] = PhaseUnwrap2px(IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd], IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd]);
            if (VoxelGroupSize[VoxelGroupFirst[g_i]] == 1)
                {if (verbose == 1) fprintf(logfile,"ref. %d %d %d %d unwrap %d %d %d %d\n", i,j,k,p, ci,cj,ck,cp);}
            else
                if (verbose == 1) fprintf(logfile,"ref. %d %d %d %d (group of %d el.) unwrap %d %d %d %d\n",i,j,k,p, VoxelGroupSize[VoxelGroupFirst[g_i]], ci,cj,ck,cp);
            // Add voxel (ci,cj,ck,cp) to the group of (i,j,k,p), update reference tables
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_i]]] = g_ci;
            VoxelGroupLast[VoxelGroupFirst[g_i]] = g_ci;
            VoxelGroupSize[VoxelGroupFirst[g_i]]++;
            VoxelGroupFirst[g_ci] = VoxelGroupFirst[g_i];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_i];
            continue;
        }

        // if (i, j, k, p) voxel have not been unwrapped before,
        // then unwrap (i, j, k, p) with respect to (ci, cj, ck, cp) and add it to its group;
        if (VoxelGroupSize[VoxelGroupFirst[g_i]] == 1)
        {
            // Unwrap (i,j,k,p) with respect to (ci,cj,ck,cp)
            IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd] = PhaseUnwrap2px(IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd], IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd]);
            if (verbose == 1) fprintf(logfile,"unwrap %d %d %d %d ref. %d %d %d %d (group of %d el.)\n", i,j,k,p, ci,cj,ck,cp, VoxelGroupSize[VoxelGroupFirst[g_ci]]);
            // Add voxel (i,j,k,p) to the group of (ci,cj,ck,cp), update reference tables
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_ci]]] = g_i;
            VoxelGroupLast[VoxelGroupFirst[g_ci]] = g_i;
            VoxelGroupSize[VoxelGroupFirst[g_ci]]++;
            VoxelGroupFirst[g_i] = VoxelGroupFirst[g_ci];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_ci]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_ci];
            continue;
        }        
       
        // if both (i, j, k, p) and (ci, cj, ck, cp) have been unwrapped before and belong to different groups,
        // unwrap and add group of (ci, cj, ck, cp) to the group of (i, j, k, p) if the group of (i, j, k, p) is bigger
        if (VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[VoxelGroupFirst[g_ci]])
        {
            if (verbose == 1) fprintf(logfile,"add to group of %d %d %d %d (%d el., other %d el.):\n",i,j,k,p,VoxelGroupSize[VoxelGroupFirst[g_i]],VoxelGroupSize[VoxelGroupFirst[g_ci]]);
            // when we unwrap the edge linking two groups, the value of the voxel in the second group is shifted by 'val'
            // we save this 'val' value to subsequently add it to all voxels in the second group (of voxel (ci, cj, ck, cp)), 
            // thus not breaking the relative phase differences between voxels in that group
            val = PhaseUnwrap2px(IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd], IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd]) - IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd];
            // Joining the group of (ci, cj, ck, cp) with the group of (i, j, k, p)
            // Link the last voxel in the group of (i, j, k, p) to the first voxel in the group of (ci, cj, ck, cp)
            VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_i]]] = VoxelGroupFirst[g_ci];
            // Update the last voxel in the group of (i, j, k, p): now it is the last voxel in the group of (ci, cj, ck, cp)
            VoxelGroupLast[VoxelGroupFirst[g_i]] = VoxelGroupLast[VoxelGroupFirst[g_ci]];
            // Update the size of the group of voxel (i, j, k, p)  
            VoxelGroupSize[VoxelGroupFirst[g_i]] += VoxelGroupSize[VoxelGroupFirst[g_ci]];
            // Update eventually the reference g_max to the largest group
            if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_i]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_i];
            // Index of the first voxel in the group of (ci, cj, ck, cp) 
            n = VoxelGroupFirst[g_ci];
            while (n != -1)
            {
                // Get ci,cj,ck,cp coordinates of the current voxel in the group of (ci, cj, ck, cp) 
                ci = ((n % (h * w * d)) % (h * w)) % h;
                cj = (int)((n % (h * w * d)) % (h * w)) / h;
                ck = (int)(n % (h * w * d)) / (h * w);
                cp = (int)n / (h * w * d);
                if (verbose == 1) fprintf(logfile," + %d %d %d %d:\n",ci,cj,ck,cp);
                // Update the reference to the first voxel in the group
                // as both groups are merged to the group of (i, j, k, p) voxel
                VoxelGroupFirst[n] = VoxelGroupFirst[g_i];
                // Add the 'val' offset to the phase value of the current voxel in the group of (ci, cj, ck, cp) 
                IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd] += val;
                // Move to next voxel in the group of (ci, cj, ck, cp) 
                n = VoxelGroupNext[n];
            }
            continue;
        }
        // unwrap and add group of (i, j, k, p) to the group of (ci, cj, ck, cp) if the group of (ci, cj, ck, cp) is bigger
        if (verbose == 1) fprintf(logfile,"add to group of %d %d %d %d (%d el., other %d el.):\n",ci,cj,ck,cp,VoxelGroupSize[VoxelGroupFirst[g_ci]],VoxelGroupSize[VoxelGroupFirst[g_i]]);
        // when we unwrap the edge linking two groups, the value of the voxel in the second group is shifted by 'val'
        // we save this 'val' value to subsequently add it to all voxels in the second group (of voxel (i, j, k)),
        // thus not breaking the relative phase differences between voxels in that group
        val = PhaseUnwrap2px(IM_phase[ci + cj * imh + ck * imh * imw + cp * imh * imw * imd], IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd]) - IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd];
        // Joining the group of (i, j, k, p) with the group of (ci, cj, ck, cp)
        // Link the last voxel in the group of (ci, cj, ck, cp) to the first voxel in the group of (i, j, k, p)
        VoxelGroupNext[VoxelGroupLast[VoxelGroupFirst[g_ci]]] = VoxelGroupFirst[g_i];
        // Update the last voxel in the group of (ci, cj, ck, cp): now it is the last voxel in the group of (i, j, k, p)
        VoxelGroupLast[VoxelGroupFirst[g_ci]] = VoxelGroupLast[VoxelGroupFirst[g_i]];
        // Update the size of the group of voxel (ci, cj, ck, cp)  
        VoxelGroupSize[VoxelGroupFirst[g_ci]] += VoxelGroupSize[VoxelGroupFirst[g_i]];
        // Update eventually the reference g_max to the largest group
        if ((g_max == -1)||(VoxelGroupSize[VoxelGroupFirst[g_ci]] > VoxelGroupSize[g_max])) g_max = VoxelGroupFirst[g_ci];
        // Index of the first voxel in the group of (i, j, k, p) 
        n = VoxelGroupFirst[g_i];
        while (n != -1)
        {
            // Get i,j,k,p coordinates of the current voxel in the group of (i, j, k, p)
            i = ((n % (h * w * d)) % (h * w)) % h;
            j = (int)((n % (h * w * d)) % (h * w)) / h;
            k = (int)(n % (h * w * d)) / (h * w);
            p = (int)n / (h * w * d);
            if (verbose == 1) fprintf(logfile," + %d %d %d %d:\n",i,j,k,p);
            // Update the reference to the first voxel in the group
            // as both groups are merged to the group of (ci, cj, ck, cp) voxel
            VoxelGroupFirst[n] = VoxelGroupFirst[g_ci];
            // Add the 'val' offset to the phase value of the current voxel in the group of (i, j, k, p) 
            IM_phase[i + j * imh + k * imh * imw + p * imh * imw * imd] += val;
            // Move to next voxel in the group of (i, j, k, p) 
            n = VoxelGroupNext[n];
        }
    }
    // If the flag is set to 1, mark in 'mask' as 1 all voxels in the biggest group of unwrapped voxels
    if ((mask_largest_unwrapped_group == 1)&&(IM_mask != nullptr))
    {
        // Reinitialize mask with 0
        fill(mask, mask + imh * imw * imd * ims, 0);
        // Go to the first voxel in the biggest group
        n = VoxelGroupFirst[g_max];
        while (n != -1)
        {
            // Get i,j,k,p coordinates of the current voxel in the biggest group
            i = ((n % (h * w * d)) % (h * w)) % h;
            j = (int)((n % (h * w * d)) % (h * w)) / h;
            k = (int)(n % (h * w * d)) / (h * w);
            p = (int)n / (h * w * d);
            // Set 'mask' to 1
            IM_mask[i + j * imh + k * imh * imw + p * imh * imw * imd] = 1;
            // Move to next voxel in the biggest group
            n = VoxelGroupNext[n];
        }   
    }
    
    if (verbose == 1) fclose(logfile);
}

#ifdef __cplusplus
};
#endif

/* this function exists so that mex may be used to compile the library
   it is not otherwise needed */

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
}

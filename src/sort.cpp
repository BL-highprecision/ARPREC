/*
 * src/mpreal4.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Additional routines, mostly those for which speed is less important.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

void mp_real::mpsort(int n, mp_real *a, int *ip)
{
  /**
   * This routine sorts the entries of the N-long
   * MP vector A into ascending order using the quicksort agorithm.
   * The permutation vector that would sort the vector is returned in IP.
   */
  int left, right, lowest, highest;
  int low_stack[50], high_stack[50];
  int stack_counter=0;
  int temp, pivot_spot;
  mp_real pivot_value;
  int i;
  if(n<=1)
    return;

  /* setup phase */
  low_stack[stack_counter] = 0;
  high_stack[stack_counter] = n-1;

  for(i=0;i<n;i++) 
    ip[i] = i;
  /* main loop */
  while(stack_counter >= 0) {
    lowest = left = low_stack[stack_counter];
    highest = right = high_stack[stack_counter];
    
    if((right - left) <= 0) {
      stack_counter--;
      continue;
    }
    if((right-left) == 1) {
      if(a[ip[left]] > a[ip[right]]) {
	//swap 
	temp = ip[left];
	ip[right] = ip[left];
	ip[left] = temp;
      }
      stack_counter--;
      continue;
    }

    /*partition phase*/
    /* pick pivot from the middle- it is likley that the a vector 
       is almost in order */
    pivot_spot = (left + right + 1) /2;
    //Swap left with pivot.
    temp = ip[left];
    ip[left] = ip[pivot_spot];
    ip[pivot_spot] = temp; 
    pivot_spot = left;
    
    pivot_value = a[ip[left++]];
    while(left < right) {
      while(a[ip[left]] < pivot_value && left<highest)
	left++;
      while(a[ip[right]] >= pivot_value && right >lowest)
	right--;
      if(left<right) {
	temp = ip[left];
	ip[left] = ip[right];
	ip[right] = temp;
      }
    }
    //swap pivot and ip[right]
    temp = ip[pivot_spot];
    ip[pivot_spot] = ip[right];
    ip[right] = temp;
    /*additional sort prep phase */

    /* write over current spot on stack */
    if(lowest < right - 1) {
      low_stack[stack_counter] = lowest;
      high_stack[stack_counter] = right-1;
      stack_counter++;
    }
    if(right+1 < highest) {
      low_stack[stack_counter] = right+1;
      high_stack[stack_counter] = highest;
    } else {
      stack_counter--;
    }
  }
  return;  
}


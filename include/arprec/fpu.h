/*
 * include/arprec/fpu.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001-02
 *
 * Contains functions necessary to set and restore the round-to-double
 * flat in the FPU.  Some processors, notably the Intel x86 processors
 * (and their clones) have floating point register wider than double, 
 * that can cause problems with some of the algorithms in the arprec 
 * package.
 *
 * Currently only x86 processors (under windows or linux) are supported, 
 * but similar fix may be needed for other processors with wider-than-double 
 * registers, e.g., Motorola 680x0 series and Intel itanium series.
 */
#ifndef ARPREC_FPU_H
#define ARPREC_FPU_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Set the round-to-double flag, and save the old control word in old_cw.
 * If old_cw is NULL, the old control word is not saved.
 */
void fpu_fix_start(unsigned int *old_cw);

/*
 * Restore the control word.
 */
void fpu_fix_end(unsigned int *old_cw);

#ifdef __cplusplus
}
#endif

#endif  
/* ARPREC_FPU_H */


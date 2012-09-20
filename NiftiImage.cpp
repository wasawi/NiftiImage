#include "NiftiImage.h"   /* typedefs, prototypes, macros, etc. */

static const std::string lib_hist = "\
----------------------------------------------------------------------\n\
History NiftiImage:\n\
0.1 18 Sep 2012 [tcw]: Started NiftiImage from niftilib version 1.43  \n\
----------------------------------------------------------------------\n";
static const std::string lib_version = "NiftiImage version 0.1 (18 Sep, 2012)";

const DTMap NiftiImage::DTypes
{
  {NIFTI_TYPE_UINT8,    {1, 0, "NIFTI_TYPE_UINT8"} },
  {NIFTI_TYPE_INT16,    {2, 2, "NIFTI_TYPE_INT16"} },
  {NIFTI_TYPE_INT32,    {4, 4, "NIFTI_TYPE_INT32"} },
  {NIFTI_TYPE_FLOAT32,   {4, 4, "NIFTI_TYPE_FLOAT32"} },
  {NIFTI_TYPE_COMPLEX64,   {8, 4, "NIFTI_TYPE_COMPLEX64"} },
  {NIFTI_TYPE_FLOAT64,   {8, 8, "NIFTI_TYPE_FLOAT64"} },
  {NIFTI_TYPE_RGB24,  {3, 0, "NIFTI_TYPE_RGB24"} },
  {NIFTI_TYPE_INT8,  {1, 0, "NIFTI_TYPE_INT8"} },
  {NIFTI_TYPE_UINT16,  {2, 2, "NIFTI_TYPE_UINT16"} },
  {NIFTI_TYPE_UINT32,  {4, 4, "NIFTI_TYPE_UINT32"} },
  {NIFTI_TYPE_INT64, {8, 8, "NIFTI_TYPE_INT64"} },
  {NIFTI_TYPE_UINT64, {8, 8, "NIFTI_TYPE_UINT64"} },
  {NIFTI_TYPE_FLOAT128, {16, 16, "NIFTI_TYPE_FLOAT128"} },
  {NIFTI_TYPE_COMPLEX128, {16,  8, "NIFTI_TYPE_COMPLEX128"} },
  {NIFTI_TYPE_COMPLEX256, {32, 16, "NIFTI_TYPE_COMPLEX256"} },
  {NIFTI_TYPE_RGBA32, {4,   0, "NIFTI_TYPE_RGBA32"} }
};

/*---------------------------------------------------------------------*/
/*! Given a NIFTI_TYPE value, such as NIFTI_TYPE_INT16, return the
 *  corresponding label as a string.  The dtype code is the
 *  macro value defined in nifti1.h.
 *//*-------------------------------------------------------------------*/
const std::string &NiftiImage::DTypeToString(const int dtype)
{
	static const std::string notfound{"DT_NOTFOUND"};
	DTMap::const_iterator it = DTypes.find(dtype);
	if (it == DTypes.end())
		return notfound;
	else
		return it->second.name;
}


/*---------------------------------------------------------------------*/
/*! Determine whether dtype is a valid NIFTI_TYPE.
 *
 *  DT_UNKNOWN is considered invalid
 *
 *  The only difference 'for_nifti' makes is that DT_BINARY
 *  should be invalid for a NIfTI dataset.
 *//*-------------------------------------------------------------------*/
const bool NiftiImage::DTypeIsValid(const int dtype)
{
    DTMap::const_iterator it = DTypes.find(dtype);
	if (it == DTypes.end())
		return false;
	else
		return true;
}

/*---------------------------------------------------------------------*/
/*! Display the nifti_type_list table.
 *//*-------------------------------------------------------------------*/
void NiftiImage::printDTypeList()
{
	std::cout << "NiftiImage Valid Datatypes" << std::endl;
	std::cout << "Name                Type    Size (bytes) Swapsize (bytes)" << std::endl;
	
	DTMap::const_iterator it;
	for (it = DTypes.begin(); it != DTypes.end(); it++)
	{
		std::cout << it->second.name << " "
				  << it->first << " "
				  << it->second.size << " "
				  << it->second.swapsize << std::endl;
	}
}

/*---------------------------------------------------------------------------*/
/* prototypes for internal functions - not part of exported library          */

/* extension routines */
/*static int  nifti_read_extensions( nifti_image *nim, znzFile fp, int remain );
static int  nifti_read_next_extension( nifti1_extension * nex, nifti_image *nim,                                       int remain, znzFile fp );
static int  nifti_check_extension(nifti_image *nim, int size,int code, int rem);
static int  nifti_add_exten_to_list(nifti1_extension *  new_ext,
                                    nifti1_extension ** list, long new_length);
static int  nifti_fill_extension(nifti1_extension * ext, const char * data,
                                 int len, int ecode);*/

/* for nifti_read_collapsed_image: */
/*static int  rci_read_data(nifti_image *nim, int *pivots, int *prods, int nprods,
						  const int dims[], char *data, znzFile fp, size_t base_offset);
static int  rci_alloc_mem(void ** data, int prods[8], int nprods, int nbyper );
static int  make_pivot_list(nifti_image * nim, const int dims[], int pivots[],
                            int prods[], int * nprods );*/
/*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI units type.
 
 \param  uu NIfTI-1 unit code
 
 \return pointer to static string for the given unit type
 
 \warning Do not free() or modify this string!
 It points to static storage.
 
 \sa     NIFTI1_UNITS group in nifti1.h
 *//*-------------------------------------------------------------------------*/
const char *nifti_units_string( int uu )
{
	switch( uu ){
		case NIFTI_UNITS_METER:  return "m" ;
		case NIFTI_UNITS_MM:     return "mm" ;
		case NIFTI_UNITS_MICRON: return "um" ;
		case NIFTI_UNITS_SEC:    return "s" ;
		case NIFTI_UNITS_MSEC:   return "ms" ;
		case NIFTI_UNITS_USEC:   return "us" ;
		case NIFTI_UNITS_HZ:     return "Hz" ;
		case NIFTI_UNITS_PPM:    return "ppm" ;
		case NIFTI_UNITS_RADS:   return "rad/s" ;
	}
	return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI transform type.
 
 \param  xx NIfTI-1 xform code
 
 \return pointer to static string describing xform code
 
 \warning Do not free() or modify this string!
 It points to static storage.
 
 \sa     NIFTI1_XFORM_CODES group in nifti1.h
 *//*-------------------------------------------------------------------------*/
const char *nifti_xform_string( int xx )
{
	switch( xx ){
		case NIFTI_XFORM_SCANNER_ANAT:  return "Scanner Anat" ;
		case NIFTI_XFORM_ALIGNED_ANAT:  return "Aligned Anat" ;
		case NIFTI_XFORM_TALAIRACH:     return "Talairach" ;
		case NIFTI_XFORM_MNI_152:       return "MNI_152" ;
	}
	return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI intent type.
 
 \param  ii NIfTI-1 intent code
 
 \return pointer to static string describing code
 
 \warning Do not free() or modify this string!
 It points to static storage.
 
 \sa     NIFTI1_INTENT_CODES group in nifti1.h
 *//*-------------------------------------------------------------------------*/
const char *nifti_intent_string( int ii )
{
	switch( ii ){
		case NIFTI_INTENT_CORREL:     return "Correlation statistic" ;
		case NIFTI_INTENT_TTEST:      return "T-statistic" ;
		case NIFTI_INTENT_FTEST:      return "F-statistic" ;
		case NIFTI_INTENT_ZSCORE:     return "Z-score"     ;
		case NIFTI_INTENT_CHISQ:      return "Chi-squared distribution" ;
		case NIFTI_INTENT_BETA:       return "Beta distribution" ;
		case NIFTI_INTENT_BINOM:      return "Binomial distribution" ;
		case NIFTI_INTENT_GAMMA:      return "Gamma distribution" ;
		case NIFTI_INTENT_POISSON:    return "Poisson distribution" ;
		case NIFTI_INTENT_NORMAL:     return "Normal distribution" ;
		case NIFTI_INTENT_FTEST_NONC: return "F-statistic noncentral" ;
		case NIFTI_INTENT_CHISQ_NONC: return "Chi-squared noncentral" ;
		case NIFTI_INTENT_LOGISTIC:   return "Logistic distribution" ;
		case NIFTI_INTENT_LAPLACE:    return "Laplace distribution" ;
		case NIFTI_INTENT_UNIFORM:    return "Uniform distribition" ;
		case NIFTI_INTENT_TTEST_NONC: return "T-statistic noncentral" ;
		case NIFTI_INTENT_WEIBULL:    return "Weibull distribution" ;
		case NIFTI_INTENT_CHI:        return "Chi distribution" ;
		case NIFTI_INTENT_INVGAUSS:   return "Inverse Gaussian distribution" ;
		case NIFTI_INTENT_EXTVAL:     return "Extreme Value distribution" ;
		case NIFTI_INTENT_PVAL:       return "P-value" ;
			
		case NIFTI_INTENT_LOGPVAL:    return "Log P-value" ;
		case NIFTI_INTENT_LOG10PVAL:  return "Log10 P-value" ;
			
		case NIFTI_INTENT_ESTIMATE:   return "Estimate" ;
		case NIFTI_INTENT_LABEL:      return "Label index" ;
		case NIFTI_INTENT_NEURONAME:  return "NeuroNames index" ;
		case NIFTI_INTENT_GENMATRIX:  return "General matrix" ;
		case NIFTI_INTENT_SYMMATRIX:  return "Symmetric matrix" ;
		case NIFTI_INTENT_DISPVECT:   return "Displacement vector" ;
		case NIFTI_INTENT_VECTOR:     return "Vector" ;
		case NIFTI_INTENT_POINTSET:   return "Pointset" ;
		case NIFTI_INTENT_TRIANGLE:   return "Triangle" ;
		case NIFTI_INTENT_QUATERNION: return "Quaternion" ;
			
		case NIFTI_INTENT_DIMLESS:    return "Dimensionless number" ;
	}
	return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI slice_code.
 
 \param  ss NIfTI-1 slice order code
 
 \return pointer to static string describing code
 
 \warning Do not free() or modify this string!
 It points to static storage.
 
 \sa     NIFTI1_SLICE_ORDER group in nifti1.h
 *//*-------------------------------------------------------------------------*/
const char *nifti_slice_string( int ss )
{
	switch( ss ){
		case NIFTI_SLICE_SEQ_INC:  return "sequential_increasing"    ;
		case NIFTI_SLICE_SEQ_DEC:  return "sequential_decreasing"    ;
		case NIFTI_SLICE_ALT_INC:  return "alternating_increasing"   ;
		case NIFTI_SLICE_ALT_DEC:  return "alternating_decreasing"   ;
		case NIFTI_SLICE_ALT_INC2: return "alternating_increasing_2" ;
		case NIFTI_SLICE_ALT_DEC2: return "alternating_decreasing_2" ;
	}
	return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI orientation.
 
 \param ii orientation code
 
 \return pointer to static string holding the orientation information
 
 \warning Do not free() or modify the return string!
 It points to static storage.
 
 \sa  NIFTI_L2R in nifti1_io.h
 *//*-------------------------------------------------------------------------*/
const char *nifti_orientation_string( int ii )
{
	switch( ii ){
		case NIFTI_L2R: return "Left-to-Right" ;
		case NIFTI_R2L: return "Right-to-Left" ;
		case NIFTI_P2A: return "Posterior-to-Anterior" ;
		case NIFTI_A2P: return "Anterior-to-Posterior" ;
		case NIFTI_I2S: return "Inferior-to-Superior" ;
		case NIFTI_S2I: return "Superior-to-Inferior" ;
	}
	return "Unknown" ;
}
/*---------------------------------------------------------------------------*/
/*! compute the (closest) orientation from a 4x4 ijk->xyz tranformation matrix
 
 <pre>
 Input:  4x4 matrix that transforms (i,j,k) indexes to (x,y,z) coordinates,
 where +x=Right, +y=Anterior, +z=Superior.
 (Only the upper-left 3x3 corner of R is used herein.)
 Output: 3 orientation codes that correspond to the closest "standard"
 anatomical orientation of the (i,j,k) axes.
 Method: Find which permutation of (x,y,z) has the smallest angle to the
 (i,j,k) axes directions, which are the columns of the R matrix.
 Errors: The codes returned will be zero.
 
 For example, an axial volume might get return values of
 *icod = NIFTI_R2L   (i axis is mostly Right to Left)
 *jcod = NIFTI_P2A   (j axis is mostly Posterior to Anterior)
 *kcod = NIFTI_I2S   (k axis is mostly Inferior to Superior)
 </pre>
 
 \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h
 
 \see nifti_quatern_to_mat44, nifti_mat44_to_quatern,
 nifti_make_orthog_mat44
 *//*-------------------------------------------------------------------------*/
void nifti_matrix4d_to_orientation( Matrix4d R , int *icod, int *jcod, int *kcod )
{
	float xi,xj,xk , yi,yj,yk , zi,zj,zk , val,detQ,detP ;
	Matrix3d P , Q , M ;
	int i,j,k=0,p,q,r , ibest,jbest,kbest,pbest,qbest,rbest ;
	float vbest ;
	
	if( icod == NULL || jcod == NULL || kcod == NULL ) return ; /* bad */
	
	*icod = *jcod = *kcod = 0 ; /* error returns, if sh*t happens */
	
	/* load column vectors for each (i,j,k) direction from matrix */
	
	/*-- i axis --*/ /*-- j axis --*/ /*-- k axis --*/
	
	xi = R(0, 0) ; xj = R(0, 1) ; xk = R(0, 2) ;
	yi = R(1, 0) ; yj = R(1, 1) ; yk = R(1, 2) ;
	zi = R(2, 0) ; zj = R(2, 1) ; zk = R(2, 2) ;
	
	/* normalize column vectors to get unit vectors along each ijk-axis */
	
	/* normalize i axis */
	
	val = sqrt( xi*xi + yi*yi + zi*zi ) ;
	if( val == 0.0 ) return ;                 /* stupid input */
	xi /= val ; yi /= val ; zi /= val ;
	
	/* normalize j axis */
	
	val = sqrt( xj*xj + yj*yj + zj*zj ) ;
	if( val == 0.0 ) return ;                 /* stupid input */
	xj /= val ; yj /= val ; zj /= val ;
	
	/* orthogonalize j axis to i axis, if needed */
	
	val = xi*xj + yi*yj + zi*zj ;    /* dot product between i and j */
	if( fabs(val) > 1.e-4 ){
		xj -= val*xi ; yj -= val*yi ; zj -= val*zi ;
		val = sqrt( xj*xj + yj*yj + zj*zj ) ;  /* must renormalize */
		if( val == 0.0 ) return ;              /* j was parallel to i? */
		xj /= val ; yj /= val ; zj /= val ;
	}
	
	/* normalize k axis; if it is zero, make it the cross product i x j */
	
	val = sqrt( xk*xk + yk*yk + zk*zk ) ;
	if( val == 0.0 ){ xk = yi*zj-zi*yj; yk = zi*xj-zj*xi ; zk=xi*yj-yi*xj ; }
	else            { xk /= val ; yk /= val ; zk /= val ; }
	
	/* orthogonalize k to i */
	
	val = xi*xk + yi*yk + zi*zk ;    /* dot product between i and k */
	if( fabs(val) > 1.e-4 ){
		xk -= val*xi ; yk -= val*yi ; zk -= val*zi ;
		val = sqrt( xk*xk + yk*yk + zk*zk ) ;
		if( val == 0.0 ) return ;      /* bad */
		xk /= val ; yk /= val ; zk /= val ;
	}
	
	/* orthogonalize k to j */
	
	val = xj*xk + yj*yk + zj*zk ;    /* dot product between j and k */
	if( fabs(val) > 1.e-4 ){
		xk -= val*xj ; yk -= val*yj ; zk -= val*zj ;
		val = sqrt( xk*xk + yk*yk + zk*zk ) ;
		if( val == 0.0 ) return ;      /* bad */
		xk /= val ; yk /= val ; zk /= val ;
	}
	
	Q(0, 0) = xi ; Q(0, 1) = xj ; Q(0, 2) = xk ;
	Q(1, 0) = yi ; Q(1, 1) = yj ; Q(1, 2) = yk ;
	Q(2, 0) = zi ; Q(2, 1) = zj ; Q(2, 2) = zk ;
	
	/* at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes */
	
	detQ = Q.determinant();
	if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */
	
	/* Build and test all possible +1/-1 coordinate permutation matrices P;
	 then find the P such that the rotation matrix M=PQ is closest to the
	 identity, in the sense of M having the smallest total rotation angle. */
	
	/* Despite the formidable looking 6 nested loops, there are
	 only 3*3*3*2*2*2 = 216 passes, which will run very quickly. */
	
	vbest = -666.0 ; ibest=pbest=qbest=rbest=1 ; jbest=2 ; kbest=3 ;
	for( i=1 ; i <= 3 ; i++ ){     /* i = column number to use for row #1 */
		for( j=1 ; j <= 3 ; j++ ){    /* j = column number to use for row #2 */
			if( i == j ) continue ;
			for( k=1 ; k <= 3 ; k++ ){  /* k = column number to use for row #3 */
				if( i == k || j == k ) continue ;
				P(0, 0) = P(0, 1) = P(0, 2) =
				P(1, 0) = P(1, 1) = P(1, 2) =
				P(2, 0) = P(2, 1) = P(2, 2) = 0.0 ;
				for( p=-1 ; p <= 1 ; p+=2 ){    /* p,q,r are -1 or +1      */
					for( q=-1 ; q <= 1 ; q+=2 ){   /* and go into rows #1,2,3 */
						for( r=-1 ; r <= 1 ; r+=2 ){
							P(0, i-1) = p ; P(1, j-1) = q ; P(2, k-1) = r ;
							detP = P.determinant() ;           /* sign of permutation */
							if( detP * detQ <= 0.0 ) continue ;  /* doesn't match sign of Q */
							M = P * Q ;
							
							/* angle of M rotation = 2.0*acos(0.5*sqrt(1.0+trace(M)))       */
							/* we want largest trace(M) == smallest angle == M nearest to I */
							
							val = M(0, 0) + M(1, 1) + M(2, 2) ; /* trace */
							if( val > vbest ){
								vbest = val ;
								ibest = i ; jbest = j ; kbest = k ;
								pbest = p ; qbest = q ; rbest = r ;
							}
						}}}}}}
	
	/* At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.
	 
	 The matrix P that corresponds is the best permutation approximation
	 to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
	 to the (i,j,k) axes.
	 
	 For example, the first row of P (which contains pbest in column ibest)
	 determines the way the i axis points relative to the anatomical
	 (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
	 which is direction P2A (if pbest > 0) or A2P (if pbest < 0).
	 
	 So, using ibest and pbest, we can assign the output code for
	 the i axis.  Mutatis mutandis for the j and k axes, of course. */
	
	switch( ibest*pbest ){
		case  1: i = NIFTI_L2R ; break ;
		case -1: i = NIFTI_R2L ; break ;
		case  2: i = NIFTI_P2A ; break ;
		case -2: i = NIFTI_A2P ; break ;
		case  3: i = NIFTI_I2S ; break ;
		case -3: i = NIFTI_S2I ; break ;
	}
	
	switch( jbest*qbest ){
		case  1: j = NIFTI_L2R ; break ;
		case -1: j = NIFTI_R2L ; break ;
		case  2: j = NIFTI_P2A ; break ;
		case -2: j = NIFTI_A2P ; break ;
		case  3: j = NIFTI_I2S ; break ;
		case -3: j = NIFTI_S2I ; break ;
	}
	
	switch( kbest*rbest ){
		case  1: k = NIFTI_L2R ; break ;
		case -1: k = NIFTI_R2L ; break ;
		case  2: k = NIFTI_P2A ; break ;
		case -2: k = NIFTI_A2P ; break ;
		case  3: k = NIFTI_I2S ; break ;
		case -3: k = NIFTI_S2I ; break ;
	}
	
	*icod = i ; *jcod = j ; *kcod = k ; return ;
}

/*----------------------------------------------------------------------*/
/*! generic: swap siz bytes at a time from the given list of n sets
 *//*--------------------------------------------------------------------*/
void nifti_swap_bytes( size_t n , int siz , void *ar )
{
	register size_t ii ;
	unsigned char * cp0 = (unsigned char *)ar, * cp1, * cp2 ;
	register unsigned char tval ;
	
	for( ii=0 ; ii < n ; ii++ ){
		cp1 = cp0;  cp2 = cp0+(siz-1);
		while ( cp2 > cp1 )
		{
			tval = *cp1 ; *cp1 = *cp2 ; *cp2 = tval ;
			cp1++; cp2--;
		}
		cp0 += siz;
	}
	return ;
}

/*-------------------------------------------------------------------------*/
/*! Byte swap NIFTI-1 file header in various places and ways.
 
 If is_nifti, swap all (even UNUSED) fields of NIfTI header.
 Else, swap as a nifti_analyze75 struct.
 *//*---------------------------------------------------------------------- */
void swap_nifti_header( struct nifti_1_header *h , int is_nifti )
{
	
	/* if ANALYZE, swap as such and return */
	if( ! is_nifti ) {
		nifti_swap_as_analyze((nifti_analyze75 *)h);
		return;
	}
	
	/* otherwise, swap all NIFTI fields */
	
	nifti_swap_bytes(1, 4, &h->sizeof_hdr);
	nifti_swap_bytes(1, 4, &h->extents);
	nifti_swap_bytes(1, 2, &h->session_error);
	
	nifti_swap_bytes(8, 2, h->dim);
	nifti_swap_bytes(1, 4, &h->intent_p1);
	nifti_swap_bytes(1, 4, &h->intent_p2);
	nifti_swap_bytes(1, 4, &h->intent_p3);
	
	nifti_swap_bytes(1, 2, &h->intent_code);
	nifti_swap_bytes(1, 2, &h->datatype);
	nifti_swap_bytes(1, 2, &h->bitpix);
	nifti_swap_bytes(1, 2, &h->slice_start);
	
	nifti_swap_bytes(8, 4, h->pixdim);
	
	nifti_swap_bytes(1, 4, &h->vox_offset);
	nifti_swap_bytes(1, 4, &h->scl_slope);
	nifti_swap_bytes(1, 4, &h->scl_inter);
	nifti_swap_bytes(1, 2, &h->slice_end);
	
	nifti_swap_bytes(1, 4, &h->cal_max);
	nifti_swap_bytes(1, 4, &h->cal_min);
	nifti_swap_bytes(1, 4, &h->slice_duration);
	nifti_swap_bytes(1, 4, &h->toffset);
	nifti_swap_bytes(1, 4, &h->glmax);
	nifti_swap_bytes(1, 4, &h->glmin);
	
	nifti_swap_bytes(1, 2, &h->qform_code);
	nifti_swap_bytes(1, 2, &h->sform_code);
	
	nifti_swap_bytes(1, 4, &h->quatern_b);
	nifti_swap_bytes(1, 4, &h->quatern_c);
	nifti_swap_bytes(1, 4, &h->quatern_d);
	nifti_swap_bytes(1, 4, &h->qoffset_x);
	nifti_swap_bytes(1, 4, &h->qoffset_y);
	nifti_swap_bytes(1, 4, &h->qoffset_z);
	
	nifti_swap_bytes(4, 4, h->srow_x);
	nifti_swap_bytes(4, 4, h->srow_y);
	nifti_swap_bytes(4, 4, h->srow_z);
	
	return ;
}

/*-------------------------------------------------------------------------*/
/*! Byte swap as an ANALYZE 7.5 header
 *
 *  return non-zero on failure
 *//*---------------------------------------------------------------------- */
int nifti_swap_as_analyze( nifti_analyze75 * h )
{
	if( !h ) return 1;
	
	nifti_swap_bytes(1, 4, &h->sizeof_hdr);
	nifti_swap_bytes(1, 4, &h->extents);
	nifti_swap_bytes(1, 2, &h->session_error);
	
	nifti_swap_bytes(8, 2, h->dim);
	nifti_swap_bytes(1, 2, &h->unused8);
	nifti_swap_bytes(1, 2, &h->unused9);
	nifti_swap_bytes(1, 2, &h->unused10);
	nifti_swap_bytes(1, 2, &h->unused11);
	nifti_swap_bytes(1, 2, &h->unused12);
	nifti_swap_bytes(1, 2, &h->unused13);
	nifti_swap_bytes(1, 2, &h->unused14);
	
	nifti_swap_bytes(1, 2, &h->datatype);
	nifti_swap_bytes(1, 2, &h->bitpix);
	nifti_swap_bytes(1, 2, &h->dim_un0);
	
	nifti_swap_bytes(8, 4, h->pixdim);
	
	nifti_swap_bytes(1, 4, &h->vox_offset);
	nifti_swap_bytes(1, 4, &h->funused1);
	nifti_swap_bytes(1, 4, &h->funused2);
	nifti_swap_bytes(1, 4, &h->funused3);
	
	nifti_swap_bytes(1, 4, &h->cal_max);
	nifti_swap_bytes(1, 4, &h->cal_min);
	nifti_swap_bytes(1, 4, &h->compressed);
	nifti_swap_bytes(1, 4, &h->verified);
	
	nifti_swap_bytes(1, 4, &h->glmax);
	nifti_swap_bytes(1, 4, &h->glmin);
	
	nifti_swap_bytes(1, 4, &h->views);
	nifti_swap_bytes(1, 4, &h->vols_added);
	nifti_swap_bytes(1, 4, &h->start_field);
	nifti_swap_bytes(1, 4, &h->field_skip);
	
	nifti_swap_bytes(1, 4, &h->omax);
	nifti_swap_bytes(1, 4, &h->omin);
	nifti_swap_bytes(1, 4, &h->smax);
	nifti_swap_bytes(1, 4, &h->smin);
	
	return 0;
}

/*--------------------------------------------------------------------------*/
/*! check whether the given type is on the "approved" list
 
 The code is valid if it is non-negative, and does not exceed
 NIFTI_MAX_FTYPE.
 
 \return 1 if nifti_type is valid, 0 otherwise
 \sa NIFTI_FTYPE_* codes in nifti1_io.h
 *//*------------------------------------------------------------------------*/
int is_valid_nifti_type( int nifti_type )
{
	if( nifti_type >= NIFTI_FTYPE_ANALYZE &&   /* smallest type, 0 */
       nifti_type <= NIFTI_MAX_FTYPE )
		return 1;
	return 0;
}

/*----------------------------------------------------------------------*/
/*! display the contents of the nifti_1_header (send to stdout)
 
 \param info if non-NULL, print this character string
 \param hp   pointer to nifti_1_header
 *//*--------------------------------------------------------------------*/
int disp_nifti_1_header( const char * info, const nifti_1_header * hp )
{
	int c;
	
	fputs( "-------------------------------------------------------\n", stdout );
	if ( info )  fputs( info, stdout );
	if ( !hp  ){ fputs(" ** no nifti_1_header to display!\n",stdout); return 1; }
	
	fprintf(stdout," nifti_1_header :\n"
			"    sizeof_hdr     = %d\n"
			"    data_type[10]  = %s\n"
			"    db_name[18]    = %s\n",
			hp->sizeof_hdr, hp->data_type, hp->db_name);
	fprintf(stdout, "    extents        = %d\n"
			"    session_error  = %d\n"
			"    regular        = 0x%x\n"
			"    dim_info       = 0x%x\n",
			hp->extents, hp->session_error, hp->regular, hp->dim_info );
	fprintf(stdout, "    dim[8]         =");
	for ( c = 0; c < 8; c++ ) fprintf(stdout," %d", hp->dim[c]);
	fprintf(stdout, "\n"
			"    intent_p1      = %f\n"
			"    intent_p2      = %f\n"
			"    intent_p3      = %f\n"
			"    intent_code    = %d\n"
			"    datatype       = %d\n"
			"    bitpix         = %d\n"
			"    slice_start    = %d\n"
			"    pixdim[8]      =",
			hp->intent_p1, hp->intent_p2, hp->intent_p3, hp->intent_code,
			hp->datatype, hp->bitpix, hp->slice_start);
	/* break pixdim over 2 lines */
	for ( c = 0; c < 4; c++ ) fprintf(stdout," %f", hp->pixdim[c]);
	fprintf(stdout, "\n                    ");
	for ( c = 4; c < 8; c++ ) fprintf(stdout," %f", hp->pixdim[c]);
	fprintf(stdout, "\n"
			"    vox_offset     = %f\n"
			"    scl_slope      = %f\n"
			"    scl_inter      = %f\n"
			"    slice_end      = %d\n"
			"    slice_code     = %d\n"
			"    xyzt_units     = 0x%x\n"
			"    cal_max        = %f\n"
			"    cal_min        = %f\n"
			"    slice_duration = %f\n"
			"    toffset        = %f\n"
			"    glmax          = %d\n"
			"    glmin          = %d\n",
			hp->vox_offset, hp->scl_slope, hp->scl_inter, hp->slice_end,
			hp->slice_code, hp->xyzt_units, hp->cal_max, hp->cal_min,
			hp->slice_duration, hp->toffset, hp->glmax, hp->glmin);
	fprintf(stdout,
			"    descrip        = '%.80s'\n"
			"    aux_file       = '%.24s'\n"
			"    qform_code     = %d\n"
			"    sform_code     = %d\n"
			"    quatern_b      = %f\n"
			"    quatern_c      = %f\n"
			"    quatern_d      = %f\n"
			"    qoffset_x      = %f\n"
			"    qoffset_y      = %f\n"
			"    qoffset_z      = %f\n"
			"    srow_x[4]      = %f, %f, %f, %f\n"
			"    srow_y[4]      = %f, %f, %f, %f\n"
			"    srow_z[4]      = %f, %f, %f, %f\n"
			"    intent_name    = '%-.16s'\n"
			"    magic          = '%-.4s'\n",
			hp->descrip, hp->aux_file, hp->qform_code, hp->sform_code,
			hp->quatern_b, hp->quatern_c, hp->quatern_d,
			hp->qoffset_x, hp->qoffset_y, hp->qoffset_z,
			hp->srow_x[0], hp->srow_x[1], hp->srow_x[2], hp->srow_x[3],
			hp->srow_y[0], hp->srow_y[1], hp->srow_y[2], hp->srow_y[3],
			hp->srow_z[0], hp->srow_z[1], hp->srow_z[2], hp->srow_z[3],
			hp->intent_name, hp->magic);
	fputs( "-------------------------------------------------------\n", stdout );
	fflush(stdout);
	
	return 0;
}
/*----------------------------------------------------------------------
 * Read the extensions into the nifti_image struct   08 Dec 2004 [rickr]
 *
 * This function is called just after the header struct is read in, and
 * it is assumed the file pointer has not moved.  The value in remain
 * is assumed to be accurate, reflecting the bytes of space for potential
 * extensions.
 *
 * return the number of extensions read in, or < 0 on error
 *----------------------------------------------------------------------*/
//static int nifti_read_extensions( nifti_image *nim, znzFile fp, int remain )
//{
//	nifti1_extender    extdr;      /* defines extension existence  */
//	nifti1_extension   extn;       /* single extension to process  */
//	nifti1_extension * Elist;      /* list of processed extensions */
//	long               posn, count;
//	
//	if( !nim || znz_isnull(fp) ) {
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"** nifti_read_extensions: bad inputs (%p,%p)\n",
//					(void *)nim, (void *)fp);
//		return -1;
//	}
//	
//	posn = znztell(fp);
//	
//	if( (posn != sizeof(nifti_1_header)) &&
//       (nim->nifti_type != NIFTI_FTYPE_ASCII) )
//		fprintf(stderr,"** WARNING: posn not header size (%ld, %d)\n",
//				posn, (int)sizeof(nifti_1_header));
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d nre: posn = %ld, offset = %d, type = %d, remain = %d\n",
//				posn, nim->iname_offset, nim->nifti_type, remain);
//	
//	if( remain < 16 ){
//		if( g_opts.debug > 2 ){
//			if( g_opts.skip_blank_ext )
//				fprintf(stderr,"-d no extender in '%s' is okay, as "
//						"skip_blank_ext is set\n",nim->fname);
//			else
//				fprintf(stderr,"-d remain=%d, no space for extensions\n",remain);
//		}
//		return 0;
//	}
//	
//	count = (int)znzread( extdr.extension, 1, 4, fp ); /* get extender */
//	
//	if( count < 4 ){
//		if( g_opts.debug > 1 )
//			fprintf(stderr,"-d file '%s' is too short for an extender\n",
//					nim->fname);
//		return 0;
//	}
//	
//	if( extdr.extension[0] != 1 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d extender[0] (%d) shows no extensions for '%s'\n",
//					extdr.extension[0], nim->fname);
//		return 0;
//	}
//	
//	remain -= 4;
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d found valid 4-byte extender, remain = %d\n", remain);
//	
//	/* so we expect extensions, but have no idea of how many there may be */
//	
//	count = 0;
//	Elist = NULL;
//	while (nifti_read_next_extension(&extn, nim, remain, fp) > 0)
//	{
//		if( nifti_add_exten_to_list(&extn, &Elist, count+1) < 0 ){
//			if( g_opts.debug > 0 )
//				fprintf(stderr,"** failed adding ext %ld to list\n", count);
//			return -1;
//		}
//		
//		/* we have a new extension */
//		if( g_opts.debug > 1 ){
//			fprintf(stderr,"+d found extension #%ld, code = 0x%x, size = %d\n",
//					count, extn.ecode, extn.esize);
//			if( extn.ecode == NIFTI_ECODE_AFNI && g_opts.debug > 2 ) /* ~XML */
//				fprintf(stderr,"   AFNI extension: %.*s\n",
//						extn.esize-8,extn.edata);
//			else if( extn.ecode == NIFTI_ECODE_COMMENT && g_opts.debug > 2 )
//				fprintf(stderr,"   COMMENT extension: %.*s\n",        /* TEXT */
//						extn.esize-8,extn.edata);
//		}
//		remain -= extn.esize;
//		count++;
//	}
//	
//	if( g_opts.debug > 2 ) fprintf(stderr,"+d found %ld extension(s)\n", count);
//	
//	nim->num_ext = (int)count;
//	nim->ext_list = Elist;
//	
//	return (int)count;
//}
//
//
///*----------------------------------------------------------------------*/
///*! nifti_add_extension - add an extension, with a copy of the data
// 
// Add an extension to the nim->ext_list array.
// Fill this extension with a copy of the data, noting the
// length and extension code.
// 
// \param nim    - nifti_image to add extension to
// \param data   - raw extension data
// \param length - length of raw extension data
// \param ecode  - extension code
// 
// \sa extension codes NIFTI_ECODE_* in nifti1_io.h
// \sa nifti_free_extensions, valid_nifti_extensions, nifti_copy_extensions
// 
// \return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//int nifti_add_extension(nifti_image *nim, const char * data, int len, int ecode)
//{
//	nifti1_extension ext;
//	
//	/* error are printed in functions */
//	if( nifti_fill_extension(&ext, data, len, ecode) )                 return -1;
//	if( nifti_add_exten_to_list(&ext, &nim->ext_list, nim->num_ext+1)) return -1;
//	
//	nim->num_ext++;  /* success, so increment */
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------*/
///* nifti_add_exten_to_list     - add a new nifti1_extension to the list
// 
// We will append via "malloc, copy and free", because on an error,
// the list will revert to the previous one (sorry realloc(), only
// quality dolphins get to become part of St@rk!st brand tunafish).
// 
// return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//static int nifti_add_exten_to_list( nifti1_extension *  new_ext,
//								   nifti1_extension ** list, long new_length )
//{
//	nifti1_extension * tmplist;
//	
//	tmplist = *list;
//	*list = (nifti1_extension *)malloc(new_length * sizeof(nifti1_extension));
//	
//	/* check for failure first */
//	if( ! *list ){
//		fprintf(stderr,"** failed to alloc %ld extension structs (%ld bytes)\n",
//				new_length, new_length*(int)sizeof(nifti1_extension));
//		if( !tmplist ) return -1;  /* no old list to lose */
//		
//		*list = tmplist;  /* reset list to old one */
//		return -1;
//	}
//	
//	/* if an old list exists, copy the pointers and free the list */
//	if( tmplist ){
//		memcpy(*list, tmplist, (new_length-1)*sizeof(nifti1_extension));
//		free(tmplist);
//	}
//	
//	/* for some reason, I just don't like struct copy... */
//	(*list)[new_length-1].esize = new_ext->esize;
//	(*list)[new_length-1].ecode = new_ext->ecode;
//	(*list)[new_length-1].edata = new_ext->edata;
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d allocated and appended extension #%ld to list\n",
//				new_length);
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------*/
///* nifti_fill_extension  - given data and length, fill an extension struct
// 
// Allocate memory for data, copy data, set the size and code.
// 
// return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//static int nifti_fill_extension( nifti1_extension *ext, const char * data,
//                                int len, int ecode)
//{
//	int esize;
//	
//	if( !ext || !data || len < 0 ){
//		fprintf(stderr,"** fill_ext: bad params (%p,%p,%d)\n",
//				(void *)ext, data, len);
//		return -1;
//	} else if( ! nifti_is_valid_ecode(ecode) ){
//		fprintf(stderr,"** fill_ext: invalid ecode %d\n", ecode);
//		return -1;
//	}
//	
//	/* compute esize, first : len+8, and take ceiling up to a mult of 16 */
//	esize = len+8;
//	if( esize & 0xf ) esize = (esize + 0xf) & ~0xf;
//	ext->esize = esize;
//	
//	/* allocate esize-8 (maybe more than len), using calloc for fill */
//	ext->edata = (char *)calloc(esize-8, sizeof(char));
//	if( !ext->edata ){
//		fprintf(stderr,"** NFE: failed to alloc %d bytes for extension\n",len);
//		return -1;
//	}
//	
//	memcpy(ext->edata, data, len);  /* copy the data, using len */
//	ext->ecode = ecode;             /* set the ecode */
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d alloc %d bytes for ext len %d, ecode %d, esize %d\n",
//				esize-8, len, ecode, esize);
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------
// * nifti_read_next_extension  - read a single extension from the file
// *
// * return (>= 0 is okay):
// *
// *     success      : esize
// *     no extension : 0
// *     error        : -1
// *----------------------------------------------------------------------*/
//static int nifti_read_next_extension( nifti1_extension * nex, nifti_image *nim,
//									 int remain, znzFile fp )
//{
//	int swap = nim->byteorder != nifti_short_order();
//	int count, size, code;
//	
//	/* first clear nex */
//	nex->esize = nex->ecode = 0;
//	nex->edata = NULL;
//	
//	if( remain < 16 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d only %d bytes remain, so no extension\n", remain);
//		return 0;
//	}
//	
//	/* must start with 4-byte size and code */
//	count = (int)znzread( &size, 4, 1, fp );
//	if( count == 1 ) count += (int)znzread( &code, 4, 1, fp );
//	
//	if( count != 2 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d current extension read failed\n");
//		znzseek(fp, -4*count, SEEK_CUR); /* back up past any read */
//		return 0;                        /* no extension, no error condition */
//	}
//	
//	if( swap ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d pre-swap exts: code %d, size %d\n", code, size);
//		
//		nifti_swap_4bytes(1, &size);
//		nifti_swap_4bytes(1, &code);
//	}
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d potential extension: code %d, size %d\n", code, size);
//	
//	if( !nifti_check_extension(nim, size, code, remain) ){
//		if( znzseek(fp, -8, SEEK_CUR) < 0 ){      /* back up past any read */
//			fprintf(stderr,"** failure to back out of extension read!\n");
//			return -1;
//		}
//		return 0;
//	}
//	
//	/* now get the actual data */
//	nex->esize = size;
//	nex->ecode = code;
//	
//	size -= 8;  /* subtract space for size and code in extension */
//	nex->edata = (char *)malloc(size * sizeof(char));
//	if( !nex->edata ){
//		fprintf(stderr,"** failed to allocate %d bytes for extension\n",size);
//		return -1;
//	}
//	
//	count = (int)znzread(nex->edata, 1, size, fp);
//	if( count < size ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"-d read only %d (of %d) bytes for extension\n",
//					count, size);
//		free(nex->edata);
//		nex->edata = NULL;
//		return -1;
//	}
//	
//	/* success! */
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d successfully read extension, code %d, size %d\n",
//				nex->ecode, nex->esize);
//	
//	return nex->esize;
//}
//
//
///*----------------------------------------------------------------------*/
///*! for each extension, check code, size and data pointer
// *//*--------------------------------------------------------------------*/
//int valid_nifti_extensions(const nifti_image * nim)
//{
//	nifti1_extension * ext;
//	int                c, errs;
//	
//	if( nim->num_ext <= 0 || nim->ext_list == NULL ){
//		if( g_opts.debug > 2 ) fprintf(stderr,"-d empty extension list\n");
//		return 0;
//	}
//	
//	/* for each extension, check code, size and data pointer */
//	ext = nim->ext_list;
//	errs = 0;
//	for ( c = 0; c < nim->num_ext; c++ ){
//		if( ! nifti_is_valid_ecode(ext->ecode) ) {
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, invalid code %d\n", c, ext->ecode);
//			errs++;
//		}
//		
//		if( ext->esize <= 0 ){
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, bad size = %d\n", c, ext->esize);
//			errs++;
//		} else if( ext->esize & 0xf ){
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, size %d not multiple of 16\n",
//						c, ext->esize);
//			errs++;
//		}
//		
//		if( ext->edata == NULL ){
//			if( g_opts.debug > 1 ) fprintf(stderr,"-d ext %d, missing data\n", c);
//			errs++;
//		}
//		
//		ext++;
//	}
//	
//	if( errs > 0 ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"-d had %d extension errors, none will be written\n",
//					errs);
//		return 0;
//	}
//	
//	/* if we're here, we're good */
//	return 1;
//}
//
//
///*----------------------------------------------------------------------*/
///*! check whether the extension code is valid
// 
// \return 1 if valid, 0 otherwise
// *//*--------------------------------------------------------------------*/
//int nifti_is_valid_ecode( int ecode )
//{
//	if( ecode < NIFTI_ECODE_IGNORE  ||   /* minimum code number (0) */
//       ecode > NIFTI_MAX_ECODE     ||   /* maximum code number     */
//       ecode & 1 )                      /* cannot be odd           */
//		return 0;
//	
//	return 1;
//}
//
//
///*----------------------------------------------------------------------
// * check for valid size and code, as well as can be done
// *----------------------------------------------------------------------*/
//static int nifti_check_extension(nifti_image *nim, int size, int code, int rem)
//{
//	/* check for bad code before bad size */
//	if( ! nifti_is_valid_ecode(code) ) {
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d invalid extension code %d\n",code);
//		return 0;
//	}
//	
//	if( size < 16 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d ext size %d, no extension\n",size);
//		return 0;
//	}
//	
//	if( size > rem ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d ext size %d, space %d, no extension\n", size, rem);
//		return 0;
//	}
//	
//	if( size & 0xf ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d nifti extension size %d not multiple of 16\n",size);
//		return 0;
//	}
//	
//	if( nim->nifti_type == NIFTI_FTYPE_ASCII && size > LNI_MAX_NIA_EXT_LEN ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d NVE, bad nifti_type 3 size %d\n", size);
//		return 0;
//	}
//	
//	return 1;
//}
//
///*--------------------------------------------------------------------------*/
///*! free the nifti extensions
// 
// - If any edata pointer is set in the extension list, free() it.
// - Free ext_list, if it is set.
// - Clear num_ext and ext_list from nim.
// 
// \return 0 on success, -1 on error
// 
// \sa nifti_add_extension, nifti_copy_extensions
// *//*------------------------------------------------------------------------*/
//int nifti_free_extensions( nifti_image *nim )
//{
//	int c ;
//	if( nim == NULL ) return -1;
//	if( nim->num_ext > 0 && nim->ext_list ){
//		for( c = 0; c < nim->num_ext; c++ )
//			if ( nim->ext_list[c].edata ) free(nim->ext_list[c].edata);
//		free(nim->ext_list);
//	}
//	/* or if it is inconsistent, warn the user (if we are not in quiet mode) */
//	else if ( (nim->num_ext > 0 || nim->ext_list != NULL) && (g_opts.debug > 0) )
//		fprintf(stderr,"** warning: nifti extension num/ptr mismatch (%d,%p)\n",
//				nim->num_ext, (void *)nim->ext_list);
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d free'd %d extension(s)\n", nim->num_ext);
//	
//	nim->num_ext = 0;
//	nim->ext_list = NULL;
//	
//	return 0;
//}
///* return number of extensions written, or -1 on error */
//static int nifti_write_extensions(znzFile fp, nifti_image *nim)
//{
//	nifti1_extension * list;
//	char               extdr[4] = { 0, 0, 0, 0 };
//	int                c, size, ok = 1;
//	
//	if( znz_isnull(fp) || !nim || nim->num_ext < 0 ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"** nifti_write_extensions, bad params\n");
//		return -1;
//	}
//	
//	/* if no extensions and user requests it, skip extender */
//	if( g_opts.skip_blank_ext && (nim->num_ext == 0 || ! nim->ext_list ) ){
//		if( g_opts.debug > 1 )
//			fprintf(stderr,"-d no exts and skip_blank_ext set, "
//					"so skipping 4-byte extender\n");
//		return 0;
//	}
//	
//	/* if invalid extension list, clear num_ext */
//	if( ! valid_nifti_extensions(nim) ) nim->num_ext = 0;
//	
//	/* write out extender block */
//	if( nim->num_ext > 0 ) extdr[0] = 1;
//	if( nifti_write_buffer(fp, extdr, 4) != 4 ){
//		fprintf(stderr,"** failed to write extender\n");
//		return -1;
//	}
//	
//	list = nim->ext_list;
//	for ( c = 0; c < nim->num_ext; c++ ){
//		size = (int)nifti_write_buffer(fp, &list->esize, sizeof(int));
//		ok = (size == (int)sizeof(int));
//		if( ok ){
//			size = (int)nifti_write_buffer(fp, &list->ecode, sizeof(int));
//			ok = (size == (int)sizeof(int));
//		}
//		if( ok ){
//			size = (int)nifti_write_buffer(fp, list->edata, list->esize - 8);
//			ok = (size == list->esize - 8);
//		}
//		
//		if( !ok ){
//			fprintf(stderr,"** failed while writing extension #%d\n",c);
//			return -1;
//		} else if ( g_opts.debug > 2 )
//			fprintf(stderr,"+d wrote extension %d of %d bytes\n", c, size);
//		
//		list++;
//	}
//	
//	if( g_opts.debug > 1 )
//		fprintf(stderr,"+d wrote out %d extension(s)\n", nim->num_ext);
//	
//	return nim->num_ext;
//}
///*----------------------------------------------------------------------*/
///*! \fn int nifti_copy_extensions(nifti_image * nim_dest, nifti_image * nim_src)
// \brief copy the nifti1_extension list from src to dest
// 
// Duplicate the list of nifti1_extensions.  The dest structure must
// be clear of extensions.
// \return 0 on success, -1 on failure
// 
// \sa nifti_add_extension, nifti_free_extensions
// */
//int nifti_copy_extensions(nifti_image * nim_dest, const nifti_image * nim_src)
//{
//	char   * data;
//	size_t   bytes;
//	int      c, size, old_size;
//	
//	if( nim_dest->num_ext > 0 || nim_dest->ext_list != NULL ){
//		fprintf(stderr,"** will not copy extensions over existing ones\n");
//		return -1;
//	}
//	
//	if( g_opts.debug > 1 )
//		fprintf(stderr,"+d duplicating %d extension(s)\n", nim_src->num_ext);
//	
//	if( nim_src->num_ext <= 0 ) return 0;
//	
//	bytes = nim_src->num_ext * sizeof(nifti1_extension);  /* I'm lazy */
//	nim_dest->ext_list = (nifti1_extension *)malloc(bytes);
//	if( !nim_dest->ext_list ){
//		fprintf(stderr,"** failed to allocate %d nifti1_extension structs\n",
//				nim_src->num_ext);
//		return -1;
//	}
//	
//	/* copy the extension data */
//	nim_dest->num_ext = 0;
//	for( c = 0; c < nim_src->num_ext; c++ ){
//		size = old_size = nim_src->ext_list[c].esize;
//		if( size & 0xf ) size = (size + 0xf) & ~0xf; /* make multiple of 16 */
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"+d dup'ing ext #%d of size %d (from size %d)\n",
//					c, size, old_size);
//		/* data length is size-8, as esize includes space for esize and ecode */
//		data = (char *)calloc(size-8,sizeof(char));      /* maybe size > old */
//		if( !data ){
//			fprintf(stderr,"** failed to alloc %d bytes for extention\n", size);
//			if( c == 0 ) { free(nim_dest->ext_list); nim_dest->ext_list = NULL; }
//			/* otherwise, keep what we have (a.o.t. deleting them all) */
//			return -1;
//		}
//		/* finally, fill the new structure */
//		nim_dest->ext_list[c].esize = size;
//		nim_dest->ext_list[c].ecode = nim_src->ext_list[c].ecode;
//		nim_dest->ext_list[c].edata = data;
//		memcpy(data, nim_src->ext_list[c].edata, old_size-8);
//		
//		nim_dest->num_ext++;
//	}
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------*/
///*! compute the total size of all extensions
// 
// \return the total of all esize fields
// 
// Note that each esize includes 4 bytes for ecode, 4 bytes for esize,
// and the bytes used for the data.  Each esize also needs to be a
// multiple of 16, so it may be greater than the sum of its 3 parts.
// *//*--------------------------------------------------------------------*/
//int nifti_extension_size(nifti_image *nim)
//{
//	int c, size = 0;
//	
//	if( !nim || nim->num_ext <= 0 ) return 0;
//	
//	if( g_opts.debug > 2 ) fprintf(stderr,"-d ext sizes:");
//	
//	for ( c = 0; c < nim->num_ext; c++ ){
//		size += nim->ext_list[c].esize;
//		if( g_opts.debug > 2 ) fprintf(stderr,"  %d",nim->ext_list[c].esize);
//	}
//	
//	if( g_opts.debug > 2 ) fprintf(stderr," (total = %d)\n",size);
//	
//	return size;
//}

/*----------------------------------------------------------------------*/
/*! get the byte order for this CPU
 
 - LSB_FIRST means least significant byte, first (little endian)
 - MSB_FIRST means most significant byte, first (big endian)
 *//*--------------------------------------------------------------------*/
int nifti_short_order(void)   /* determine this CPU's byte order */
{
	union { unsigned char bb[2] ;
		short         ss    ; } fred ;
	
	fred.bb[0] = 1 ; fred.bb[1] = 0 ;
	
	return (fred.ss == 1) ? LSB_FIRST : MSB_FIRST ;
}

NiftiImage::NiftiImage() :
	_dim(),
	_voxdim(),
	_mode(NIFTI_CLOSED),
	_gz(false)
{
	_qform.setIdentity(); _sform.setIdentity();
}

NiftiImage::NiftiImage(const NiftiImage &clone) :
	_nvox(clone._nvox),
	_qform(clone._qform), _sform(clone._sform), _inverse(clone._inverse),
	_datatype(clone._datatype),	_mode(NIFTI_CLOSED), _gz(false),
	_voxoffset(0),
	scaling_slope(clone.scaling_slope), scaling_inter(clone.scaling_inter),
	calibration_min(clone.calibration_min), calibration_max(clone.calibration_max),
	qform_code(clone.qform_code), sform_code(clone.sform_code),
	freq_dim(clone.freq_dim), phase_dim(clone.phase_dim),
	slice_dim(clone.slice_dim), slice_code(clone.slice_code),
	slice_start(clone.slice_start), slice_end(clone.slice_end),
	slice_duration(clone.slice_duration),
	toffset(clone.toffset),
	xyz_units(clone.xyz_units), time_units(clone.time_units),
	intent_code(clone.intent_code),
	intent_p1(clone.intent_p1), intent_p2(clone.intent_p2), intent_p3(clone.intent_p3),
	intent_name(clone.intent_name),
	description(clone.description),
	aux_file(clone.aux_file)
{
	for (int i = 0; i < 8; i++)
	{
		_dim[i] = clone._dim[i];
		_voxdim[i] = clone._voxdim[i];
	}
	//int                num_ext ;  /*!< number of extensions in ext_list       */
	//nifti1_extension * ext_list ; /*!< array of extension structs (with data) */
	//analyze_75_orient_code analyze75_orient; /*!< for old analyze files, orient */
}

bool isGZippedFile(const std::string &fname)
{
	if (fname.find_last_of(".") != std::string::npos)
	{
		std::string ext = fname.substr(fname.find_last_of(".") + 1);
		std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
		if (ext == "gz")
			return true;
	}
	return false;
}

void NiftiImage::setFilenames(const std::string &fname)
{	
	std::string ext = fname.substr(fname.find_last_of(".") + 1);
	_basename = fname.substr(0, fname.find_last_of("."));
	_gz = false;
	if (ext == "gz")
	{
		_gz = true;
		ext = _basename.substr(_basename.find_last_of(".") + 1);
		_basename = _basename.substr(0, _basename.find_last_of("."));
	}
	if (ext == "hdr" || ext == "img")
	{
		_imgname = _basename + ".img";
		_hdrname = _basename + ".hdr";
	}
	else if (ext == "nii")
	{
		_imgname = _hdrname = _basename + ".nii";
	}
	else
	{
		std::cerr << "Extension " << ext << " is not a valid NIfTI extension." << std::endl;
		abort();
	}
	if (_gz)
	{
		_imgname += ".gz";
		_hdrname += ".gz";
	}
}

const std::string &NiftiImage::basename() { return _basename; }

/* we already assume ints are 4 bytes */
#undef ZNZ_MAX_BLOCK_SIZE
#define ZNZ_MAX_BLOCK_SIZE (1<<30)

size_t NiftiImage::read(void *buff, size_t size, size_t nmemb)
{
	if (_gz)
	{
		unsigned long remain = size*nmemb;
		char *cbuf = (char *)buff;
		unsigned long n2read;
		int nread;
		/* gzread/write take unsigned int length, so maybe read in int pieces
		 (noted by M Hanke, example given by M Adler)   6 July 2010 [rickr] */
		while(remain > 0) {
			n2read = (remain < ZNZ_MAX_BLOCK_SIZE) ? remain : ZNZ_MAX_BLOCK_SIZE;
			nread = gzread(_gzFile, (void *)cbuf, (unsigned int)n2read);
			if( nread < 0 ) return nread; /* returns -1 on error */
			
			remain -= nread;
			cbuf += nread;
			
			/* require reading n2read bytes, so we don't get stuck */
			if( nread < (int)n2read ) break;  /* return will be short */
		}
		
		/* warn of a short read that will seem complete */
		if( remain > 0 && remain < size )
			std::cerr << "NiftiImage: Zipped read short by " << remain << " bytes." << std::endl;
		return nmemb - remain/size;   /* return number of members processed */
	}
	else
		return fread(buff, size, nmemb, _file);
}

size_t NiftiImage::write(const void *buff, size_t size, size_t nmemb)
{
	if (_gz) {
		unsigned long remain = size*nmemb;
		char     * cbuf = (char *)buff;
		unsigned long n2write;
		int        nwritten;
		while( remain > 0 ) {
			n2write = (remain < ZNZ_MAX_BLOCK_SIZE) ? remain : ZNZ_MAX_BLOCK_SIZE;
			nwritten = gzwrite(_gzFile, (void *)cbuf, (unsigned int)n2write);
			
			/* gzread returns 0 on error, but in case that ever changes... */
			if( nwritten < 0 ) return nwritten;
			
			remain -= nwritten;
			cbuf += nwritten;
			
			/* require writing n2write bytes, so we don't get stuck */
			if( nwritten < (int)n2write ) break;
		}
		
		/* warn of a short write that will seem complete */
		if( remain > 0 && remain < size )
			std::cerr << "NiftiImage: Zipped write short by " << remain << " bytes." << std::endl;
		
		return nmemb - remain/size;   /* return number of members processed */
	}
	return fwrite(buff, size, nmemb, _file);
}

long NiftiImage::seek(long offset, int whence)
{
	long newpos, currpos;
	if (_gz)
		currpos = gztell(_gzFile);
	else
		currpos = ftell(_file);
	if (_mode == NIFTI_READ) {
		if (whence == SEEK_SET)
			newpos = offset;
		else if (whence == SEEK_CUR)
			newpos = currpos + offset;
		else
		{
			std::cerr << "NiftiImage: Reading with SEEK_END is not implemented." << std::endl;
			abort();
		}
		
		if (newpos < _voxoffset)
		{
			std::cerr << "NiftiImage: Attempted to seek into the header while reading from " << _imgname << "." << std::endl;
			abort();
		}
	} else if (_mode == NIFTI_WRITE) {
		if (whence == SEEK_SET)
			newpos = offset;
		else if (whence == SEEK_CUR)
			newpos = currpos + offset;
		else
		{
			std::cerr << "NiftiImage: Writing with SEEK_END is not implemented." << std::endl;
			abort();
		}
		
		if (newpos < _voxoffset)
		{
			std::cerr << "NiftiImage: Attempted to seek into the header while writing to " << _imgname << "." << std::endl;
			abort();
		}
		else if (_gz && newpos < currpos)
		{
			std::cerr << "NiftiImage: Cannot seek backwards while writing a file with libz." << std::endl;
			abort();
		}
	} else {
		std::cerr << "NiftiImage: Cannot seek in a closed file." << std::endl;
		abort();
	}
	long error;
	if (_gz)
	{
		error = gzseek(_gzFile, offset, whence);
		currpos = gztell(_gzFile);
	}
	else
	{
		error = fseek(_file, offset, whence);
		currpos = ftell(_file);
	}
	if (currpos != newpos)
	{
		std::cerr << "NiftiImage: Seek moved to an unexpected location ("
				  << currpos << ", not " << newpos << ")." << std::endl;
	}
	return error;
}

/*----------------------------------------------------------------------
 * check whether byte swapping is needed
 *
 * dim[0] should be in [0,7], and sizeof_hdr should be accurate
 *
 * \returns  > 0 : needs swap
 *             0 : does not need swap
 *           < 0 : error condition
 *----------------------------------------------------------------------*/
int NiftiImage::needs_swap( short dim0, int hdrsize )
{
	short d0    = dim0;     /* so we won't have to swap them on the stack */
	int   hsize = hdrsize;
	
	if( d0 != 0 ){     /* then use it for the check */
		if( d0 > 0 && d0 <= 7 ) return 0;
		
		nifti_swap_bytes(1, 2, &d0);        /* swap? */
		if( d0 > 0 && d0 <= 7 ) return 1;
		
		return -1;        /* bad, naughty d0 */
	}
	/* dim[0] == 0 should not happen, but could, so try hdrsize */
	if( hsize == sizeof(nifti_1_header) ) return 0;
	nifti_swap_bytes(1, 4, &hsize);     /* swap? */
	if( hsize == sizeof(nifti_1_header) ) return 1;
	
	return -2;     /* bad, naughty hsize */
}

inline float NiftiImage::fixFloat(const float f)
{
	if (std::isfinite(f))
		return f;
	else
		return 0.;
}

void NiftiImage::readHeader(std::string path)
{
	struct nifti_1_header  nhdr;
	size_t                 bytes_read;

	if (_gz)
	{
		_gzFile = gzopen(path.c_str(), "rb");
		_file = NULL;
	}
	else
	{
		_file = fopen(path.c_str(), "rb");
		_gzFile = NULL;
	}
	if (!(_gzFile || _file))
	{
		std::cerr << "Failed to open header from file: " << path << std::endl;
		abort();
	}
	
	if (_gz)
		bytes_read = gzread(_gzFile, &nhdr, sizeof(nhdr));  /* read the thing */
	else
		bytes_read = fread(&nhdr, sizeof(nhdr), 1, _file);
	
	/* keep file open so we can check for exts. after nifti_convert_nhdr2nim() */
	if (bytes_read < sizeof(nhdr))
	{
		std::cerr << "Only read " << bytes_read << " bytes from " << _hdrname <<
		             ", should have been " << sizeof(nhdr) << "." << std::endl;
		abort();
	}
	
	/**- check if we must swap bytes */
	int doswap = needs_swap(nhdr.dim[0], nhdr.sizeof_hdr); /* swap data flag */
	if(doswap < 0)
	{
		std::cerr << "Could not determine byte order of header " <<
					 _hdrname << "." << std::endl;
		abort();
	}
	
	// Check the magic string is set to one of the possible NIfTI values,
	// otherwise process as an ANALYZE file
	int is_nifti = ((nhdr.magic[0]=='n' && nhdr.magic[3]=='\0') && 
                    (nhdr.magic[1]=='i' || nhdr.magic[1]=='+') &&
                    (nhdr.magic[2]>='1' && nhdr.magic[2]<='9'))
				   ? nhdr.magic[2]-'0' : 0;
	if (doswap)
		swap_nifti_header(&nhdr, is_nifti);
	
	if(nhdr.datatype == DT_BINARY || nhdr.datatype == DT_UNKNOWN  )
		std::cerr << "Bad datatype in header " << _hdrname << std::endl;
	if(nhdr.dim[1] <= 0)
		std::cerr << "Bad first dimension in header " << _hdrname << std::endl;
	
	/* fix bad dim[] values in the defined dimension range */
	for(int i=2; i <= nhdr.dim[0]; i++)
		if(nhdr.dim[i] <= 0) nhdr.dim[i] = 1;
	
	/* fix any remaining bad dim[] values, so garbage does not propagate */
	/* (only values 0 or 1 seem rational, otherwise set to arbirary 1)   */
	for(int i=nhdr.dim[0]+1; i < 8; i++)
		if(nhdr.dim[i] != 1 && nhdr.dim[i] != 0) nhdr.dim[i] = 1 ;
	
	_byteorder = nifti_short_order();
	if (doswap)
		_byteorder = REVERSE_ORDER(_byteorder);
	
	/**- Set dimensions arrays */
	for (int i = 0; i < 8; i++)
	{
		_dim[i] = nhdr.dim[i];
		_voxdim[i] = nhdr.pixdim[i];
	}
	
	_nvox = 1;
	for(int i=1; i <= _dim[0]; i++ )
		_nvox *= _dim[i];
	
	/**- set the type of data in voxels and how many bytes per voxel */
	_datatype = nhdr.datatype;
	
	/**- compute qto_xyz transformation from pixel indexes (i,j,k) to (x,y,z) */
	Affine3d S; S = Scaling<double>(_voxdim[1], _voxdim[2], _voxdim[3]);
	if( !is_nifti || nhdr.qform_code <= 0 ){
		/**- if not nifti or qform_code <= 0, use grid spacing for qto_xyz */
		_qform = S.matrix();
		qform_code = NIFTI_XFORM_UNKNOWN ;
	} else {
		/**- else NIFTI: use the quaternion-specified transformation */
		double b = fixFloat(nhdr.quatern_b);
		double c = fixFloat(nhdr.quatern_c);
		double d = fixFloat(nhdr.quatern_d);
		
		double x = fixFloat(nhdr.qoffset_x);
		double y = fixFloat(nhdr.qoffset_y);
		double z = fixFloat(nhdr.qoffset_z);
		double qfac = (nhdr.pixdim[0] < 0.0) ? -1.0 : 1.0 ;  /* left-handedness? */
		double a = sqrt(1 - (b*b + c*c + d*d));
		Quaterniond Q(a, b, c, d);
		Affine3d T; T = Translation3d(x, y, z);
		_qform = T*S*Q;
		_inverse = _qform.inverse();
		if (qfac < 0.)
			_qform.matrix().block(0, 2, 3, 1) *= -1.;
		qform_code = nhdr.qform_code;
	}
	/**- load sto_xyz affine transformation, if present */
	if( !is_nifti || nhdr.sform_code <= 0 )
	{	/**- if not nifti or sform_code <= 0, then no sto transformation */
		sform_code = NIFTI_XFORM_UNKNOWN ;
	} else {
		/**- else set the sto transformation from srow_*[] */
		_sform.setIdentity();
		for (int i = 0; i < 4; i++)
		{
			_sform(0, i) = nhdr.srow_x[i];
			_sform(1, i) = nhdr.srow_y[i];
			_sform(2, i) = nhdr.srow_z[i];
		}
		_inverse = _sform.inverse();
		sform_code = nhdr.sform_code ;
	}
	
	/**- set miscellaneous NIFTI stuff */
	if (is_nifti) {
		scaling_slope   = fixFloat(nhdr.scl_slope);
		if (scaling_slope == 0.)
			scaling_slope = 1.;
		scaling_inter   = fixFloat(nhdr.scl_inter);
		intent_code = nhdr.intent_code;
		intent_p1 = fixFloat(nhdr.intent_p1);
		intent_p2 = fixFloat(nhdr.intent_p2);
		intent_p3 = fixFloat(nhdr.intent_p3);
		toffset   = fixFloat(nhdr.toffset);
		
		intent_name = std::string(nhdr.intent_name);
		xyz_units = XYZT_TO_SPACE(nhdr.xyzt_units);
		time_units = XYZT_TO_TIME(nhdr.xyzt_units);
		freq_dim  = DIM_INFO_TO_FREQ_DIM (nhdr.dim_info);
		phase_dim = DIM_INFO_TO_PHASE_DIM(nhdr.dim_info);
		slice_dim = DIM_INFO_TO_SLICE_DIM(nhdr.dim_info);
		slice_code     = nhdr.slice_code;
		slice_start    = nhdr.slice_start;
		slice_end      = nhdr.slice_end;
		slice_duration = fixFloat(nhdr.slice_duration);
	}
	
	/**- set Miscellaneous ANALYZE stuff */
	
	calibration_min = fixFloat(nhdr.cal_min);
	calibration_max = fixFloat(nhdr.cal_max);
	
	description = std::string(nhdr.descrip);
	aux_file    = std::string(nhdr.aux_file);
	
	/**- set ioff from vox_offset (but at least sizeof(header)) */
	if (_hdrname == _imgname) {
		_voxoffset = (int)nhdr.vox_offset;
		if (_voxoffset < (int)sizeof(nhdr)) _voxoffset = (int)sizeof(nhdr);
	} else {
		_voxoffset = (int)nhdr.vox_offset ;
	}
	
	/* clear extension fields */
	_num_ext = 0;
	_ext_list = NULL;
	/**- check for extensions (any errors here means no extensions) */
	//if( NIFTI_ONEFILE(nhdr) ) remaining = nim->iname_offset - sizeof(nhdr);
	//else                      remaining = filesize - sizeof(nhdr);
	
	//(void)nifti_read_extensions(nim, fp, remaining);
	
	if (_hdrname != _imgname)
	{	// Need to close the header and open the image
		if (_gz)
		{
			gzclose(_gzFile);
			_gzFile = gzopen(_imgname.c_str(), "rb");
		}
		else
		{
			fclose(_file);
			_file = fopen(_imgname.c_str(), "rb");
		}
	}
	// Now we have a valid transform, store the inverse for future
}

void NiftiImage::writeHeader(std::string path)
{
	struct nifti_1_header nhdr;
	memset(&nhdr,0,sizeof(nhdr)) ;  /* zero out header, to be safe */
	/**- load the ANALYZE-7.5 generic parts of the header struct */
	nhdr.sizeof_hdr = sizeof(nhdr);
	nhdr.regular    = 'r';             /* for some stupid reason */
	
	for (int i = 0; i < 8; i++)
	{	// Copy this way so types can be changed
		nhdr.dim[i] = _dim[i];
		nhdr.pixdim[i] = _voxdim[i];
	}
	
	nhdr.datatype = _datatype;
	nhdr.bitpix   = 8 * DTypes.find(_datatype)->second.size;
	
	if(calibration_max > calibration_min) {
		nhdr.cal_max = calibration_max;
		nhdr.cal_min = calibration_min;
	}
	
	if(scaling_slope != 0.0) {
		nhdr.scl_slope = scaling_slope;
		nhdr.scl_inter = scaling_inter;
	}
	
	strncpy(nhdr.descrip, description.c_str(), 80);
	strncpy(nhdr.aux_file, aux_file.c_str(), 24);
	

	if(_imgname == _hdrname)
		strcpy(nhdr.magic,"n+1");
	else
		strcpy(nhdr.magic,"ni1");
	for (int i = 1; i < 8; i++)
		nhdr.pixdim[i] = fabs(nhdr.pixdim[i]);
	
	nhdr.intent_code = intent_code;
	nhdr.intent_p1   = intent_p1;
	nhdr.intent_p2   = intent_p2;
	nhdr.intent_p3   = intent_p3;
	strncpy(nhdr.intent_name, intent_name.c_str(), 16);
	
	// Check that _voxoffset is sensible
	if (_imgname == _hdrname && _voxoffset < nhdr.sizeof_hdr)
		_voxoffset = 352;
	nhdr.vox_offset = _voxoffset ;
	nhdr.xyzt_units = SPACE_TIME_TO_XYZT(xyz_units, time_units);
	nhdr.toffset    = toffset ;
	
	if(qform_code > 0) {
		nhdr.qform_code = qform_code ;
		Quaterniond Q(_qform.rotation());
		Translation3d T(_qform.translation());
		Affine3d S; S = Scaling<double>(_voxdim[1], _voxdim[2], _voxdim[3]);

		// NIfTI REQUIRES a (or w) >= 0. Because Q and -Q represent the same
		// rotation, if w < 0 simply store -Q
		if (Q.w() < 0)
		{
			nhdr.quatern_b  = -Q.x();
			nhdr.quatern_c  = -Q.y();
			nhdr.quatern_d  = -Q.z();
		}
		else
		{
			nhdr.quatern_b = Q.x();
			nhdr.quatern_c = Q.y();
			nhdr.quatern_d = Q.z();
		}
		nhdr.qoffset_x  = T.x();
		nhdr.qoffset_y  = T.y();
		nhdr.qoffset_z  = T.z();
	}
	
	if(sform_code > 0) {
		nhdr.sform_code = sform_code;
		for (int i = 0; i < 4; i++)
		{
			nhdr.srow_x[i]  = _sform(0, i);
			nhdr.srow_y[i]  = _sform(1, i);
			nhdr.srow_z[i]  = _sform(2, i);
		}
	}
	
	nhdr.dim_info = FPS_INTO_DIM_INFO(freq_dim, phase_dim, slice_dim);
	nhdr.slice_code     = slice_code;
	nhdr.slice_start    = slice_start;
	nhdr.slice_end      = slice_end;
	nhdr.slice_duration = slice_duration;
	
	if (_gz)
	{
		_gzFile = gzopen(_hdrname.c_str(), "wb");
		_file = NULL;
	}
	else
	{
		_file =	fopen(_hdrname.c_str(), "wb");
		_gzFile = NULL;
	}
		
	if(!(_gzFile || _file)) {
		std::cerr << "NiftiImage: Cannot open header file " << _hdrname << " for writing." << std::endl;
		abort();
	}
	
	/* write the header and extensions */
	size_t bytesWritten;
	if (_gz)
		bytesWritten = gzwrite(_gzFile, &nhdr, sizeof(nhdr));
	else
		bytesWritten = fwrite(&nhdr, sizeof(nhdr), 1, _file);
	/* partial file exists, and errors have been printed, so ignore return */
	//if( nim->nifti_type != NIFTI_FTYPE_ANALYZE )
	//	(void)nifti_write_extensions(fp,nim);
	if(bytesWritten < sizeof(nhdr)) {
		std::cerr << "NiftiImage: Could not write header to file " << _hdrname << "." << std::endl;
		abort();
	}
	if (_hdrname != _imgname)
	{	// Close header and open image file
		if (_gz)
		{
			gzclose(_gzFile);
			_gzFile = gzopen(_imgname.c_str(), "wb");
		}
		else
		{
			fclose(_file);
			_file = fopen(_imgname.c_str(), "wb");
		}
		if (!(_file || _gzFile))
		{
			std::cerr << "Could not open image file " << _imgname << " for writing." << std::endl;
			abort();
		}
	}
}

void *NiftiImage::readBuffer(size_t start, size_t length)
{
	if (_mode == NIFTI_CLOSED)
	{
		std::cerr << "NiftiImage: Cannot read from a closed file." << std::endl;
		return NULL;
	}
	if (_mode == NIFTI_WRITE)
	{
		std::cerr << "NiftiImage: Cannot read from a file opened for writing." << std::endl;
		return NULL;
	}
	
	void *raw = malloc(length);
	seek(_voxoffset + start, SEEK_SET);
	size_t bytesRead;
	if (_gz)
		bytesRead = gzread(_gzFile, raw, static_cast<unsigned int>(length));
	else
		bytesRead = fread(raw, length, 1, _file);
	
	if (bytesRead < length)
	{
		std::cerr << "NiftiImage: Read buffer was short by " << (length - bytesRead)
				  << " bytes." << std::endl;
		free(raw);
		return NULL;
	}
	int swapsize = DTypes.find(_datatype)->second.swapsize;
	if (swapsize > 1 && _byteorder != nifti_short_order())
		nifti_swap_bytes(length / swapsize, swapsize, raw);
	return raw;
}

void NiftiImage::writeBuffer(void *data, size_t start, size_t length)
{
	if (_mode == NIFTI_CLOSED)
	{
		std::cerr << "NiftiImage: Cannot write to a closed file." << std::endl;
		return;
	}
	if (_mode == NIFTI_READ)
	{
		std::cerr << "NiftiImage: Cannot write to a file opened for writing." << std::endl;
		return;
	}
	seek(_voxoffset + start, SEEK_SET);
	size_t bytesWritten;
	if (_gz)
		bytesWritten = gzwrite(_gzFile, data, static_cast<unsigned int>(length));
	else
		bytesWritten = fwrite(data, length, 1, _file);
	
	if (bytesWritten < length)
	{
		std::cerr << "NiftiImage: Write buffer was short by " << (length - bytesWritten)
				  << " bytes." << std::endl;
	}
}

void *NiftiImage::readRawVolume(const int vol)
{
	size_t bytesPerVolume = voxelsPerVolume() *
	                        DTypes.find(_datatype)->second.size;
	void *raw = readBuffer(vol * bytesPerVolume, bytesPerVolume);
	return raw;
}
void *NiftiImage::readRawAllVolumes()
{
	void *raw =	readBuffer(0, nvox() * DTypes.find(_datatype)->second.swapsize);
	return raw;
}
		
void NiftiImage::open(std::string filename, char mode)
{
	setFilenames(filename);
	if (_mode != NIFTI_CLOSED)
	{
		std::cerr << "NiftiImage: Attempted to open file " << filename
				  << " using NiftiImage that is already open with file "
				  << _imgname << "." << std::endl;
		abort();
	}
	if (mode == NIFTI_READ) {
		_mode = NIFTI_READ;
		readHeader(_hdrname); // readHeader leaves _file pointing to image file
		seek(_voxoffset, SEEK_SET);
	} else if (mode == NIFTI_WRITE) {
		_mode = NIFTI_WRITE;
		writeHeader(_hdrname); // writeHeader ensures file is opened to the correct file
		seek(_voxoffset, SEEK_SET);
	} else {
		std::cerr << "Invalid NiftImage mode '" << mode << "'." << std::endl;
		abort();
	}
}

void NiftiImage::close()
{
	if (_gz)
	{
		gzflush(_gzFile, Z_FINISH);
		gzclose(_gzFile);
		_gzFile = NULL;
	}
	else
	{
		fflush(_file);
		fclose(_file);
		_file = NULL;
	}
	_mode = NIFTI_CLOSED;
}

NiftiImage::~NiftiImage()
{
	close();
}

int NiftiImage::ndim() const { return _dim[0]; }
int NiftiImage::nx() const { return _dim[1]; }
int NiftiImage::ny() const { return _dim[2]; }
int NiftiImage::nz() const { return _dim[3]; }
int NiftiImage::nt() const { return _dim[4]; }
void NiftiImage::setnt(const int nt)
{
	if (_mode == NIFTI_CLOSED)
	{
		if (nt > 1)
		{
			_dim[4] = nt;
			_dim[0] = 4;
		}
		else
		{
			_dim[4] = 1;
			_dim[0] = 3;
		}
	}
	else
	{
		std::cerr << "NiftiImage: Cannot change the dimensions of an image once opened." << std::endl;
		abort();
	}	
}

void NiftiImage::setDims(const int nx, const int ny, const int nz, const int nt)
{
	//setnt will check if the file is closed
	setnt(nt);
	_dim[1] = nx;
	_dim[2] = ny;
	_dim[3] = nz;
}
int NiftiImage::voxelsPerVolume() const { return _dim[1]*_dim[2]*_dim[3]; };
int NiftiImage::nvox() const { return _dim[1]*_dim[2]*_dim[3]*_dim[4]; };
float NiftiImage::dx() const { return _voxdim[1]; }
float NiftiImage::dy() const { return _voxdim[2]; }
float NiftiImage::dz() const { return _voxdim[3]; }

int NiftiImage::datatype() const { return _datatype; }
void NiftiImage::setDatatype(const int dt)
{
	if (_mode == NIFTI_READ)
	{
		std::cerr << "NiftiImage: Cannot set the datatype of a file opened for reading." << std::endl;
		return;
	}
	if (DTypeIsValid(dt))
		_datatype = dt;
	else
		std::cerr << "NiftiImage: Attempted to set invalid datatype " << dt << std::endl;
}

const Matrix4d &NiftiImage::qform() const { return _qform.matrix(); }
const Matrix4d &NiftiImage::sform() const { return _sform.matrix(); }
const Matrix4d &NiftiImage::ijk_to_xyz() const
{
	if ((sform_code > 0) && (sform_code >= qform_code))
		return _sform.matrix();
	else // There is always a _qform matrix
		return _qform.matrix();
}
const Matrix4d &NiftiImage::xyz_to_ijk() const
{
	return _inverse.matrix();
}

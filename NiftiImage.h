/** \file NiftiImage.h
 \brief Declaration for NiftiImage class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef _NIFTI_IO_HEADER_3
#define _NIFTI_IO_HEADER_3

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include <string>
#include <iostream>
#include <algorithm>
#include <complex>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "nifti1.h"                  /*** NIFTI-1 header specification ***/

using namespace Eigen;

/*! \enum analyze_75_orient_code
 *  \brief Old-style analyze75 orientation
 *         codes.
 */
typedef enum _analyze75_orient_code {
	a75_transverse_unflipped = 0,
	a75_coronal_unflipped = 1,
	a75_sagittal_unflipped = 2,
	a75_transverse_flipped = 3,
	a75_coronal_flipped = 4,
	a75_sagittal_flipped = 5,
	a75_orient_unknown = 6
} analyze_75_orient_code;

/* struct for return from nifti_image_read_bricks() */
typedef struct {
	int       nbricks;    /* the number of allocated pointers in 'bricks' */
	size_t    bsize;      /* the length of each data block, in bytes      */
	void   ** bricks;     /* array of pointers to data blocks             */
} nifti_brick_list;


/*****************************************************************************/
/*------------------ NIfTI version of ANALYZE 7.5 structure -----------------*/

/* (based on fsliolib/dbh.h, but updated for version 7.5) */

typedef struct {
	/* header info fields - describes the header    overlap with NIfTI */
	/*                                              ------------------ */
	int sizeof_hdr;                  /* 0 + 4        same              */
	char data_type[10];              /* 4 + 10       same              */
	char db_name[18];                /* 14 + 18      same              */
	int extents;                     /* 32 + 4       same              */
	short int session_error;         /* 36 + 2       same              */
	char regular;                    /* 38 + 1       same              */
	char hkey_un0;                   /* 39 + 1                40 bytes */
	
	/* image dimension fields - describes image sizes */
	short int dim[8];                /* 0 + 16       same              */
	short int unused8;               /* 16 + 2       intent_p1...      */
	short int unused9;               /* 18 + 2         ...             */
	short int unused10;              /* 20 + 2       intent_p2...      */
	short int unused11;              /* 22 + 2         ...             */
	short int unused12;              /* 24 + 2       intent_p3...      */
	short int unused13;              /* 26 + 2         ...             */
	short int unused14;              /* 28 + 2       intent_code       */
	short int datatype;              /* 30 + 2       same              */
	short int bitpix;                /* 32 + 2       same              */
	short int dim_un0;               /* 34 + 2       slice_start       */
	float pixdim[8];                 /* 36 + 32      same              */
	
	float vox_offset;                /* 68 + 4       same              */
	float funused1;                  /* 72 + 4       scl_slope         */
	float funused2;                  /* 76 + 4       scl_inter         */
	float funused3;                  /* 80 + 4       slice_end,        */
	/* slice_code,       */
	/* xyzt_units        */
	float cal_max;                   /* 84 + 4       same              */
	float cal_min;                   /* 88 + 4       same              */
	float compressed;                /* 92 + 4       slice_duration    */
	float verified;                  /* 96 + 4       toffset           */
	int glmax,glmin;                 /* 100 + 8              108 bytes */
	
	/* data history fields - optional */
	char descrip[80];                /* 0 + 80       same              */
	char aux_file[24];               /* 80 + 24      same              */
	char orient;                     /* 104 + 1      NO GOOD OVERLAP   */
	char originator[10];             /* 105 + 10     FROM HERE DOWN... */
	char generated[10];              /* 115 + 10                       */
	char scannum[10];                /* 125 + 10                       */
	char patient_id[10];             /* 135 + 10                       */
	char exp_date[10];               /* 145 + 10                       */
	char exp_time[10];               /* 155 + 10                       */
	char hist_un0[3];                /* 165 + 3                        */
	int views;                       /* 168 + 4                        */
	int vols_added;                  /* 172 + 4                        */
	int start_field;                 /* 176 + 4                        */
	int field_skip;                  /* 180 + 4                        */
	int omax, omin;                  /* 184 + 8                        */
	int smax, smin;                  /* 192 + 8              200 bytes */
} nifti_analyze75;                                   /* total:  348 bytes */


/*****************************************************************************/
/*--------------- Prototypes of functions defined in this file --------------*/
const char *nifti_units_string      ( int uu ) ;
const char *nifti_intent_string     ( int ii ) ;
const char *nifti_xform_string      ( int xx ) ;
const char *nifti_slice_string      ( int ss ) ;
const char *nifti_orientation_string( int ii ) ;

int   nifti_is_inttype(int dt);
void  nifti_swap_bytes(size_t n, int siz, void *ar);

void  swap_nifti_header ( struct nifti_1_header *h , int is_nifti ) ;
int   nifti_swap_as_analyze( nifti_analyze75 *h );


void nifti_disp_lib_hist( void ) ;     /* to display library history */
void nifti_disp_lib_version( void ) ;  /* to display library version */
int nifti_disp_matrix_orient( const char * mesg, Matrix4d mat );
int nifti_disp_type_list( int which );

int    disp_nifti_1_header(const char * info, const nifti_1_header * hp ) ;
void   nifti_set_debug_level( int level ) ;
void   nifti_set_skip_blank_ext( int skip ) ;
void   nifti_set_allow_upper_fext( int allow ) ;
int nifti_short_order(void) ;              /* CPU byte order */


/* Orientation codes that might be returned from nifti_mat44_to_orientation().*/

#define NIFTI_L2R  1    /* Left to Right         */
#define NIFTI_R2L  2    /* Right to Left         */
#define NIFTI_P2A  3    /* Posterior to Anterior */
#define NIFTI_A2P  4    /* Anterior to Posterior */
#define NIFTI_I2S  5    /* Inferior to Superior  */
#define NIFTI_S2I  6    /* Superior to Inferior  */

void nifti_matrix4d_to_orientation( Matrix4d R , int *icod, int *jcod, int *kcod ) ;

/*--------------------- Low level IO routines ------------------------------*/
/* other routines */
int    nifti_hdr_looks_good        (const nifti_1_header * hdr);
int    nifti_is_valid_ecode        (int ecode);
/*int    nifti_add_extension(nifti_image * nim, const char * data, int len,
                           int ecode );
int    nifti_compiled_with_zlib    (void);
int    nifti_copy_extensions (nifti_image *nim_dest,const nifti_image *nim_src);
int    nifti_free_extensions (nifti_image *nim);
int  * nifti_get_intlist     (int nvals , const char *str);
char * nifti_strdup          (const char *str);
int    valid_nifti_extensions(const nifti_image *nim);*/


/*-------------------- Some C convenience macros ----------------------------*/

/* NIfTI-1.1 extension codes:
 see http://nifti.nimh.nih.gov/nifti-1/documentation/faq#Q21 */

#define NIFTI_ECODE_IGNORE           0  /* changed from UNKNOWN, 29 June 2005 */

#define NIFTI_ECODE_DICOM            2  /* intended for raw DICOM attributes  */

#define NIFTI_ECODE_AFNI             4  /* Robert W Cox: rwcox@nih.gov
http://afni.nimh.nih.gov/afni      */

#define NIFTI_ECODE_COMMENT          6  /* plain ASCII text only              */

#define NIFTI_ECODE_XCEDE            8  /* David B Keator: dbkeator@uci.edu
http://www.nbirn.net/Resources
/Users/Applications/
/xcede/index.htm              */

#define NIFTI_ECODE_JIMDIMINFO      10  /* Mark A Horsfield:
mah5@leicester.ac.uk
http://someplace/something         */

#define NIFTI_ECODE_WORKFLOW_FWDS   12  /* Kate Fissell: fissell@pitt.edu
http://kraepelin.wpic.pitt.edu
/~fissell/NIFTI_ECODE_WORKFLOW_FWDS
/NIFTI_ECODE_WORKFLOW_FWDS.html   */

#define NIFTI_ECODE_FREESURFER      14  /* http://surfer.nmr.mgh.harvard.edu  */

#define NIFTI_ECODE_PYPICKLE        16  /* embedded Python objects
http://niftilib.sourceforge.net
/pynifti                     */

/* LONI MiND codes: http://www.loni.ucla.edu/twiki/bin/view/Main/MiND */
#define NIFTI_ECODE_MIND_IDENT      18  /* Vishal Patel: vishal.patel@ucla.edu*/
#define NIFTI_ECODE_B_VALUE         20
#define NIFTI_ECODE_SPHERICAL_DIRECTION 22
#define NIFTI_ECODE_DT_COMPONENT    24
#define NIFTI_ECODE_SHC_DEGREEORDER 26  /* end LONI MiND codes                */

#define NIFTI_ECODE_VOXBO           28  /* Dan Kimberg: www.voxbo.org         */

#define NIFTI_ECODE_CARET           30  /* John Harwell: john@brainvis.wustl.edu
http://brainvis.wustl.edu/wiki
/index.php/Caret:Documentation
:CaretNiftiExtension             */

#define NIFTI_MAX_ECODE             30  /******* maximum extension code *******/

/* nifti_type file codes */
#define NIFTI_FTYPE_ANALYZE   0
#define NIFTI_FTYPE_NIFTI1_1  1
#define NIFTI_FTYPE_NIFTI1_2  2
#define NIFTI_FTYPE_ASCII     3
#define NIFTI_MAX_FTYPE       3    /* this should match the maximum code */

#undef  MSB_FIRST
#undef  LSB_FIRST
#undef  REVERSE_ORDER
#define LSB_FIRST 1
#define MSB_FIRST 2
#define REVERSE_ORDER(x) (3-(x))    /* convert MSB_FIRST <--> LSB_FIRST */

#define LNI_MAX_NIA_EXT_LEN 100000  /* consider a longer extension invalid */

/*! NIfTI header class */

enum NIFTIIMAGE_MODES
{
	NIFTI_CLOSED = 0,
	NIFTI_READ = 'r',
	NIFTI_WRITE = 'w'
};

struct NiftiDataTypeInfo
{
	int size, swapsize;
	std::string name;
};
typedef std::map<int, NiftiDataTypeInfo> DTMap;

class NiftiImage
{
	private:
		int _dim[8], _nvox;               /*!< number of voxels = nx*ny*nz*...*nw   */
		float _voxdim[8];
		Affine3d _qform, _sform, _inverse;
		
		std::string _basename, _imgname, _hdrname;     /* Paths to header and image files*/
		int _voxoffset, _byteorder; /* Offset and byte swap info */
		int _datatype;                      /* Datatype on disk */
		int _num_ext;
		nifti1_extension *_ext_list;
		char _mode;
		bool _gz;
		
		FILE *_file;
		gzFile _gzFile;
		
		void setFilenames(const std::string &filename);
		static int needs_swap(short dim0, int hdrsize);
		static float fixFloat(const float f);
		
		size_t read(void *buff, size_t size, size_t nmemb);
		size_t write(const void *buff, size_t size, size_t nmemb);
		long seek(long offset, int whence);
		int rewind();
		
		void readHeader(std::string path);
		void writeHeader(std::string path);
		void *readBuffer(size_t start, size_t length);
		void writeBuffer(void *buffer, size_t start, size_t length);
		template<typename T> T *convertBuffer(const void *raw, const size_t bufferSize)
		{
			T *converted = (T *)malloc(bufferSize * sizeof(T));
			for (int i = 0; i < bufferSize; i++) {
				switch (_datatype) {
					// NOTE: C++11 specifies that C++ 'complex<type>' and C 'type complex'
					// should be interchangeable even at pointer level, so the following
					// should work.
					case NIFTI_TYPE_INT8:      converted[i] = (T)((char *)raw)[i]; break;
					case NIFTI_TYPE_UINT8:     converted[i] = (T)((unsigned char *)raw)[i]; break;
					case NIFTI_TYPE_INT16:     converted[i] = (T)((short *)raw)[i]; break;
					case NIFTI_TYPE_UINT16:    converted[i] = (T)((unsigned short *)raw)[i]; break;
						
					case NIFTI_TYPE_RGB24:     break ;
					case NIFTI_TYPE_RGBA32:    break ;
						
					case NIFTI_TYPE_INT32:     converted[i] = (T)((int *)raw)[i]; break;
					case NIFTI_TYPE_UINT32:    converted[i] = (T)((unsigned int *)raw)[i]; break;
					
					case NIFTI_TYPE_FLOAT32:   converted[i] = (T)((float *)raw)[i]; break;
					//case DT_COMPLEX64: converted[i] = (T)((std::complex<float> *)raw)[i]; break;
						
					case NIFTI_TYPE_FLOAT64:   converted[i] = (T)((double *)raw)[i]; break;
					case NIFTI_TYPE_INT64:     converted[i] = (T)((long *)raw)[i]; break;
					case NIFTI_TYPE_UINT64:    converted[i] = (T)((unsigned long *)raw)[i]; break;
						
					case NIFTI_TYPE_FLOAT128:  converted[i] = (T)((long double *)raw)[i]; break;
						
					//case DT_COMPLEX128: converted[i] = (T)((std::complex<double> *)raw)[i]; break;
						
					//case DT_COMPLEX256: converted[i] = (T)((std::complex<long double> *)raw)[i]; break;
				}
			}
			return converted;
		}
		template<typename T> void *convertToBuffer(const T *raw, const size_t bufferSize)
		{
			void *converted = malloc(DTypes.find(_datatype)->second.size * bufferSize);
			for (int i = 0; i < bufferSize; i++) {
				switch (_datatype) {
					// NOTE: C++11 specifies that C++ 'complex<type>' and C 'type complex'
					// should be interchangeable even at pointer level, so the following
					// should work.
					case NIFTI_TYPE_INT8:      ((char *)converted)[i] =          raw[i]; break;
					case NIFTI_TYPE_UINT8:     ((unsigned char *)converted)[i] =  raw[i]; break;
					case NIFTI_TYPE_INT16:     ((short *)converted)[i] =          raw[i]; break;
					case NIFTI_TYPE_UINT16:    ((unsigned short *)converted)[i] = raw[i]; break;
						
					case NIFTI_TYPE_RGB24:     break ;
					case NIFTI_TYPE_RGBA32:    break ;
						
					case NIFTI_TYPE_INT32:     ((int *)converted)[i] =          raw[i]; break;
					case NIFTI_TYPE_UINT32:    ((unsigned int *)converted)[i] = raw[i]; break;
					
					case NIFTI_TYPE_FLOAT32:   ((float *)converted)[i] =        raw[i]; break;
					//case DT_COMPLEX64: converted[i] = (T)((std::complex<float> *)raw)[i]; break;
						
					case NIFTI_TYPE_FLOAT64:   ((double *)converted)[i] =        raw[i]; break;
					case NIFTI_TYPE_INT64:     ((long *)converted)[i] =          raw[i]; break;
					case NIFTI_TYPE_UINT64:    ((unsigned long *)converted)[i] = raw[i]; break;
						
					case NIFTI_TYPE_FLOAT128:  ((long double *)converted)[i] =   raw[i]; break;
						
					//case DT_COMPLEX128: converted[i] = (T)((std::complex<double> *)raw)[i]; break;
						
					//case DT_COMPLEX256: converted[i] = (T)((std::complex<long double> *)raw)[i]; break;
				}
			}
			return converted;
		}
		static const DTMap DTypes;
		static const std::string &DTypeToString(const int dtype);
		static const bool DTypeIsValid(const int dtype);
	public:
		~NiftiImage();
		NiftiImage();
		NiftiImage(const NiftiImage &clone);
		NiftiImage(const int nx, const int ny, const int nz, const int nt,
		           const float dx, const float dy, const float dz, const float dt,
				   const int datatype);
		NiftiImage(const std::string filename);
		NiftiImage &operator=(const NiftiImage &other);
		static void printDTypeList();
		
		bool open(std::string filename, char mode);
		void close();
		const std::string &basename();
		void *readRawVolume(const int vol);
		void *readRawAllVolumes();
		template<typename T> T *readVolume(const int vol)
		{
			size_t bytesPerVolume = voxelsPerVolume() *
			                        DTypes.find(_datatype)->second.size;
			void *raw = readBuffer(vol * bytesPerVolume, bytesPerVolume);
			T *converted = convertBuffer<T>(raw, voxelsPerVolume());
			free(raw);
			return converted;
		}
		
		template<typename T> T *readAllVolumes()
		{
			void *raw =	readBuffer(0, nvox() * DTypes.find(_datatype)->second.size);
			T *converted = convertBuffer<T>(raw, nvox());
			free(raw);
			return converted;
		}
		
		template<typename T> void writeVolume(const int vol, T *data)
		{
			size_t bytesPerVolume = voxelsPerVolume() * DTypes.find(_datatype)->second.size;
			void *converted = convertToBuffer<T>(data, voxelsPerVolume());
			writeBuffer(converted, vol * bytesPerVolume, bytesPerVolume);
			free(converted);
		}
		
		template<typename T> void writeAllVolumes(T *data)
		{
			void *converted = convertBuffer<T>(data, nvox());
			writeBuffer(converted, 0, nvox() * DTypes.find(_datatype)->second.size);
			free(converted);
		}
		
		int ndim() const;
		int nx() const;
		int ny() const;
		int nz() const;
		int nt() const;
		int voxelsPerVolume() const;
		int nvox() const;
		void setnt(const int nt);
		void setDims(int nx, int ny, int nz, int nt);
		
		int datatype() const;
		void setDatatype(const int dt);
		bool volumesCompatible(const NiftiImage &other) const; /// Check whether a volume from this header can be used in calculations with another header.
		
		float dx() const;
		float dy() const;
		float dz() const;
		float dt() const;
		void setVoxDims(double dx, double dy, double dz, double dt);
		
		float scaling_slope;
		float scaling_inter;
		float calibration_min;
		float calibration_max;
		
		int qform_code;
		int sform_code;
		const Matrix4d &qform() const;
		const Matrix4d &sform() const;
		void set_qform(const Matrix4d& new_qform);
		void set_qform(const Affine3d &R, const Affine3d &T);
		void set_qform(const double b, const double c, const double d,
					   const double x, const double y, const double z);
		void set_sform(const Matrix4d& new_sform);
		const Matrix4d &ijk_to_xyz() const;
		const Matrix4d &xyz_to_ijk() const;
		
		int freq_dim ;               /*!< indexes (1,2,3, or 0) for MRI    */
		int phase_dim;               /*!< directions in dim[]/pixdim[]     */
		int slice_dim;               /*!< directions in dim[]/pixdim[]     */
		
		int   slice_code;             /*!< code for slice timing pattern    */
		int   slice_start;            /*!< index for start of slices        */
		int   slice_end;              /*!< index for end of slices          */
		float slice_duration;         /*!< time between individual slices   */
		float toffset;                /*!< time coordinate offset */
	
		int xyz_units;                /*!< dx,dy,dz units: NIFTI_UNITS_* code  */
		int time_units;               /*!< dt       units: NIFTI_UNITS_* code  */
		int   intent_code ;           /*!< statistic type (or something)       */
		float intent_p1 ;             /*!< intent parameters                   */
		float intent_p2 ;             /*!< intent parameters                   */
		float intent_p3 ;             /*!< intent parameters                   */
		std::string intent_name;      /*!< optional description of intent data */
	
		std::string description;      /*!< optional text to describe dataset   */
		std::string aux_file;         /*!< auxiliary filename                  */
		
		int                num_ext ;  /*!< number of extensions in ext_list       */
		nifti1_extension * ext_list ; /*!< array of extension structs (with data) */
		analyze_75_orient_code analyze75_orient; /*!< for old analyze files, orient */
};
#endif /* _NIFTI_IO_HEADER_ */

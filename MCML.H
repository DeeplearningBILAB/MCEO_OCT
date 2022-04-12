/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2016.
*
*	Monte Carlo simulation of photon propogation in 
**	multi-layered turbid media with embedded object
**	for Class I and Class II signal of optical 
**	coherence tomography in ANSI Standard C.
****
*	Starting Date:		07/2015.
*	Completion Date:	02/2016.
*
**	Vijitha Periyasamy, B.E.
**	Manojit Pramanik, Ph.D.
**	Biomedical Imaging Laboratory,
**  School of Chemical and Biomedical Engineering,
**	Nanyang Technological University, Singapore 637459.
*
*	This program was based on:
*	(1) The Pascal code written by Marleen Keijzer and 
*	Steven L. Jacques in this laboratory in 1989, which
*	deals with multi-layered turbid media.
*
*	(2) Algorithm for semi-infinite turbid medium by 
*	S.A. Prahl, M. Keijzer, S.L. Jacques, A.J. Welch, 
*	SPIE Institute Series Vol. IS 5 (1989), and by 
*	A.N. Witt, The Astrophysical journal Supplement
*	Series 35, 1-6 (1977).
*	
*	(3)	MCML-Monte Carlo modeling of photon transport  
*   in multi-layered tissues, L.-H. Wang, S. L. Jacques,
*	and L.-Q. Zheng Computer Methods and programs in 
*	Biomedicine, 47, 131-146 (1995)
*
*   (4) Monte Carlo simulation of light transport in turbid 
*   medium with embedded object - spherical, cylindrical, 
*   ellipsoidal, or cuboidal object embedded within 
*   multilayered tissues, V. Periyasamy and M. Pramanik, 
*   Journal of Biomedical Optics 19(4), 045003 (2014).
*
*   (5) Importance sampling based Monte Carlo simulation of time
*	domain optical coherence tomography with embedded objects, 
*	V. Periyasamy, and M. Pramanik, Applied Optics, (2016).
*
*   (6) Improved importance sampling for Monte Carlo simulation
*	of time-domain optical coherence tomography, I. T. Lima, 
*	A. Kalra, and S. S. Sherif, Biomedical Optics Express, 2(5),
*	1069 (2011).
**	

*	Major modifications include:
*		. Conform to ANSI Standard C.
**		. Object embedded completely in a layer, 
**        with refractive index different from that of 
**        surrounding tissue
**		and
**		Improved importance sampling for class I and
**		class II photons of OCT.
****
*	General Naming Conventions:
*	Preprocessor names: all capital letters, 
*		e.g. #define PREPROCESSORS
*	Globals: first letter of each word is capital, no 
*		underscores, 
*		e.g. short GlobalVar;
*	Dummy variables:  first letter of each word is capital,
*		and words are connected by underscores, 
*		e.g. void NiceFunction(char Dummy_Var);
*	Local variables:  all lower cases, words are connected 
*		by underscores,
*		e.g. short local_var;
*	Function names or data types:  same as Globals.
*
****
*	Dimension of length: cm.
*
****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define PI 3.1415926
#define WEIGHT 1E-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */
#define STRLEN 256		/* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

/****************** Stuctures *****************************/

/****
*	Structure used to describe a photon packet.
****/
typedef struct {
	double x, y ,z;	/* Cartesian coordinates.[cm] */
	double ux, uy, uz;/* directional cosines of a photon. */
	double w;			/* weight. */
	Boolean dead;		/* 1 if photon is terminated. */
	short layer;		/* index to layer where the photon */
	/* packet resides. */
	double s;			/* current step size. [cm]. */
	double sleft;		/* step size left. dimensionless [-]. */
	Boolean inObj;    /** Is photon within object (1) or out of it (0) **/
	double zMax;		/** Maximum z reached by the photon **/
	double biasfact;	/** Keep track of the biasing factor **/
	double opLen;		/** Optical path lenght **/
	long int first_bias_scatter; /** track the first biased scattering **/
	double launchLoc;	/** Store the x location of launch during B-scan **/
} PhotonStruct;

/****
*	Structure used to describe the geometry and optical
*	properties of a layer.
*	z0 and z1 are the z coordinates for the upper boundary
*	and lower boundary respectively.
*
*	cos_crit0 and cos_crit1 are the cosines of the 
*	critical angle of total internal reflection for the
*	upper boundary and lower boundary respectively.
*	They are set to zero if no total internal reflection
*	exists.
*	They are used for computation speed.
****/
typedef struct {
	double z0, z1;	/* z coordinates of a layer. [cm] */
	double n;			/* refractive index of a layer. */
	double mua;	    /* absorption coefficient. [1/cm] */
	double mus;	    /* scattering coefficient. [1/cm] */
	double g;		    /* anisotropy. */
	double cos_crit0,	cos_crit1;	
} LayerStruct;

/****
**	Structure used to describe the geometry and optical
**	properties of a object.
**	cx, cy, cz are the x, y, and the z coordinates of the
**  center of the object, respectively.
**  crx, cry and crz are the radius of the object in 
**  in x, y, and z coordinates in cm.
**  Min and Max in x, y, and z coordinates are to define
**  the planes of cuboid.
*
**	cos_crit0 and cos_crit1 are the cosines of the 
**	critical angle of total internal reflection for the
**	surrounding medium to the object and object to the 
**  surrounding medium, respectively.
**	They are set to zero if no total internal reflection
**	exists.
**	Included by Vijitha et al
****/
typedef struct {
	double cx, cy, cz;		/** x,y,z coordinates of center of object. [cm] **/
	double crx, cry, crz;		 /** dimensions of object **/
	double XMin, YMin, ZMin;		/** Minmum cordinates of cuboid **/
	double XMax, YMax, ZMax;		/** Maxmum cordinates of cuboid **/
	double n;			 /** refractive index of a layer. **/
	double mua;	     /** absorption coefficient. [1/cm] **/
	double mus;	     /** scattering coefficient. [1/cm] **/
	double g;		     /** anisotropy. **/
	double cos_crit0,	cos_crit1;	
} ObjStruct;

/****
*	Input parameters for each independent run.
*
*	z and r are for the cylindrical coordinate system. [cm]
*	a is for the angle alpha between the photon exiting 
*	direction and the surface normal. [radian]
*
*	The grid line separations in z, r, and alpha
*	directions are dz, dr, and da respectively.  The numbers 
*	of grid lines in z, r, and alpha directions are
*	nz, nr, and na respectively.
*
*	The member layerspecs will point to an array of 
*	structures which store parameters of each layer. 
*	This array has (number_layers + 2) elements. One
*	element is for a layer.
*	The layers 0 and (num_layers + 1) are for top ambient 
*	medium and the bottom ambient medium respectively.
*
**	Layer in which object is embedded and its optical properties
**	are read from the template file.
****/
typedef struct {
	char	 out_fname[STRLEN];	/* output file name. */
	char	 out_fformat;		/* output file format. */
	/* 'A' for ASCII, */
	/* 'B' for binary. */
	double	 num_photons; 		/* to be traced. */
	double Wth; 				/* play roulette if photon */
	/* weight < Wth.*/

	double dz;				/* z grid separation.[cm] */ 
	double dr;				/* r grid separation.[cm] */
	double da;				/* alpha grid separation. */
	/* [radian] */
	short nz;					/* array range 0..nz-1. */
	short nr;					/* array range 0..nr-1. */
	short na;					/* array range 0..na-1. */

	short	num_layers;			/* number of layers. */
	short objLayer;				/** 0 if the object doesn't exist else layer in which the object is **/
	short objCode;				/**(0) No embedded object (1) Sphere (2) Cylinder (3) Ellipsoid (4) Cuboid **/
	LayerStruct * layerspecs;	/* layer parameters. */	
	ObjStruct * ObjSpecs;  		/** object parameters. **/

	double latRes;				/** lateral step size for B-scan in cm **/
	short num_scan;				/** Number of lateral scans **/ 
	double x_start;				/** x-start for B-scan **/
	double x_end;				/** x-end for B-scan **/
	int launchCnt;				/** Launch counter **/
	double lc;					/** coherence length **/
	double theta_max;			/** detector receive angle in degree **/
	double r_max;				/** collecting radius **/


	double bias_p;				/** probability of biasing towards detector **/
	double bias_a;				/** biasing coefficient **/

} InputStruct;

/****
*	Structures for scoring physical quantities. 
*	z and r represent z and r coordinates of the 
*	cylindrical coordinate system. [cm]
*	a is the angle alpha between the photon exiting 
*	direction and the normal to the surfaces. [radian]
*	See comments of the InputStruct.
*	See manual for the physcial quantities.
**	Also included is the absorbance within object.
****/
typedef struct {
	double    Rsp;	/* specular reflectance. [-] */
	double ** Rd_ra;	/* 2D distribution of diffuse */
	/* reflectance. [1/(cm2 sr)] */
	double *  Rd_r;	/* 1D radial distribution of diffuse */
	/* reflectance. [1/cm2] */
	double *  Rd_a;	/* 1D angular distribution of diffuse */
	/* reflectance. [1/sr] */
	double    Rd;		/* total diffuse reflectance. [-] */

	double ** A_rz;	/* 2D probability density in turbid */
	/* media over r & z. [1/cm3] */

	double *  A_z;	/* 1D probability density over z. */
	/* [1/cm] */
	double *  A_l;	/* each layer's absorption */
	/* probability. [-] */
	double    A;		/* total absorption probability. [-] */

	double    A_obj;		/** absorption probability in the object. [-] **/

	double ** Tt_ra;		/* 2D distribution of total */
	/* transmittance. [1/(cm2 sr)] */
	double *  Tt_r;	/* 1D radial distribution of */
	/* transmittance. [1/cm2] */
	double *  Tt_a;	/* 1D angular distribution of */
	/* transmittance. [1/sr] */
	double    Tt;		/* total transmittance. [-] */

	double *  Rd_z1;	/** Depth reached by type 1 diffussed reflected photons **/
	double *  Rd_z2;	/** Depth reached by type 2 diffussed reflected photons **/

	double ** Rd_z1_bs;	/** Class 1 photons for B-scan  **/
	double ** Rd_z2_bs;	/** Class 2 photons for B-scan  **/
} OutStruct;

/***********************************************************
*	Routine prototypes for dynamic memory allocation and 
*	release of arrays and matrices.
*	Modified from Numerical Recipes in C.
****/
double  *AllocVector(short, short);
double  **AllocMatrix(short, short,short, short);
void 	FreeVector(double *, short, short);
void 	FreeMatrix(double **, short, short, short, short);
void 	nrerror(char *);
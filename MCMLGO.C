/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2016.
*
*	Launch, move, and record photon weight.
*
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al.
****/

#include "mcml.h"

#define STANDARDTEST 0
/* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0     
/* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)	
/* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6		
/* cosine of about 1.57 - 1e-6 rad. */


/***********************************************************
*	A random number generator from Numerical Recipes in C.
****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

float ran3(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
			inext=0;
			inextp=31;
			*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return (double)mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
*	Generate a random number between 0 and 1.  Take a 
*	number as seed the first time entering the function.  
*	The seed is limited to 1<<15.  
*	We found that when idum is too large, ran3 may return 
*	numbers beyond 0 and 1.
****/
double RandomNum(void)
{
	static Boolean first_time=1;
	static int idum;	/* seed for ran3. */

	if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
		idum = - 1;
#else
		idum = -(int)time(NULL)%(1<<15);
		/* use 16-bit integer as the seed. */
#endif
		ran3(&idum);
		first_time = 0;
		idum = 1;
	}

	return( (double)ran3(&idum) );
}

/***********************************************************
*	Compute the specular reflection. 
*
*	If the first layer is a turbid medium, use the Fresnel
*	reflection from the boundary of the first layer as the 
*	specular reflectance.
*
*	If the first layer is glass, multiple reflections in
*	the first layer is considered to get the specular
*	reflectance.
*
*	The subroutine assumes the Layerspecs array is correctly 
*	initialized.
****/
double Rspecular(LayerStruct * Layerspecs_Ptr)
{
	double r1, r2;
	/* direct reflections from the 1st and 2nd layers. */
	double temp;

	temp =(Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
		/(Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
	r1 = temp*temp;

	if((Layerspecs_Ptr[1].mua == 0.0) 
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n)
				/(Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);
			r2 = temp*temp;
			r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
	}

	return (r1);	
}

/***********************************************************
*	Compute the Fresnel reflectance.
*
*	Make sure that the cosine of the incident angle a1
*	is positive, and the case when the angle is greater 
*	than the critical angle is ruled out.
*
* 	Avoid trigonometric function operations as much as
*	possible, because they are computation-intensive.
****/
double RFresnel(double n1,	/* incident refractive index.*/
	double n2,	/* transmit refractive index.*/
	double ca1,	/* cosine of the incident */
	/* angle. 0<a1<90 degrees. */
	double * ca2_Ptr)  /* pointer to the */
	/* cosine of the transmission */
	/* angle. a2>0. */
{
	double r;

	if(n1==n2) {			  	/** matched boundary. **/
		*ca2_Ptr = ca1;
		r = 0.0;
	}
	else if(ca1>COSZERO) {	/** normal incident. **/
		*ca2_Ptr = ca1;
		r = (n2-n1)/(n2+n1);
		r *= r;
	}
	else if(ca1<COS90D)  {	/** very slant. **/
		*ca2_Ptr = 0.0;
		r = 1.0;
	}
	else  {			  		/** general. **/
		double sa1, sa2;	
		/* sine of the incident and transmission angles. */
		double ca2;

		sa1 = sqrt(1-ca1*ca1);
		sa2 = n1*sa1/n2;
		if(sa2>=1.0) {	
			/* double check for total internal reflection. */
			*ca2_Ptr = 0.0;
			r = 1.0;
		}
		else  {
			double cap, cam;	/* cosines of the sum ap or */
			/* difference am of the two */
			/* angles. ap = a1+a2 */
			/* am = a1 - a2. */
			double sap, sam;	/* sines. */

			*ca2_Ptr = ca2 = sqrt(1-sa2*sa2);

			cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
			cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
			sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
			sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
			r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
			/* rearranged for speed. */
		}
	}
	return(r);
}

/***********************************************************
*	Initialize a photon packet.
****/
void LaunchPhoton(double Rspecular,
	LayerStruct  * Layerspecs_Ptr,
	PhotonStruct * Photon_Ptr)
{
	Photon_Ptr->w	 	= 1.0 - Rspecular;	
	Photon_Ptr->dead 	= 0;
	Photon_Ptr->layer = 1;
	Photon_Ptr->s	= 0;
	Photon_Ptr->sleft= 0;

	Photon_Ptr->x 	= 0.0;	
	Photon_Ptr->y	= 0.0;	
	Photon_Ptr->z	= 0.0;	
	Photon_Ptr->ux	= 0.0;	
	Photon_Ptr->uy	= 0.0;	
	Photon_Ptr->uz	= 1;
	Photon_Ptr->inObj	= 0;

	Photon_Ptr->biasfact    = 1.0; /**For OCT---- **/
	Photon_Ptr->opLen       = 0.0;	/**For OCT **/
	Photon_Ptr->zMax 	      = 0.0;	/**For OCT **/
	Photon_Ptr->first_bias_scatter = 0;	/**For OCT **/

	if((Layerspecs_Ptr[1].mua == 0.0) 
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			Photon_Ptr->layer 	= 2;
			Photon_Ptr->z	= Layerspecs_Ptr[2].z0;	
	}
}


/***********************************************************
**	Initialize a photon packet when there is a object.
****/
void LaunchPhotonObj(double Rspecular,
					 InputStruct * In_Ptr,
					 PhotonStruct * Photon_Ptr,
					 long i_photon)
{
	FILE *file;
	double uz1, ni, nt, uz, r;
	short launchConfig;
	double numOfPhotons;
	ni = In_Ptr->layerspecs[0].n;
	nt = In_Ptr->layerspecs[1].n;
	numOfPhotons = In_Ptr->num_photons;

	Photon_Ptr->inObj   = 0;  /** Initialize that the photon is outside object. **/
	Photon_Ptr->w	 	= 1.0 - Rspecular;	
	Photon_Ptr->dead 	= 0;
	Photon_Ptr->layer   = 1;
	Photon_Ptr->s	    = 0;
	Photon_Ptr->sleft   = 0;

	Photon_Ptr->x 	= In_Ptr->x_start+In_Ptr->latRes*In_Ptr->launchCnt;	/** For B-scan **/	
	Photon_Ptr->launchLoc = Photon_Ptr->x;	/**Location of photon in x-axis for B-scan is saved **/
	Photon_Ptr->y	= 0.0;	
	Photon_Ptr->z	= 0.0;	
	Photon_Ptr->ux	= 0.0;	
	Photon_Ptr->uy	= 0.0;	
	Photon_Ptr->uz	= 1;
	Photon_Ptr->biasfact    = 1.0; /**For OCT **/
	Photon_Ptr->opLen       = 0.0;	/**For OCT **/
	Photon_Ptr->zMax 	      = 0.0;	/**For OCT **/
	Photon_Ptr->first_bias_scatter = 0;	/**For OCT **/

	if((In_Ptr->layerspecs[1].mua == 0.0) 
		&& (In_Ptr->layerspecs[1].mus == 0.0))  { /* glass layer. */
			Photon_Ptr->layer = 2;
			Photon_Ptr->z = In_Ptr->layerspecs[2].z0;	
	}
}


/***********************************************************
*	Choose (sample) a new theta angle for photon propagation
*	according to the anisotropy.
*
*	If anisotropy g is 0, then
*		cos(theta) = 2*rand-1.
*	otherwise
*		sample according to the Henyey-Greenstein function.
*
*	Returns the cosine of the polar deflection angle theta.
****/
double SpinTheta(double g, double randNum)
{
	double cost;

	if(g == 0.0) 
		cost = 2*RandomNum() -1;
	else {
		double temp = (1-g*g)/(1-g+2*g*randNum);
		cost = (1+g*g - temp*temp)/(2*g);
		if(cost < -1) cost = -1;
		else if(cost > 1) cost = 1;
	}
	return(cost);
}


/***********************************************************
*	Choose a new direction for photon propagation by 
*	sampling the polar deflection angle theta and the 
*	azimuthal angle psi.
*
*	Note:
*  	theta: 0 - pi so sin(theta) is always positive 
*  	feel free to use sqrt() for cos(theta).
* 
*  	psi:   0 - 2pi 
*  	for 0-pi  sin(psi) is + 
*  	for pi-2pi sin(psi) is - 
****/
void Spin(InputStruct  *	In_Ptr,
	PhotonStruct * Photon_Ptr)
{
	double cost, sint;	/* cosine and sine of the */
	/* polar deflection angle theta. */
	double cosp, sinp;	/* cosine and sine of the */
	/* azimuthal angle psi. */
	double ux = Photon_Ptr->ux;
	double uy = Photon_Ptr->uy;
	double uz = Photon_Ptr->uz;
	double psi;
	double g=0.9;
	cost = SpinTheta(g,RandomNum());
	sint = sqrt(1.0 - cost*cost);	
	/* sqrt() is faster than sin(). */

	psi = 2.0*PI*RandomNum(); /* spin psi 0-2pi. */
	cosp = cos(psi);
	if(psi<PI)
		sinp = sqrt(1.0 - cosp*cosp);	
	/* sqrt() is faster than sin(). */
	else
		sinp = - sqrt(1.0 - cosp*cosp);	

	if(fabs(uz) > COSZERO)  { 	/* normal incident. */
		Photon_Ptr->ux = sint*cosp;
		Photon_Ptr->uy = sint*sinp;
		Photon_Ptr->uz = cost*SIGN(uz);	
		/* SIGN() is faster than division. */
	}
	else  {		/* regular incident. */
		double temp = sqrt(1.0 - uz*uz);
		Photon_Ptr->ux = sint*(ux*uz*cosp - uy*sinp)
			/temp + ux*cost;
		Photon_Ptr->uy = sint*(uy*uz*cosp + ux*sinp)
			/temp + uy*cost;
		Photon_Ptr->uz = -sint*cosp*temp + uz*cost;
	}
}
/***********************************************************
**	theta computation for biased scattering
****/
double SpinTheta_forA(double a,double rand_num)
{
	double cost_a;
	double temp = rand_num*((1/(1-a))-(1/sqrt((a*a)+1))) + (1/sqrt((a*a)+1));	/**Biomedical optics express 3, 692 (2012): equation 2**/
	cost_a = (1/(2*a))*((a*a)+1-(1/(temp*temp)));
	return(cost_a);
}
/***********************************************************
*	Biased scattering for OCT 
****/
void SpinBiased(InputStruct  *  In_Ptr,
	PhotonStruct * Photon_Ptr)
{
	double cost, sint;      /* cosine and sine of the, polar deflection angle theta. */
	double cosp, sinp;      /* cosine and sine of the, azimuthal angle psi. */

	double ux, uy, uz, ux1, uy1, uz1, uxf, uyf, uzf;
	double psi;
	double alpha1, beta1, gamma1, mag, w1;
	double g, a, p, costS;
	double biasVal;
	double cost_a, cost_g, cost_g_bias;
	double rand_num_a, rand_num_p;
	double f_cost_g, f_cost_a, temp1;

	ux1 = Photon_Ptr->ux;	/** Direction cosies of the photon before scattering **/
	uy1 = Photon_Ptr->uy;
	uz1 = Photon_Ptr->uz;

	g = In_Ptr->layerspecs[Photon_Ptr->layer].g;
	a = In_Ptr->bias_a;
	p = In_Ptr->bias_p;

	w1 = 0; /**Getting the apparent location of detector in z **/
	mag = sqrt(Photon_Ptr->x*Photon_Ptr->x + Photon_Ptr->y*Photon_Ptr->y + (w1-Photon_Ptr->z)*(w1-Photon_Ptr->z));
	if(mag < 1E-12) {
		alpha1 = 0;
		beta1 = 0;
		gamma1 = -1;
	} else {
		alpha1 = -Photon_Ptr->x/mag;	/** getting the direction from the photon current location towards the detector **/
		beta1 = -Photon_Ptr->y/mag;
		gamma1 = (w1-Photon_Ptr->z)/mag;
	}

	uxf = alpha1;		/** Direction cosines between the photon and fiber before scattering **/
	uyf = beta1;
	uzf = gamma1;

	rand_num_a = RandomNum();
	cost_g = SpinTheta(In_Ptr->layerspecs[Photon_Ptr->layer].g,rand_num_a);	/** HG function	**/
	cost_a = SpinTheta_forA(a,rand_num_a,cost_g);	/** IS function; Biomedical optics express 3, 692 (2012): Equation 9 **/

	rand_num_p = RandomNum();
	if (rand_num_p<p) {	/** Biased scattering **/
		cost = cost_a;
		ux = alpha1;
		uy = beta1;
		uz = gamma1;

		if(Photon_Ptr->first_bias_scatter == 0) 	Photon_Ptr->first_bias_scatter = 1;
		
	} else {	/** Unbiased scattering	**/
		cost = cost_g;
		ux = Photon_Ptr->ux;
		uy = Photon_Ptr->uy;
		uz = Photon_Ptr->uz;
	}
	sint = sqrt(1.0 - cost*cost);	/* sqrt() is faster than sin(). */

	psi = 2.0 * PI * RandomNum(); /* spin psi 0-2pi. */
	cosp = cos(psi);

	if(psi<PI) sinp = sqrt(1.0 - cosp*cosp);  /* sqrt() is faster than sin(). */
	else 	sinp = - sqrt(1.0 - cosp*cosp);

	if(fabs(uz) > COSZERO) {        /* normal incident. */
		Photon_Ptr->ux = sint*cosp;
		Photon_Ptr->uy = sint*sinp;
		Photon_Ptr->uz = cost*SIGN(uz);  /* SIGN() is faster than division. */
	} else {                /* regular incident. */
		double temp = sqrt(1.0 - uz*uz);
		Photon_Ptr->ux = sint*(ux*uz*cosp - uy*sinp)/temp + ux*cost;
		Photon_Ptr->uy = sint*(uy*uz*cosp + ux*sinp)/temp + uy*cost;
		Photon_Ptr->uz = -sint*cosp*temp + uz*cost;
	}

	cost_a = Photon_Ptr->ux*uxf + Photon_Ptr->uy*uyf +Photon_Ptr->uz*uzf;		/** recompute angle between scattered direction and fiber	**/
	temp1 = 1 + a*a - 2*a*cost_a;
	f_cost_a = (1/(1-((1-a)/sqrt(a*a+1))))*((a*(1-a))/sqrt(temp1*temp1*temp1));	/** Biomedical optics express 3, 692 (2012): equation 2	**/
	
	cost_g = Photon_Ptr->ux*ux1 + Photon_Ptr->uy*uy1 +Photon_Ptr->uz*uz1;	/** recompute angle between scattered direction and direction before scattering	**/
	temp1 = 1 + g*g - 2*g*cost_g;
	f_cost_g = (1-g*g)/(2*sqrt(temp1*temp1*temp1));
	
	biasVal = f_cost_g/(p*f_cost_a + (1-p)*f_cost_g);	/** Biomedical optics express 3, 692 (2012): equation 4 **/

	Photon_Ptr->biasfact = Photon_Ptr->biasfact*biasVal;
}



/***********************************************************
*	Move the photon s away in the current layer of medium.  
****/
void Hop(PhotonStruct *	Photon_Ptr)
{
	double s = Photon_Ptr->s;

	Photon_Ptr->x += s*Photon_Ptr->ux;
	Photon_Ptr->y += s*Photon_Ptr->uy;
	Photon_Ptr->z += s*Photon_Ptr->uz;
	Photon_Ptr->opLen += s;		/**Increment the optical path length of the photon**/
	if (Photon_Ptr->zMax<Photon_Ptr->z) {	/**Check if the photon's current z location is maximum depth reached. **/
		Photon_Ptr->zMax=Photon_Ptr->z;	/**Save the current location as maximum depth when the current location is maximum**/
	}
}			

/***********************************************************
*	If uz != 0, return the photon step size in glass, 
*	Otherwise, return 0.
*
*	The step size is the distance between the current 
*	position and the boundary in the photon direction.
*
*	Make sure uz !=0 before calling this function.
****/
void StepSizeInGlass(PhotonStruct *  Photon_Ptr,
	InputStruct  *  In_Ptr)
{
	double dl_b;	/* step size to boundary. */
	short  layer = Photon_Ptr->layer;
	double uz = Photon_Ptr->uz;

	/* Stepsize to the boundary. */	
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z)
		/uz;
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z)
		/uz;
	else
		dl_b = 0.0;

	Photon_Ptr->s = dl_b;
}

/***********************************************************
*	Pick a step size for a photon packet when it is in 
*	tissue.
*	If the member sleft is zero, make a new step size 
*	with: -log(rnd)/(mua+mus).
*	Otherwise, pick up the leftover in sleft.
*
*	Layer is the index to layer.
*	In_Ptr is the input parameters.
****/
void StepSizeInTissue(PhotonStruct * Photon_Ptr,
	InputStruct  * In_Ptr)
{
	short  layer = Photon_Ptr->layer;
	double mua;
	double mus;
	mua = In_Ptr->layerspecs[layer].mua;
	mus = In_Ptr->layerspecs[layer].mus;

	if(Photon_Ptr->sleft == 0.0) {  /* make a new step. */
		double rnd;

		do rnd = RandomNum(); 
		while( rnd <= 0.0 );    /* avoid zero. */
		Photon_Ptr->s = -log(rnd)/(mua+mus);
	}
	else {	/* take the leftover. */
		Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
		Photon_Ptr->sleft = 0.0;
	}
}

/***********************************************************
**	Pick a step size for a photon packet when it is in object.
**	If the member sleft is zero, make a new step size 
**	with: -log(rnd)/(mua+mus).
**	Otherwise, pick up the leftover in sleft.
*
**	In_Ptr is the input parameters.
****/
void StepSizeInObject(PhotonStruct * Photon_Ptr,
	InputStruct  * In_Ptr)
{
	double mua;
	double mus;
	mua = In_Ptr->ObjSpecs[0].mua;
	mus = In_Ptr->ObjSpecs[0].mus;		

	if(Photon_Ptr->sleft == 0.0) {  /** make a new step. **/
		double rnd;

		do rnd = RandomNum(); 
		while( rnd <= 0.0 );    /** avoid zero. **/
		Photon_Ptr->s = -log(rnd)/(mua+mus);
	}
	else {	/** take the leftover. **/
		Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
		Photon_Ptr->sleft = 0.0;
	}
}


/***********************************************************
*	Check if the step will hit the boundary.
*	Return 1 if hit boundary.
*	Return 0 otherwise.
*
* 	If the projected step hits the boundary, the members
*	s and sleft of Photon_Ptr are updated.
****/
Boolean HitBoundary(PhotonStruct *  Photon_Ptr,
					InputStruct  *  In_Ptr)
{
	double dl_b;  /* length to boundary. */
	short  layer = Photon_Ptr->layer;
	double uz = Photon_Ptr->uz;
	Boolean hit;

	/* Distance to the boundary. */
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1 
				- Photon_Ptr->z)/uz;	/* dl_b>0. */
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0 
				- Photon_Ptr->z)/uz;	/* dl_b>0. */

	if(uz != 0.0 && Photon_Ptr->s > dl_b) {
		/* not horizontal & crossing. */
		double mut = In_Ptr->layerspecs[layer].mua 
					 + In_Ptr->layerspecs[layer].mus;

		Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
		Photon_Ptr->s     = dl_b;
		hit = 1;
	}
	else
		hit = 0;

	return(hit);
}

/************************************
** Check if the photon hits the sphere boundary
*******/

Boolean HitSph(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, cr;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double mut;
	double err =  1E-9;		/** Very small error **/

	double s = Photon_Ptr->s;
	double newx,newy,newz;

	cx = In_Ptr->ObjSpecs[0].cx; 
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm **/

	A = ux*ux + uy*uy + uz*uz;
	B = 2*(ux*(rx-cx) + uy*(ry-cy) + uz*(rz - cz));
	C = (rx-cx)*(rx-cx) + (ry-cy)*(ry-cy) + (rz - cz)*(rz - cz) - cr*cr;

	del = B*B - 4*A*C;

	if (del > 0) { /** If delta is less than 0 then photon path is not intersecting the sphere **/ 

		srDel = sqrt(del);
		rot2 = (-B - srDel)/(2*A);

		if (rot2 >= err) {
			dl_b = rot2;
			if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
				if (Photon_Ptr->inObj == 1)
					mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
				else
					mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

			Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
			Photon_Ptr->s = dl_b;
			return(1);
			}
		} else {
			rot1 = (-B + srDel)/(2*A);
			if (rot1 >= err) {
				dl_b = rot1;
				if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
					if (Photon_Ptr->inObj == 1)
						mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
					else
						mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

				Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
				Photon_Ptr->s = dl_b;
				return(1);
				}
			}
		}

		newx = Photon_Ptr->x + s*Photon_Ptr->ux;
        newy = Photon_Ptr->y + s*Photon_Ptr->uy;
        newz = Photon_Ptr->z + s*Photon_Ptr->uz;

        if(Photon_Ptr->inObj == 0) {//photon currently outside
			if ((newx-cx)*(newx-cx)+(newy-cy)*(newy-cy)+(newz-cz)*(newz-cz) < cr*cr) { //future position inside
				Photon_Ptr->inObj = 1;
            }
        } else {//photon currently inside
			if ((newx-cx)*(newx-cx)+(newy-cy)*(newy-cy)+(newz-cz)*(newz-cz) > cr*cr) {//future position outside
				Photon_Ptr->inObj = 0;
            }
        }
		return(0);
	}

	return(0);	
}


/************************************
** Check if the photon hits the cylinder boundary
*******/

Boolean HitCyl(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, cr, cux, cuy, cuz;
	double d, e_x, e_y, e_z, f, g_x, g_y, g_z, h_x, h_y, h_z, delP_x, delP_y, delP_z;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double X1, Y1, Z1, X2, Y2, Z2;
	double dB1, dB2;
	double alpha1 = 0, alpha2 = 0.0, beta1 = 0, beta2 = 0.0, gamma1 = 0.0, gamma2 = 0.0;
	double u1, v1, w1, mag1, u2, v2, w2, mag2;
	double mut;
	double err =  1E-6;		/** Very small error **/
	double large = 1E24;	/** Large number**/
	
	double s = Photon_Ptr->s;
	double newx, newy, newz;

	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	cux = 1;
	cuy = 0;
	cuz = 0;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference  http://www.mrl.nyu.edu/~dzorin/rendering/lectures/lecture3/lecture3-6pp.pdf **/
	d = ux*cux + uy*cuy + uz*cuz;
	e_x = ux-d*cux;
	e_y = uy-d*cuy;
	e_z = uz-d*cuz;

	A = e_x*e_x + e_y*e_y + e_z*e_z ;

	delP_x = rx-cx;
	delP_y = ry-cy;
	delP_z = rz-cz;

	f = delP_x*cux + delP_y*cuy + delP_z*cuz;
	g_x = delP_x - f*cux;
	g_y = delP_y - f*cuy;
	g_z = delP_z - f*cuz;

	B = 2*(e_x*g_x + e_y*g_y + e_z*g_z);

	C = g_x*g_x + g_y*g_y + g_z*g_z - cr*cr;

	del = B*B - 4*A*C;
 
	if (del > 0) { /** If delta is less than 0 then photon path is not intersecting the object **/ 
		srDel = sqrt(del);
		rot1 = (-B + srDel)/(2*A);
		rot2 = (-B - srDel)/(2*A);

		X1 = rx + rot1*ux;
		Y1 = ry + rot1*uy;
		Z1 = rz + rot1*uz;

		X2 = rx + rot2*ux;
		Y2 = ry + rot2*uy;
		Z2 = rz + rot2*uz;

		u1 = X1 - rx;
		v1 = Y1 - ry;
		w1 = Z1 - rz;
		mag1 = sqrt(u1*u1 + v1*v1 + w1*w1);
		if (mag1 != 0)	{ /** If the photon is not at the intersection1, find angle and distance between photon and intersection1 **/ 
			alpha1 = u1/mag1;
			beta1 = v1/mag1;
			gamma1 = w1/mag1;
			dB1 = (X1-rx)*(X1-rx) + (Y1-ry)*(Y1-ry) + (Z1-rz)*(Z1-rz);
		}
		else	 {
			dB1 = large;
		}

		u2 = X2 - rx;
		v2 = Y2 - ry;
		w2 = Z2 - rz;
		mag2 = sqrt(u2*u2 + v2*v2 + w2*w2);
		if (mag2 != 0)	{ /** If the photon is not at the intersection2, find angle and distance between photon and intersection2 **/ 
			alpha2 = u2/mag2;
			beta2 = v2/mag2;
			gamma2 = w2/mag2;
			dB2 = (X2-rx)*(X2-rx) + (Y2-ry)*(Y2-ry) + (Z2-rz)*(Z2-rz);
		}
		else {
			dB2 = large;	
		}

		if (fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err 
			&& fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** Both the intersection points are in the direction of photon **/
		{
			if (Photon_Ptr->inObj == 0) dl_b = sqrt(min(dB1,dB2)); /** For photon outside object **/
			else
				dl_b = sqrt(max(dB1,dB2)); /** For small difference in the boundary and photon location **/
		}
		else if (mag1 !=0 && fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err) /** photon within the cylinder **/
		{
			if (Photon_Ptr->inObj != 0) dl_b = sqrt(dB1);
		}
		else if (mag2 !=0 && fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** photon within the cylinder **/
		{
			if (Photon_Ptr->inObj != 0) dl_b = sqrt(dB2);
		}
	
		if(dl_b != -1) {
			if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
				if (Photon_Ptr->inObj == 1)
					mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
				else
					mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

				Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
				Photon_Ptr->s = dl_b;
				return(1);
			}
		}
		newx = Photon_Ptr->x + s*Photon_Ptr->ux;
        newy = Photon_Ptr->y + s*Photon_Ptr->uy;
        newz = Photon_Ptr->z + s*Photon_Ptr->uz;

        if(Photon_Ptr->inObj == 0) {//photon currently outside
			if ((newy-cy)*(newy-cy)+(newz-cz)*(newz-cz) < cr*cr) { //future position inside
				Photon_Ptr->inObj = 1;
            }
        } else {//photon currently inside
			if ((newy-cy)*(newy-cy)+(newz-cz)*(newz-cz) > cr*cr) {//future position outside
				Photon_Ptr->inObj = 0;
            }
        }
		return(0);
	}

	return(0);
}

/************************************
** Check if the photon hits the ellipsoid boundary
*******/

Boolean HitEll(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, erx, ery, erz;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double X1, Y1, Z1, X2, Y2, Z2;
	double dB1, dB2;
	double alpha1 = 0, alpha2 = 0.0, beta1 = 0, beta2 = 0.0, gamma1 = 0.0, gamma2 = 0.0;
	double u1, v1, w1, mag1, u2, v2, w2, mag2;
	double mut;
	double err =  1E-6;		/** Very small error **/
	double large = 100000;	/** Large number**/

	double s = Photon_Ptr->s;
	double newx, newy, newz;

	cx = In_Ptr->ObjSpecs[0].cx; 
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	erx = In_Ptr->ObjSpecs[0].crx;
	ery = In_Ptr->ObjSpecs[0].cry;
	erz = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference http://www.ogre3d.org/forums/viewtopic.php?f=2&t=26442 **/

	A = (ux*ux)/(erx*erx) + (uy*uy)/(ery*ery) + (uz*uz)/(erz*erz);
	B = ((2*ux*(rx-cx)/(erx*erx)) + (2*uy*(ry-cy)/(ery*ery)) + (2*uz*(rz - cz)/(erz*erz)));
	C = (rx-cx)*(rx-cx)/(erx*erx) + (ry-cy)*(ry-cy)/(ery*ery) + (rz - cz)*(rz - cz)/(erz*erz) - 1;

	del = B*B - 4*A*C;

	if (del > 0) { /** If delta is less than 0 then photon path is not intersecting the object **/ 

		srDel = sqrt(del);
		rot1 = (-B + srDel)/(2*A);
		rot2 = (-B - srDel)/(2*A);

		X1 = rx + rot1*ux;
		Y1 = ry + rot1*uy;
		Z1 = rz + rot1*uz;

		X2 = rx + rot2*ux;
		Y2 = ry + rot2*uy;
		Z2 = rz + rot2*uz;

		u1 = X1 - rx;
		v1 = Y1 - ry;
		w1 = Z1 - rz;
		mag1 = sqrt(u1*u1 + v1*v1 + w1*w1);
		if (mag1 != 0)	{ /** If the photon is not at the intersection1, find angle and distance between photon and intersection1 **/ 
			alpha1 = u1/mag1;
			beta1 = v1/mag1;
			gamma1 = w1/mag1;
			dB1 = (X1-rx)*(X1-rx) + (Y1-ry)*(Y1-ry) + (Z1-rz)*(Z1-rz);
		}
		else {
			dB1 = large;
		}

		u2 = X2 - rx;
		v2 = Y2 - ry;
		w2 = Z2 - rz;
		mag2 = sqrt(u2*u2 + v2*v2 + w2*w2);
		if (mag2 != 0)	{ /** If the photon is not at the intersection2, find angle and distance between photon and intersection2 **/ 
			alpha2 = u2/mag2;
			beta2 = v2/mag2;
			gamma2 = w2/mag2;
			dB2 = (X2-rx)*(X2-rx) + (Y2-ry)*(Y2-ry) + (Z2-rz)*(Z2-rz);
		}
		else
		{
			dB2 = large;	
		}

		if (fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err 
			&& fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** Both the intersection points are in the direction of photon **/
		{
			if (Photon_Ptr->inObj == 0) dl_b = sqrt(min(dB1,dB2)); /** For photon outside ellipsoid **/
			else dl_b = sqrt(max(dB1,dB2)); /** For small difference in the boundary and photon location **/
		}
		else if (mag1 !=0 && fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err) /** photon within the ellipsoid **/
		{
			if (Photon_Ptr->inObj != 0) dl_b = sqrt(dB1);
		}
		else if (mag2 !=0 && fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** photon within the ellipsoid **/
		{
			if (Photon_Ptr->inObj != 0) dl_b = sqrt(dB2);
		}

		if(dl_b != -1) {
			if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
				if (Photon_Ptr->inObj == 1)
					mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
				else
					mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

				Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
				Photon_Ptr->s = dl_b;
				return(1);
			}
		}
		
		newx = Photon_Ptr->x + s*Photon_Ptr->ux;
        newy = Photon_Ptr->y + s*Photon_Ptr->uy;
        newz = Photon_Ptr->z + s*Photon_Ptr->uz;

        if(Photon_Ptr->inObj == 0) {//photon currently outside
			if ((((newx-cx)*(newx-cx)/(erx*erx))+((newy-cy)*(newy-cy)/(ery*ery))+((newz-cz)*(newz-cz)/((erz*erz)))-1) < 0.0) { //future position inside
				Photon_Ptr->inObj = 1;
            }
        } else {//photon currently inside
			if ((((newx-cx)*(newx-cx)/(erx*erx))+((newy-cy)*(newy-cy)/(ery*ery))+((newz-cz)*(newz-cz)/((erz*erz)))-1) > 0.0) {//future position outside
				Photon_Ptr->inObj = 0;
            }
        }
		
		return(0);
	}

	return(0);
}


/***********************************************************
**	Check if the given distance moved by a photon 
**	lies on the surface of box.
**	Return 1 if hits surface.
**	Return 0 otherwise.
*****/

Boolean checkOnPlane(double bXMin, double bYMin, double bZMin, double bXMax, double bYMax, double bZMax,
					 double x, double y, double z, double ux, double uy, double uz, double dist) 
{
	double newX, newY, newZ;
	newX = (x + ux*dist);	/*Find the point after the photon travels the distance */
	newY = (y + uy*dist);
	newZ = (z + uz*dist);
	if((bXMin <= newX && newX <= bXMax) && /* Check if the point is on box surface */
	   (bYMin <= newY && newY <= bYMax) &&
	   (bZMin <= fabs(newZ) && newZ <= bZMax)) {
		   return(1);
	}
	else return(0);
}
/***********************************************************
**	Check if the step will hit the surface of box.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
**	The minimum distance to boundary is saved.
****/

Boolean hitPla(int cubNum, InputStruct  *  In_Ptr, 
			   PhotonStruct *  Photon_Ptr, double * minDistBox) 
{
	double cXMin, cYMin, cZMin, cXMax, cYMax, cZMax;
	double x, y, z, ux, uy, uz;
	double distX1, distX2, distY1, distY2, distZ1, distZ2;
	double minDistX = 9999, minDistY = 9999, minDistZ = 9999;
	float newX, newY, newZ;
	double minDist, tempMin;
	double ERR = 1E-10, HIGHVAL = 99999;

	/** Get the box dimensions **/
	cXMin = In_Ptr->ObjSpecs[cubNum].XMin;
	cYMin = In_Ptr->ObjSpecs[cubNum].YMin;
	cZMin = In_Ptr->ObjSpecs[cubNum].ZMin;
	cXMax = In_Ptr->ObjSpecs[cubNum].XMax;
	cYMax = In_Ptr->ObjSpecs[cubNum].YMax;
	cZMax = In_Ptr->ObjSpecs[cubNum].ZMax;

	x = Photon_Ptr->x;
	y = Photon_Ptr->y;
	z = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/**	Find the distance of photon with respect to 6 directions. Discard the negative distance	**/
	if(ux != 0) {
		distX1 = (cXMin-x)/ux;
		if(distX1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distX1)){
			distX1 = HIGHVAL;		/** Setting a high value because point doesn't lie on surface with the distance **/
		}
		distX2 = (cXMax-x)/ux;
		if(distX2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distX2)){
			distX2=HIGHVAL;
		}
		minDistX = min(distX1, distX2);
	}

	if(uy != 0) {
		distY1 = (cYMin-y)/uy;
		if(distY1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distY1)){
			distY1 = HIGHVAL;
		}
		distY2 = (cYMax-y)/uy;
		if(distY2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distY2)){
			distY2 = HIGHVAL;
		}
		minDistY = min(distY1, distY2);
	}

	if (uz !=0) {
		distZ1 = (cZMin-z)/uz;
		if(distZ1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distZ1)){
			distZ1 = HIGHVAL;
		}
		distZ2 = (cZMax-z)/uz;
		if(distZ2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distZ2)){
			distZ2 = HIGHVAL;
		}
		minDistZ = min(distZ1, distZ2);
	}

	if (minDistX == HIGHVAL && minDistY == HIGHVAL && minDistZ == HIGHVAL) {
		return 0; /** No hit when there is no valid distance. **/
	}

	minDist=min(min(minDistX, minDistY), minDistZ); /** Minimum of the distance in x, y, and z directions **/

	if (minDist <= Photon_Ptr->s) {
		*minDistBox = minDist;
		return(1);		/** Hit if minimum distance is less than step-size **/
	}
	return(0);
}

/***********************************************************
**	Check if the step will hit the surface of cuboid.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
*
**	If the projected step hits the boundary, the members
**	s and sleft of Photon_Ptr are updated.
****/
Boolean HitCub(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	Boolean hitPlane = 0;
	double hitDist;
	double mut;

	hitPlane = hitPla(0, In_Ptr, Photon_Ptr, &hitDist);		/** Check if photon hits cuboid **/

	/**	Update the photon pointer when there is a hit. **/
	if(hitPlane == 1) {
		if(hitDist <= Photon_Ptr->s) {
			if (Photon_Ptr->inObj == 1) {
				mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
			}
			else {
				mut = In_Ptr->layerspecs[Photon_Ptr->layer].mua + In_Ptr->layerspecs[Photon_Ptr->layer].mus;
			}

			Photon_Ptr->sleft = (Photon_Ptr->s - hitDist)*mut;
			Photon_Ptr->s = hitDist;
			return(1);
		}
	}
	else{
		return(0);
	}
}

/***********************************************************
*	Drop photon weight inside the tissue (not glass).
*
*  The photon is assumed not dead. 
*
*	The weight drop is dw = w*mua/(mua+mus).
*
*	The dropped weight is assigned to the absorption array 
*	elements.
****/
void Drop(InputStruct  *	In_Ptr, 
		  PhotonStruct *	Photon_Ptr,
		  OutStruct * Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double izd, ird;	/* LW 5/20/98. To avoid out of short range.*/
	short  ix, iy, iz, ir;	/* index to z & r. */
	short  layer = Photon_Ptr->layer;
	double mua, mus;		

	/* compute array indices. */

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd>In_Ptr->nz-1) iz=In_Ptr->nz-1;
	else iz = izd;

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	/* update photon weight. */
	mua = In_Ptr->layerspecs[layer].mua;
	mus = In_Ptr->layerspecs[layer].mus;
	dwa = Photon_Ptr->w * mua/(mua+mus);
	Photon_Ptr->w -= dwa;

	/* assign dwa to the absorption array element. */
	Out_Ptr->A_rz[ir][iz] += dwa;
}

/***********************************************************
**	Drop photon weight inside the object.
*
**  The photon is assumed not dead. 
*
**	The weight drop is dw = w*mua/(mua+mus).
*
**	The dropped weight is assigned to the absorption variable 
**	elements in rz coordinates and object absorbance variable.
****/
void DropObj(InputStruct  *	In_Ptr, 
			 PhotonStruct *	Photon_Ptr,
			 OutStruct *	 Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double izd, ird;	/* LW 5/20/98. To avoid out of short range.*/
	short  iz, ir;	/* index to z & r. */
	short  layer = Photon_Ptr->layer;
	double mua, mus;	

	/* compute array indices. */

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz - 1) iz = In_Ptr->nz - 1;
	else iz = izd;

	ird = sqrt(x*x + y*y)/In_Ptr->dr;
	if(ird > In_Ptr->nr - 1) ir = In_Ptr->nr - 1;
	else ir = ird;

	/** update photon weight. **/
	mua = In_Ptr->ObjSpecs[0].mua;
	mus = In_Ptr->ObjSpecs[0].mus;
	dwa = Photon_Ptr->w * mua/(mua+mus);
	Photon_Ptr->w -= dwa;

	Out_Ptr->A_rz[ir][iz] += dwa;	/* assign dwa to the absorption array element. */

	Out_Ptr->A_obj += dwa;	/** assign dwa to the absorption variable within object. **/
}

/***********************************************************
*	The photon weight is small, and the photon packet tries 
*	to survive a roulette.
****/
void Roulette(PhotonStruct * Photon_Ptr)
{
	if(Photon_Ptr->w == 0.0)	
		Photon_Ptr->dead = 1;
	else if(RandomNum() < CHANCE) /* survived the roulette.*/
		Photon_Ptr->w /= CHANCE;
	else 
		Photon_Ptr->dead = 1;
}


/***********************************************************
*	Record the photon weight exiting the first layer(uz<0), 
*	no matter whether the layer is glass or not, to the 
*	reflection array.
**	Also check if the photon is within collecting fiber
**	If it is withing collecting fiber, then record it as
**	class I or class II photon.
*
*	Update the photon weight as well.
****/
void RecordR(double                     Refl,   /* reflectance. */
	InputStruct  * In_Ptr,
	PhotonStruct * Photon_Ptr,
	OutStruct *    Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double z;
	int i;
	short  ir, ia, iz;      /* index to r & angle. */
	double ird, iad, izd;   /* LW 5/20/98. To avoid out of short range.*/
	double lc = In_Ptr->lc; /** coherence length **/
	double theta_max = In_Ptr->theta_max; /** detector receive angle in degree **/
	double r_max = In_Ptr->r_max; /** collecting radius **/
	double probe_wind; /** Computed probing window as per Dr Wang **/

	probe_wind=(lc/(2*In_Ptr->layerspecs[Photon_Ptr->layer].n)); /** Included by Vijitha et al; **/

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	iad = acos(-Photon_Ptr->uz)/In_Ptr->da;
	if(iad>In_Ptr->na-1) ia=In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the reflection array element. */
	Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);

	for(i=0;i<In_Ptr->nz;i++) {
		z = (i+1)*In_Ptr->dz-In_Ptr->dz/2;
		if(sqrt((Photon_Ptr->launchLoc-x)*(Photon_Ptr->launchLoc-x) + y*y)<r_max && acos(-Photon_Ptr->uz) < (theta_max/180*PI) && fabs(Photon_Ptr->opLen-2*z) < probe_wind) {			
			if(fabs(Photon_Ptr->opLen-(2*Photon_Ptr->zMax)) < probe_wind) { /** Type 1 photons **/
				Out_Ptr->Rd_z1[i] += Photon_Ptr->biasfact*Photon_Ptr->w; /** Recording of A-line class I photons **/
				Out_Ptr->Rd_z1_bs[i][In_Ptr->launchCnt] = Out_Ptr->Rd_z1_bs[i][In_Ptr->launchCnt] + Photon_Ptr->biasfact*Photon_Ptr->w;	/** Recording for B-scan **/
			} else if (fabs(Photon_Ptr->opLen-(2*Photon_Ptr->zMax)) > probe_wind) { /** Type 2 photons **/
				Out_Ptr->Rd_z2[i] += Photon_Ptr->biasfact*Photon_Ptr->w; /** Recording for A-scan **/
				Out_Ptr->Rd_z2_bs[i][In_Ptr->launchCnt] += Photon_Ptr->biasfact*Photon_Ptr->w;	/** Recording for B-scan **/
			}
		}
	}

	Photon_Ptr->w *= Refl;
}




/***********************************************************
*	Record the photon weight exiting the last layer(uz>0), 
*	no matter whether the layer is glass or not, to the 
*	transmittance array.
*
*	Update the photon weight as well.
****/
void RecordT(double 		    Refl,
			 InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 OutStruct    *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ir, ia;	/* index to r & angle. */
	double ird, iad;	/* LW 5/20/98. To avoid out of short range.*/

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	iad = acos(Photon_Ptr->uz)/In_Ptr->da; /* LW 1/12/2000. Removed -. */
	if(iad>In_Ptr->na-1) ia=In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the transmittance array element. */
	Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);

	Photon_Ptr->w *= Refl;
}

/***********************************************************
*	Decide whether the photon will be transmitted or 
*	reflected on the upper boundary (uz<0) of the current 
*	layer.
*
*	If "layer" is the first layer, the photon packet will 
*	be partially transmitted and partially reflected if 
*	PARTIALREFLECTION is set to 1,
*	or the photon packet will be either transmitted or 
*	reflected determined statistically if PARTIALREFLECTION 
*	is set to 0.
*
*	Record the transmitted photon weight as reflection.  
*
*	If the "layer" is not the first layer and the photon 
*	packet is transmitted, move the photon to "layer-1".
*
*	Update the photon parmameters.
****/
void CrossUpOrNot(InputStruct  *	In_Ptr, 
				  PhotonStruct *	Photon_Ptr,
				  OutStruct	   *Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. always */
	/* positive. */
	double r=0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer-1].n;

	/* Get r. */
	if( - uz <= In_Ptr->layerspecs[layer].cos_crit0) 
		r=1.0;		      /* total internal reflection. */
	else r = RFresnel(ni, nt, -uz, &uz1);

#if PARTIALREFLECTION
	if(layer == 1 && r<1.0) {	/* partially transmitted. */
		Photon_Ptr->uz = -uz1;	/* transmitted photon. */
		RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr);
		Photon_Ptr->uz = -uz;	/* reflected photon. */
	}		
	else if(RandomNum() > r) {/* transmitted to layer-1. */
		Photon_Ptr->layer--;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = -uz1;
	}
	else			      		/* reflected. */
		Photon_Ptr->uz = -uz;
#else
	if(RandomNum() > r) {		/* transmitted to layer-1. */
		if(layer==1)  {
			Photon_Ptr->uz = -uz1;
			RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
			Photon_Ptr->dead = 1;
		}
		else {
			Photon_Ptr->layer--;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = -uz1;
		}
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#endif
}



/***********************************************************
*	Decide whether the photon will be transmitted  or be 
*	reflected on the bottom boundary (uz>0) of the current 
*	layer.
*
*	If the photon is transmitted, move the photon to 
*	"layer+1". If "layer" is the last layer, record the 
*	transmitted weight as transmittance. See comments for 
*	CrossUpOrNot.
*
*	Update the photon parmameters.
****/
void CrossDnOrNot(InputStruct  *	In_Ptr, 
				  PhotonStruct *	Photon_Ptr,
				  OutStruct	   *Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. */
	double r=0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer+1].n;

	/* Get r. */
	if( uz <= In_Ptr->layerspecs[layer].cos_crit1) 
		r=1.0;		/* total internal reflection. */
	else r = RFresnel(ni, nt, uz, &uz1);

#if PARTIALREFLECTION	
	if(layer == In_Ptr->num_layers && r<1.0) {
		Photon_Ptr->uz = uz1;
		RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr);
		Photon_Ptr->uz = -uz;
	}
	else if(RandomNum() > r) {/* transmitted to layer+1. */
		Photon_Ptr->layer++;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = uz1;
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#else
	if(RandomNum() > r) {		/* transmitted to layer+1. */
		if(layer == In_Ptr->num_layers) {
			Photon_Ptr->uz = uz1;
			RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
			Photon_Ptr->dead = 1;
		}
		else {
			Photon_Ptr->layer++;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = uz1;
		}
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#endif
}

/***********************************************************
* If uz>0 check for boundary below and if uz<0 check for boundary above
****/
void CrossOrNot(InputStruct  *	In_Ptr, 
			    PhotonStruct *	Photon_Ptr,
				OutStruct    *	Out_Ptr)
{
	if(Photon_Ptr->uz < 0.0)
		CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
	else
		CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
}

/***********************************************************
**	Decide whether the photon will be transmitted or 
**	reflected from sphere to layer or vise-versa
*
**	The photon packet will be either transmitted or 
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotSph(InputStruct  *	In_Ptr, 
				   PhotonStruct *	Photon_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni; 
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, cr;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;


	cx = In_Ptr->ObjSpecs[0].cx; 
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Change the xyz coordinate system to local coordinate system whose z axis align along the surface normal at the point of intersection **/
	if (Photon_Ptr->inObj == 0){
		u1 = cx-rx;
		v1 = cy-ry;
		w1 = cz-rz;
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = rx-cx;
		v1 = ry-cy;
		w1 = rz-cz;
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;


	if(fabs(gamma1) < COSZERO){ 
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit) 
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}

/***********************************************************
**	Decide whether the photon will be transmitted or 
**	reflected from cylinder to layer or vise-versa
*
**	The photon packet will be either transmitted or 
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotCyl(InputStruct  *	In_Ptr, 
				   PhotonStruct *	Photon_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni; 
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, cr;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;
	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;
	if (Photon_Ptr->inObj == 0){
		u1 = rx-rx;
		v1 = cy-ry;
		w1 = cz-rz;
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = rx-rx;
		v1 = ry-cy;
		w1 = rz-cz;
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;

	if(fabs(gamma1) < COSZERO){ 
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit) 
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}

/***********************************************************
**	Decide whether the photon will be transmitted or 
**	reflected from ellipsoid to layer or vise-versa
*
**	The photon packet will be either transmitted or 
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotEll(InputStruct  *	In_Ptr, 
				   PhotonStruct *	Photon_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni; 
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, erx, ery, erz, d;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;


	cx = In_Ptr->ObjSpecs[0].cx; 
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	erx = In_Ptr->ObjSpecs[0].crx;
	ery = In_Ptr->ObjSpecs[0].cry;
	erz = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Change the xyz coordinate system to local coordinate system whose z axis align along the surface normal at the point of intersection **/
	if (Photon_Ptr->inObj == 0){
		u1 = ((rx-cx)/(erx*erx));
		v1 = ((ry-cy)/(ery*ery));
		w1 = ((rz-cz)/(erz*erz));
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = ((rx-cx)/(erx*erx));
		v1 = ((ry-cy)/(ery*ery));
		w1 = ((rz-cz)/(erz*erz));
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;

	if(fabs(gamma1) < COSZERO){ 
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit) 
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}


/***********************************************************
**	Decide whether the photon will be transmitted or 
**	reflected from cuboid to layer or vise-versa
*
**	The photon packet will be either transmitted or 
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotCub(InputStruct  *	In_Ptr, 
				   PhotonStruct *	Photon_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double ux = Photon_Ptr->ux; /* x directional cosine. */
	double uy = Photon_Ptr->uy; /* y directional cosine. */
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double plaAng, critAng, uz1;	/* cosines of transmission alpha. */
	int plane;
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni;
	double nt;

	if (Photon_Ptr->inObj == 0) {
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
	}
	else {
		nt = In_Ptr->layerspecs[layer].n;
		ni = In_Ptr->ObjSpecs[0].n;
	}
	if(Photon_Ptr->x == In_Ptr->ObjSpecs[0].XMin) { 
		plaAng = Photon_Ptr->ux;
		plane = 1;
	}
	else if(Photon_Ptr->x == In_Ptr->ObjSpecs[0].XMax) {
		plaAng = Photon_Ptr->ux;
		plane = 2;
	}
	else if(Photon_Ptr->y == In_Ptr->ObjSpecs[0].YMin) { 
		plaAng = Photon_Ptr->uy;
		plane = 3;
	}
	else if(Photon_Ptr->y == In_Ptr->ObjSpecs[0].YMax) {
		plaAng = Photon_Ptr->uy;
		plane = 4;
	}
	else if(Photon_Ptr->z == In_Ptr->ObjSpecs[0].ZMin) { 
		plaAng = Photon_Ptr->uz;
		plane = 5;
	}
	else {
		plaAng = Photon_Ptr->uz;
		plane = 6;
	}

	if (Photon_Ptr->inObj == 1) {
		critAng = In_Ptr->ObjSpecs[0].cos_crit1;
	}
	else {
		critAng = In_Ptr->ObjSpecs[0].cos_crit0;
	}

	/* Get r. */

	if( fabs(plaAng) <= critAng) 
		r = 1.0;		/* total internal reflection. */
	else r = RFresnel(ni, nt, fabs(plaAng), &uz1);

	if(RandomNum() > r) {
		if (Photon_Ptr->inObj == 1)	 {	/*If photon is in material change to container*/
			Photon_Ptr->inObj = 0;
		}
		else {
			Photon_Ptr->inObj = 1;
		}
		switch (plane)
		{
		case 1: 
			Photon_Ptr->ux = SIGN(ux)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 2:
			Photon_Ptr->ux = SIGN(ux)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 3: 
			Photon_Ptr->uy = SIGN(uy)*uz1;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 4: 
			Photon_Ptr->uy = SIGN(uy)*uz1;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 5: 
			Photon_Ptr->uz = SIGN(uz)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->ux *= ni/nt;
			break;

		case 6: 
			Photon_Ptr->uz = SIGN(uz)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->ux *= ni/nt;
			break;

		default:
			break;
		}
	}
	else {
		switch (plane)
		{
		case 1: 
			Photon_Ptr->ux = -ux;
			break;

		case 2:
			Photon_Ptr->ux = -ux;
			break;

		case 3: 
			Photon_Ptr->uy = -uy;
			break;

		case 4: 
			Photon_Ptr->uy = -uy;
			break;

		case 5: 
			Photon_Ptr->uz = -uz;
			break;

		case 6: 
			Photon_Ptr->uz = -uz;
			break;

		default:
			break;
		}
	}
}


/******************
** Check which object the CrossOrNot should check for
*******************/
void CrossOrNotObj(InputStruct  *  In_Ptr,
				   PhotonStruct *  Photon_Ptr,
				   OutStruct    *  Out_Ptr)
{
	switch(In_Ptr->objCode) {
	case 0:
		CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
		break;
	case 1:
		CrossOrNotSph(In_Ptr, Photon_Ptr, Out_Ptr);
		break;
	case 2:
		CrossOrNotCyl(In_Ptr, Photon_Ptr, Out_Ptr);
		break;
	case 3:
		CrossOrNotEll(In_Ptr, Photon_Ptr, Out_Ptr);
		break;
	case 4:
		CrossOrNotCub(In_Ptr, Photon_Ptr, Out_Ptr);
		break;
	}
}

/***********************************************************
*	Move the photon packet in glass layer.
*	Horizontal photons are killed because they will
*	never interact with tissue again.
****/
void HopInGlass(InputStruct  * In_Ptr,
				PhotonStruct * Photon_Ptr,
				OutStruct    * Out_Ptr)
{
	if(Photon_Ptr->uz == 0.0) { 
		/* horizontal photon in glass is killed. */
		Photon_Ptr->dead = 1;
	}
	else {
		StepSizeInGlass(Photon_Ptr, In_Ptr);
		Hop(Photon_Ptr);
		CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
	}
}

/***********************************************************
*	Set a step size, move the photon, drop some weight, 
*	choose a new photon direction for propagation.  
*
*	When a step size is long enough for the photon to 
*	hit an interface, this step is divided into two steps. 
*	First, move the photon to the boundary free of 
*	absorption or scattering, then decide whether the 
*	photon is reflected or transmitted.
*	Then move the photon in the current or transmission 
*	medium with the unfinished stepsize to interaction 
*	site.  If the unfinished stepsize is still too long, 
*	repeat the above process.  
****/
void HopDropSpinInTissue(InputStruct  *  In_Ptr,
						 PhotonStruct *  Photon_Ptr,
						 OutStruct    *  Out_Ptr)
{
	StepSizeInTissue(Photon_Ptr, In_Ptr);

	if(HitBoundary(Photon_Ptr, In_Ptr)) {
		Hop(Photon_Ptr);	/* move to boundary plane. */
		CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
	}
	else {
		Hop(Photon_Ptr);
		Drop(In_Ptr, Photon_Ptr, Out_Ptr);
		SpinBiased(In_Ptr, 
			Photon_Ptr);
	}
}


/******************
** Check which object the Hit has to be should checked for
*******************/
Boolean HitObj(PhotonStruct * Photon_Ptr, 
			   InputStruct  * In_Ptr)
{
	Boolean hit;

	switch (In_Ptr->objCode) {
	case 0:
		hit=HitBoundary(Photon_Ptr, In_Ptr);
		break;
	case 1:
		hit=HitSph(Photon_Ptr, In_Ptr);
		break;
	case 2:
		hit=HitCyl(Photon_Ptr, In_Ptr);
		break;
	case 3:
		hit=HitEll(Photon_Ptr, In_Ptr);
		break;
	case 4:
		hit=HitCub(Photon_Ptr, In_Ptr);
		break;
	}

	return(hit);
}

/***********************************************************
*	Set a step size, move the photon, drop some weight, 
*	choose a new photon direction for propagation.  
*
*	When a step size is long enough for the photon to 
**	hit an interface or object, this step is divided into two steps. 
*	First, move the photon to the boundary free of 
*	absorption or scattering, then decide whether the 
*	photon is reflected or transmitted.
*	Then move the photon in the current or transmission 
*	medium with the unfinished stepsize to interaction 
*	site.  If the unfinished stepsize is still too long, 
*	repeat the above process.  
****/
void HopDropSpinInTissueObj(InputStruct  *  In_Ptr,
							PhotonStruct *  Photon_Ptr,
							OutStruct    *  Out_Ptr)
{

	if (Photon_Ptr->inObj == 1)	
		StepSizeInObject(Photon_Ptr, In_Ptr); /* Stepsize by mut of object */
	else
		StepSizeInTissue(Photon_Ptr, In_Ptr); /* Stepsize by mut of Photon->layer  */

	if (Photon_Ptr->layer == In_Ptr->objLayer){ /** Check if the photon is in the layer where the object is **/
		if (HitObj(Photon_Ptr, In_Ptr)){ /** Check if photon hits object boundary **/
			Hop(Photon_Ptr);	/* move to object surface. */
			CrossOrNotObj(In_Ptr, Photon_Ptr, Out_Ptr);	
		} else { /** Photon is in the layer of object and doesn't hit object **/
			if(Photon_Ptr->inObj == 0) { /** Check if photon hits layer boundary if it is outside object **/ 
				if(HitBoundary(Photon_Ptr, In_Ptr)) { /* If photon hits layer boundary */
					Hop(Photon_Ptr);	/* move to boundary plane. */
					CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
				} else { /* Hop,drop,spin in layer */
					Hop(Photon_Ptr);
					Drop(In_Ptr, Photon_Ptr, Out_Ptr);
					SpinBiased(In_Ptr, Photon_Ptr);
				}
			} else { /** Photon in object. Hop, drop, spin in object.**/ 
				Hop(Photon_Ptr);
				DropObj(In_Ptr, Photon_Ptr, Out_Ptr);
				SpinBiased(In_Ptr, Photon_Ptr);	
			}
		}
	} else { /** Photon not in the layer of object. Conventional boundary check and hop, drop, spin **/
		if(HitBoundary(Photon_Ptr, In_Ptr)) {
			Hop(Photon_Ptr);	/* move to boundary plane. */
			CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
		} else {
			Hop(Photon_Ptr);
			Drop(In_Ptr, Photon_Ptr, Out_Ptr);
			SpinBiased(In_Ptr, 
				Photon_Ptr);
		}
	}		
}


/***********************************************************
****/
void HopDropSpin(InputStruct  *  In_Ptr,
				 PhotonStruct *  Photon_Ptr,
				 OutStruct    *  Out_Ptr)
{
	short layer = Photon_Ptr->layer;

	if((In_Ptr->layerspecs[layer].mua == 0.0) 
		&& (In_Ptr->layerspecs[layer].mus == 0.0)) 
		/* glass layer. */
		HopInGlass(In_Ptr, Photon_Ptr, Out_Ptr);
	else
		if (In_Ptr->objLayer == 0) /* Conventional boundary check and hop drop spin if the object doesn't exist*/
			HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr);
		else /** Boundary check involves object boundaries. Hop drop spin depends on Photon_Ptr->inObj **/
			HopDropSpinInTissueObj(In_Ptr, Photon_Ptr, Out_Ptr);

	if( Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead) 
		Roulette(Photon_Ptr);
}
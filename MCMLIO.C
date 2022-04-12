/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2016.
*
*	Input/output of data.
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al.
****/

#include "mcml.h"

/***********************************************************
*	Structure used to check against duplicated file names.
****/
struct NameList {
	char name[STRLEN];
	struct NameList * next;
};

typedef struct NameList NameNode;
typedef NameNode * NameLink;


/***********************************************************
*	Center a string according to the column width.
****/
char * CenterStr(short  Wid,
	char * InStr,
	char * OutStr)
{
	size_t nspaces;	/* number of spaces to be filled */
	/* before InStr. */

	nspaces = (Wid - strlen(InStr))/2;
	if(nspaces<0) nspaces = 0;

	strcpy(OutStr, "");
	while(nspaces--)  strcat(OutStr, " ");

	strcat(OutStr, InStr);

	return(OutStr);
}

/**************************************************************************
*	Center a string according to the column width.
****/
#define COLWIDTH 80
void CtrPuts(char *InStr)
{
	short       nspaces;		/* number of spaces to be left-filled. */
	char        outstr[STRLEN];

	nspaces = (COLWIDTH - strlen(InStr)) / 2;
	if (nspaces < 0)
		nspaces = 0;

	strcpy(outstr, "");
	while (nspaces--)
		strcat(outstr, " ");

	strcat(outstr, InStr);

	puts(outstr);
}

/***********************************************************
*	Print some messages before starting simulation.  
*	e.g. author, address, program version, year.
**	Modified by Vijitha et al.
****/
#define COLWIDTH 80
void ShowVersion(char *version)
{
	char str[STRLEN];

	CtrPuts(" ");
	CtrPuts("Monte Carlo Simulation of optical coherence tomography"); 
	CtrPuts("in Multi-layered Turbid Media with Embedded Object");
	CtrPuts(version);
	CtrPuts(" ");
	/** Updated in 2016 **/
	CtrPuts("Vijitha Periyasamy, B.E. (vijitha@ntu.edu.sg)");
	CtrPuts("Manojit Pramanik, Ph.D. (manojit@ntu.edu.sg)");
	CtrPuts("Biomedical Imaging Laboratory, School of Chemical and Biomedical Engineering,");
	CtrPuts("Nanyang Technological University, Singapore 637459\n");
	CtrPuts("The program can be obtained from");
	CtrPuts("https://sites.google.com/site/bilab2012/research/software-download");
	CtrPuts(" ");
	CtrPuts("Please cite the following articles in your publications:");
	CtrPuts("V. Periyasamy and M. Pramanik,");
	CtrPuts("\tImportance sampling-based Monte Carlo simulation of"); 
	CtrPuts("time-domain optical coherence tomography with embedded objects,");
	CtrPuts("Applied Optics 55(15), 2016 (In Press)");
	CtrPuts("and");
	CtrPuts("\tL.-H. Wang, S. L. Jacques, and L.-Q. Zheng, MCML - Monte");
	CtrPuts("\tCarlo modeling of photon transport in multi-layered tissues,");
	CtrPuts("\tComputer Methods and Programs in Biomedicine, 47, 131-146 (1995)\n");
}
#undef COLWIDTH

/***********************************************************
*	Get a filename and open it for reading, retry until 
*	the file can be opened.  '.' terminates the program.
*      
*	If Fname != NULL, try Fname first.
****/
FILE *GetFile(char *Fname)
{
	FILE * file=NULL;
	Boolean firsttime=1;

	do {
		if(firsttime && Fname[0]!='\0') { 
			/* use the filename from command line */
			firsttime = 0;
		}
		else {
			printf("Input filename(or . to exit):");
			scanf("%s", Fname);
			firsttime = 0;
		}

		if(strlen(Fname) == 1 && Fname[0] == '.') 
			exit(1);			/* exit if no filename entered. */

		file = fopen(Fname, "r");
	}  while(file == NULL);

	return(file);
}

/***********************************************************
*	Kill the ith char (counting from 0), push the following 
*	chars forward by one.
****/
void KillChar(size_t i, char * Str)
{
	size_t sl = strlen(Str);

	for(;i<sl;i++) Str[i] = Str[i+1];
}

/***********************************************************
*	Eliminate the chars in a string which are not printing 
*	chars or spaces.
*
*	Spaces include ' ', '\f', '\t' etc.
*
*	Return 1 if no nonprinting chars found, otherwise 
*	return 0.
****/
Boolean CheckChar(char * Str)
{
	Boolean found = 0;	/* found bad char. */
	size_t sl = strlen(Str);
	size_t i=0;

	while(i<sl) 
		if (Str[i]<0 || Str[i]>255)
			nrerror("Non-ASCII file\n");
		else if(isprint(Str[i]) || isspace(Str[i])) 
			i++;
		else {
			found = 1;
			KillChar(i, Str);
			sl--;
		}

		return(found);	
}

/***********************************************************
*	Return 1 if this line is a comment line in which the 
*	first non-space character is "#".
*
*	Also return 1 if this line is space line.
****/
Boolean CommentLine(char *Buf)
{
	size_t spn, cspn;

	spn = strspn(Buf, " \t");	
	/* length spanned by space or tab chars. */

	cspn = strcspn(Buf, "#\n");
	/* length before the 1st # or return. */

	if(spn == cspn) 	/* comment line or space line. */
		return(1);
	else				/* the line has data. */	 
		return(0);		
}

/***********************************************************
*	Skip space or comment lines and return a data line only.
****/
char * FindDataLine(FILE *File_Ptr)
{
	static char buf[STRLEN];	/* LW 1/11/2000. Added static. */

	buf[0] = '\0';
	do {	/* skip space or comment lines. */
		if(fgets(buf, 255, File_Ptr) == NULL)  {
			printf("Incomplete data\n");
			buf[0]='\0';
			break;
		}
		else
			CheckChar(buf);
	} while(CommentLine(buf));

	return(buf);
}

/***********************************************************
*	Skip file version, then read number of runs.
****/
short ReadNumRuns(FILE* File_Ptr)
{
	char buf[STRLEN];
	short n=0;

	FindDataLine(File_Ptr); /* skip file version. */

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading number of runs\n");
	sscanf(buf, "%hd",&n);	        
	return(n);
}


/***********************************************************
*	Read the file name and the file format.
*
*	The file format can be either A for ASCII or B for
*	binary.
****/
void ReadFnameFormat(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in file name and format. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading file name and format.\n");
	sscanf(buf, "%s %c", 
		In_Ptr->out_fname, &(In_Ptr->out_fformat) );
	if(toupper(In_Ptr->out_fformat) != 'B') 
		In_Ptr->out_fformat = 'A';
}


/***********************************************************
*	Read the number of photons.
****/
void ReadNumPhotons(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of photons. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of photons.\n");
	sscanf(buf, "%lfd", &In_Ptr->num_photons);
	if(In_Ptr->num_photons<=0) 
		nrerror("Nonpositive number of photons.\n");
}

/***********************************************************
*	Read the members dz and dr.
****/
void ReadDzDr(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in dz, dr. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading dz, dr.\n");
	sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);
	if(In_Ptr->dz<=0) nrerror("Nonpositive dz.\n");
	if(In_Ptr->dr<=0) nrerror("Nonpositive dr.\n");
}

/***********************************************************
**	Read the members lc, collecting angle and collecting diameter.
****/
void ReadLC_CA_CD(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in dz, dr. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading dz, dr.\n");
	sscanf(buf, "%lf%lf%lf", &In_Ptr->lc, &In_Ptr->theta_max, &In_Ptr->r_max);
	if(In_Ptr->lc<=0) nrerror("Nonpositive coherence length.\n");
	if(In_Ptr->theta_max<=0) nrerror("Nonpositive collecting angle.\n");
	if(In_Ptr->r_max<=0) nrerror("Nonpositive collecting diameter.\n");
	In_Ptr->bias_a = 0.925;	/** As per Biomedical Optics Express 3, 692 (2012) **/
	In_Ptr->bias_p = 0.5;	/** As per Biomedical Optics Express 3, 692 (2012) **/

}

/***********************************************************
**	Read the members x-start, x-end and number of scans.
****/
void ReadXStr_XEnd_numOfScan(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in dz, dr. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading dz, dr.\n");
	sscanf(buf, "%lf%lf%hd", &In_Ptr->x_start, &In_Ptr->x_end, &In_Ptr->num_scan);
	if(In_Ptr->num_scan<=0) nrerror("Nonpositive number of scans.\n");
	In_Ptr->latRes= (fabs(In_Ptr->x_start)+fabs(In_Ptr->x_end))/In_Ptr->num_scan;	/** Compute lateral step size **/	
	In_Ptr->launchCnt=0;	/** Track the launch count. **/
	if (In_Ptr->num_scan == 1) In_Ptr->x_start=0.0;	/** When there is only 1 launch, launch at 0.0 **/
}

/***********************************************************
*	Read the members nz, nr, na.
****/
void ReadNzNrNa(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of dz, dr, da. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of dz, dr, da's.\n");
	sscanf(buf, "%hd%hd%hd", 
		&In_Ptr->nz, &In_Ptr->nr, &In_Ptr->na);
	if(In_Ptr->nz<=0) 
		nrerror("Nonpositive number of dz's.\n");
	if(In_Ptr->nr<=0) 
		nrerror("Nonpositive number of dr's.\n");
	if(In_Ptr->na<=0) 
		nrerror("Nonpositive number of da's.\n");
	In_Ptr->da = 0.5*PI/In_Ptr->na;
}


/***********************************************************
*	Read the number of layers.
****/
void ReadNumLayers(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of layers. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of layers.\n");
	sscanf(buf, "%hd", &In_Ptr->num_layers);
	if(In_Ptr->num_layers<=0) 
		nrerror("Nonpositive number of layers.\n");
}

/***********************************************************
**	Read the object code.
**  (0) No embedded object (1) Sphere (2)Cylinder (3) Ellipsoid (4) Cuboid
**  Included by Vijitha et al.
****/
void ReadObjCode(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in object code. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading object code.\n");
	sscanf(buf, "%hd", &In_Ptr->objCode);
	if(In_Ptr->objCode<0) 
		nrerror("Wrong object code.\n");
	if(In_Ptr->objCode>4) 
		nrerror("Wrong object code.\n");
}

/***********************************************************
**	Read the layer in which object is present.
****/
void ReadObjLayer(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read whether object is present. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading object presence.\n");
	sscanf(buf, "%hd", &In_Ptr->objLayer);
	if(In_Ptr->objLayer > In_Ptr->num_layers) 
		nrerror("Error in object status.\n");
}

/***********************************************************
*	Read the refractive index n of the ambient.
****/
void ReadAmbient(FILE *File_Ptr, 
	LayerStruct * Layer_Ptr,
	char *side)
{
	char buf[STRLEN], msg[STRLEN];
	double n;

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') {
		sprintf(msg, "Rading n of %s ambient.\n", side);
		nrerror(msg);
	}

	sscanf(buf, "%lf", &n );
	if(n<=0) nrerror("Wrong n.\n");
	Layer_Ptr->n = n;
}


/***********************************************************
*	Read the parameters of one layer.
*
*	Return 1 if error detected.
*	Return 0 otherwise.
*
*	*Z_Ptr is the z coordinate of the current layer, which
*	is used to convert thickness of layer to z coordinates
*	of the two boundaries of the layer.
****/
Boolean ReadOneLayer(FILE *File_Ptr, 
	LayerStruct * Layer_Ptr,
	double *Z_Ptr)
{
	char buf[STRLEN], msg[STRLEN];
	double d, n, mua, mus, g;	/* d is thickness. */

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') return(1);	/* error. */

	sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
	if(d<0 || n<=0 || mua<0 || mus<0 || g<0 || g>1) 
		return(1);			/* error. */

	Layer_Ptr->n	= n;
	Layer_Ptr->mua = mua;	
	Layer_Ptr->mus = mus;	
	Layer_Ptr->g   = g;
	Layer_Ptr->z0	= *Z_Ptr;
	*Z_Ptr += d;
	Layer_Ptr->z1	= *Z_Ptr;

	return(0);
}

/***********************************************************
*	Read the parameters of one layer at a time.
****/
void ReadLayerSpecs(FILE *File_Ptr, 
	short Num_Layers,
	LayerStruct ** Layerspecs_PP)
{
	char msg[STRLEN];
	short i=0;
	double z = 0.0;	/* z coordinate of the current layer. */

	/* Allocate an array for the layer parameters. */
	/* layer 0 and layer Num_Layers + 1 are for ambient. */
	*Layerspecs_PP = (LayerStruct *)
		malloc((unsigned) (Num_Layers+2)*sizeof(LayerStruct));
	if (!(*Layerspecs_PP)) 
		nrerror("allocation failure in ReadLayerSpecs()");

	ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top"); 
	for(i=1; i<=Num_Layers; i++)  
		if(ReadOneLayer(File_Ptr, &((*Layerspecs_PP)[i]), &z)) {
			sprintf(msg, "Error reading %hd of %hd layers\n", 
				i, Num_Layers);
			nrerror(msg);
		}
		ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom"); 
}

/***********************************************************
*	Compute the critical angles for total internal
*	reflection according to the relative refractive index
*	of the layer.
*	All layers are processed.
****/
void CriticalAngle( short Num_Layers, 
	LayerStruct ** Layerspecs_PP)
{
	short i=0;
	double n1, n2;

	for(i=1; i<=Num_Layers; i++)  {
		n1 = (*Layerspecs_PP)[i].n;
		n2 = (*Layerspecs_PP)[i-1].n;
		(*Layerspecs_PP)[i].cos_crit0 = n1>n2 ? 
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;

		n2 = (*Layerspecs_PP)[i+1].n;
		(*Layerspecs_PP)[i].cos_crit1 = n1>n2 ? 
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
	}
}

/***********************************************************
**	Read the parameters of one object at a time.
****/
Boolean ReadObj(FILE *File_Ptr, 
	ObjStruct * Obj_Ptr, int objCode)
{
	char buf[STRLEN], msg[STRLEN];
	double n, mua, mus, g;	/** properties of object. **/
	double cx, cy, cz, crx, cry, crz;	/** dimensions of the object **/

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') return(1);	/* error. */

	sscanf(buf, "%lf%lf%lf%lf", &n, &mua, &mus, &g);
	if(n<=0 || mua<0 || mus<0 || g<0 || g>1) 
		return(1);			/** error. **/

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') return(1);	/** error. **/

	sscanf(buf, "%lf", &cz);
	strcpy(buf, FindDataLine(File_Ptr));

	if(objCode <= 2 ) {
		sscanf(buf, "%lf", &crz);
		crx=0;
		cry=0;
	}
	else {
		sscanf(buf, "%lf%lf%lf", &crx, &cry, &crz);
	}

	if(cz<0 || crz<0) 
		return(1);			/** error. **/

	Obj_Ptr->n   = n;
	Obj_Ptr->mua = mua;	
	Obj_Ptr->mus = mus;	
	Obj_Ptr->g   = g;
	Obj_Ptr->cx  = 0;
	Obj_Ptr->cy  = 0;	
	Obj_Ptr->cz  = cz;	
	Obj_Ptr->crx  = crx;
	Obj_Ptr->cry  = cry;
	Obj_Ptr->crz  = crz;

	return(0);
}

/***********************************************************
**	Read the parameters of objects.
****/
void ReadObjSpecs(FILE *File_Ptr, 
	ObjStruct ** objspecs_PP, int objCode)
{
	char msg[STRLEN];
	short i = 0;
	short Num_Obj = 1;

	/** Allocate an array for the object parameters. **/

	*objspecs_PP = (ObjStruct *)
		malloc((unsigned) (1)*sizeof(ObjStruct));
	if (!(*objspecs_PP)) 
		nrerror("allocation failure in ReadObjSpecs()");


	for(i=0; i<Num_Obj; i++)  
		if(ReadObj(File_Ptr, &((*objspecs_PP)[i]), objCode)) {
			sprintf(msg, "Error reading %hd of %hd object\n", 
				i, Num_Obj);
			nrerror(msg);
		}
}


/***********************************************************
**	Compute the critical angles for total internal
**	reflection according to the relative refractive index
**	of the object and its layer.
****/
void CriticalAngleObj(InputStruct * In_ptr, LayerStruct ** Layerspecs_PP, ObjStruct ** Objspecs_PP)
{
	short i = 0;
	short Num_Obj = 1;
	double n1, n2;
	short objLay;

	objLay = In_ptr->objLayer;
	for(i=0; i<Num_Obj; i++)  {
		n1 = (*Layerspecs_PP)[objLay].n;
		n2 = (*Objspecs_PP)[i].n;
		(*Objspecs_PP)[i].cos_crit0 = n1>n2 ? /** cos_crit0 is the critical angle for medium to object **/
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;

		n2 = (*Layerspecs_PP)[objLay].n;
		n1 = (*Objspecs_PP)[i].n;
		(*Objspecs_PP)[i].cos_crit1 = n1>n2 ? /** cos_crit1 is the critical angle for object to medium **/
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
	}
}

/***********************************************************
*	Read in the input parameters for one run.
****/
void ReadParm(FILE* File_Ptr, InputStruct * In_Ptr, ObjStruct * Obj_Ptr)

{
	char buf[STRLEN];
	int objCode;
	In_Ptr->Wth = WEIGHT;

	ReadFnameFormat(File_Ptr, In_Ptr);
	ReadNumPhotons(File_Ptr, In_Ptr);
	ReadDzDr(File_Ptr, In_Ptr);
	ReadNzNrNa(File_Ptr, In_Ptr);
	ReadLC_CA_CD(File_Ptr, In_Ptr); /* for OCT */
	ReadXStr_XEnd_numOfScan(File_Ptr, In_Ptr); /* for OCT B-scan */
	ReadNumLayers(File_Ptr, In_Ptr);
	ReadObjCode(File_Ptr, In_Ptr);
	if(In_Ptr->objCode > 0) {
		ReadObjLayer(File_Ptr, In_Ptr);
	}
	ReadLayerSpecs(File_Ptr, In_Ptr->num_layers, 
		&In_Ptr->layerspecs);
	CriticalAngle(In_Ptr->num_layers, &In_Ptr->layerspecs);

	if (In_Ptr->objLayer != 0 && In_Ptr->objCode != 0) { /** Read the parameters  and calculate critical angle of object if object is present. **/
		objCode=In_Ptr->objCode;
		ReadObjSpecs(File_Ptr, &In_Ptr->ObjSpecs, objCode);
		CriticalAngleObj(In_Ptr, &In_Ptr->layerspecs, &In_Ptr->ObjSpecs);
	}
}

/***********************************************************
*	Return 1, if the name in the name list.
*	Return 0, otherwise.
****/
Boolean NameInList(char *Name, NameLink List)
{
	while (List != NULL) {
		if(strcmp(Name, List->name) == 0) 
			return(1);
		List = List->next;
	};
	return(0);
}

/***********************************************************
*	Add the name to the name list.
****/
void AddNameToList(char *Name, NameLink * List_Ptr)
{
	NameLink list = *List_Ptr;

	if(list == NULL) {	/* first node. */
		*List_Ptr = list = (NameLink)malloc(sizeof(NameNode));
		strcpy(list->name, Name);
		list->next = NULL;
	}
	else {				/* subsequent nodes. */
		/* Move to the last node. */
		while(list->next != NULL)
			list = list->next;

		/* Append a node to the list. */
		list->next = (NameLink)malloc(sizeof(NameNode));
		list = list->next;
		strcpy(list->name, Name);
		list->next = NULL;
	}
}

/***********************************************************
*	Check against duplicated file names.
*
*	A linked list is set up to store the file names used
*	in this input data file.
****/
Boolean FnameTaken(char *fname, NameLink * List_Ptr)
{
	if(NameInList(fname, *List_Ptr))
		return(1);
	else {
		AddNameToList(fname, List_Ptr);
		return(0);
	}
}

/***********************************************************
*	Free each node in the file name list.
****/
void FreeFnameList(NameLink List)
{
	NameLink next;

	while(List != NULL) {
		next = List->next;
		free(List);
		List = next;
	}
}

/***********************************************************
*	Check the input parameters for each run.
****/
void CheckParm(FILE* File_Ptr, InputStruct * In_Ptr, ObjStruct * Obj_Ptr)
{
	short i_run;
	short num_runs;	/* number of independent runs. */
	NameLink head = NULL;
	Boolean name_taken;/* output files share the same */
	/* file name.*/
	char msg[STRLEN];

	num_runs = ReadNumRuns(File_Ptr);
	for(i_run=1; i_run<=num_runs; i_run++)  {
		printf("Checking input data for run %hd\n", i_run);
		ReadParm(File_Ptr, In_Ptr, Obj_Ptr);

		name_taken = FnameTaken(In_Ptr->out_fname, &head);
		if(name_taken) 
			sprintf(msg, "file name %s duplicated.\n", 
			In_Ptr->out_fname);

		free(In_Ptr->layerspecs);
		if(name_taken) nrerror(msg);
	}
	FreeFnameList(head);
	rewind(File_Ptr);
}


/***********************************************************
*	Allocate the arrays in OutStruct for one run, and 
*	array elements are automatically initialized to zeros.
****/
void InitOutputData(InputStruct In_Parm, 
	OutStruct * Out_Ptr, PhotonStruct Photon_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short nl = In_Parm.num_layers;	
	/* remember to use nl+2 because of 2 for ambient. */

	if(nz<=0 || nr<=0 || na<=0 || nl<=0) 
		nrerror("Wrong grid parameters.\n");

	/* Init pure numbers. */
	Out_Ptr->Rsp   = 0.0;
	Out_Ptr->Rd    = 0.0;
	Out_Ptr->A     = 0.0;
	Out_Ptr->A_obj = 0.0;	/** Initialize absorbance within object **/
	Out_Ptr->Tt    = 0.0;

	/* Allocate the arrays and the matrices. */
	Out_Ptr->Rd_ra = AllocMatrix(0,nr-1,0,na-1);
	Out_Ptr->Rd_r  = AllocVector(0,nr-1);
	Out_Ptr->Rd_a  = AllocVector(0,na-1);
	Out_Ptr->Rd_z1_bs  = AllocMatrix(0,nz-1,0,In_Parm.num_scan-1);	/** For OCT B-scan **/
	Out_Ptr->Rd_z2_bs  = AllocMatrix(0,nz-1,0,In_Parm.num_scan-1);	/** For OCT B-scan **/

	Out_Ptr->A_rz  = AllocMatrix(0,nr-1,0,nz-1);
	Out_Ptr->A_z   = AllocVector(0,nz-1);
	Out_Ptr->A_l   = AllocVector(0,nl+1);

	Out_Ptr->Tt_ra = AllocMatrix(0,nr-1,0,na-1);
	Out_Ptr->Tt_r  = AllocVector(0,nr-1);
	Out_Ptr->Tt_a  = AllocVector(0,na-1);

	/** Compute the planes of the cuboid if the object code is 4 or else initialize it to 0. **/
	if (In_Parm.objCode == 4) {
		In_Parm.ObjSpecs[0].XMin = -In_Parm.ObjSpecs[0].crx/2 + In_Parm.ObjSpecs[0].cx;
		In_Parm.ObjSpecs[0].YMin = -In_Parm.ObjSpecs[0].cry/2 + In_Parm.ObjSpecs[0].cy;
		In_Parm.ObjSpecs[0].ZMin = -In_Parm.ObjSpecs[0].crz/2 + In_Parm.ObjSpecs[0].cz;
		In_Parm.ObjSpecs[0].XMax = In_Parm.ObjSpecs[0].crx/2 + In_Parm.ObjSpecs[0].cx;
		In_Parm.ObjSpecs[0].YMax = In_Parm.ObjSpecs[0].cry/2 + In_Parm.ObjSpecs[0].cy;
		In_Parm.ObjSpecs[0].ZMax = In_Parm.ObjSpecs[0].crz/2 + In_Parm.ObjSpecs[0].cz;
	}
	else if (In_Parm.objCode != 0) {
		In_Parm.ObjSpecs[0].XMin = 0.0;
		In_Parm.ObjSpecs[0].YMin = 0.0;
		In_Parm.ObjSpecs[0].ZMin = 0.0;
		In_Parm.ObjSpecs[0].XMax = 0.0;
		In_Parm.ObjSpecs[0].YMax = 0.0;
		In_Parm.ObjSpecs[0].ZMax = 0.0;
	}

	Out_Ptr->Rd_z1  = AllocVector(0,nz-1); /**For OCT type1 photons **/
	Out_Ptr->Rd_z2  = AllocVector(0,nz-1); /**For OCT type2 photons **/
}

/***********************************************************
*	Undo what InitOutputData did.
*  i.e. free the data allocations.
****/
void FreeData(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short nl = In_Parm.num_layers;	
	/* remember to use nl+2 because of 2 for ambient. */

	free(In_Parm.layerspecs);

	FreeMatrix(Out_Ptr->Rd_ra, 0, nr-1, 0, na-1);
	FreeVector(Out_Ptr->Rd_r, 0, nr-1);
	FreeVector(Out_Ptr->Rd_a, 0, na-1);
	FreeVector(Out_Ptr->Rd_z1, 0,nz-1);	/** OCT class I photons A-line **/
	FreeVector(Out_Ptr->Rd_z2, 0,nz-1);	/** OCT class II photons  A-line **/
	FreeMatrix(Out_Ptr->Rd_z1_bs, 0, nz-1, 0, In_Parm.num_scan-1);	/** OCT class I photons B-scan **/
	FreeMatrix(Out_Ptr->Rd_z2_bs, 0, nz-1, 0, In_Parm.num_scan-1);	/** OCT class II photons B-scan **/
	FreeMatrix(Out_Ptr->A_rz, 0, nr-1, 0, nz-1);
	FreeVector(Out_Ptr->A_z, 0, nz-1);
	FreeVector(Out_Ptr->A_l, 0, nl+1);

	FreeMatrix(Out_Ptr->Tt_ra, 0, nr-1, 0, na-1);
	FreeVector(Out_Ptr->Tt_r, 0, nr-1);
	FreeVector(Out_Ptr->Tt_a, 0, na-1);
}

/***********************************************************
*	Get 1D array elements by summing the 2D array elements.
****/
void Sum2DRd(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short ir,ia;
	double sum;

	for(ir=0; ir<nr; ir++)  {
		sum = 0.0;
		for(ia=0; ia<na; ia++) sum += Out_Ptr->Rd_ra[ir][ia];
		Out_Ptr->Rd_r[ir] = sum;
	}

	for(ia=0; ia<na; ia++)  {
		sum = 0.0;
		for(ir=0; ir<nr; ir++) sum += Out_Ptr->Rd_ra[ir][ia];
		Out_Ptr->Rd_a[ia] = sum;
	}

	sum = 0.0;
	for(ir=0; ir<nr; ir++) sum += Out_Ptr->Rd_r[ir];
	Out_Ptr->Rd = sum;
}

/***********************************************************
*	Return the index to the layer according to the index
*	to the grid line system in z direction (Iz).
*
*	Use the center of box.
****/
short IzToLayer(short Iz, InputStruct In_Parm)
{
	short i=1;	/* index to layer. */
	short num_layers = In_Parm.num_layers;
	double dz = In_Parm.dz;

	while( (Iz+0.5)*dz >= In_Parm.layerspecs[i].z1 
		&& i<num_layers) i++;

	return(i);
}

/***********************************************************
*	Get 1D array elements by summing the 2D array elements.
****/
void Sum2DA(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	short ix,iy,iz,ir;
	double sum;

	for(iz=0; iz<nz; iz++)  {
		sum = 0.0;
		for(ir=0; ir<nr; ir++) sum += Out_Ptr->A_rz[ir][iz];
		Out_Ptr->A_z[iz] = sum;
	}

	sum = 0.0;
	for(iz=0; iz<nz; iz++) {
		sum += Out_Ptr->A_z[iz];
		Out_Ptr->A_l[IzToLayer(iz, In_Parm)] 
		+= Out_Ptr->A_z[iz];
	}
	Out_Ptr->A = sum;
}

/***********************************************************
*	Get 1D array elements by summing the 2D array elements.
****/
void Sum2DTt(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short ir,ia;
	double sum;

	for(ir=0; ir<nr; ir++)  {
		sum = 0.0;
		for(ia=0; ia<na; ia++) sum += Out_Ptr->Tt_ra[ir][ia];
		Out_Ptr->Tt_r[ir] = sum;
	}

	for(ia=0; ia<na; ia++)  {
		sum = 0.0;
		for(ir=0; ir<nr; ir++) sum += Out_Ptr->Tt_ra[ir][ia];
		Out_Ptr->Tt_a[ia] = sum;
	}

	sum = 0.0;
	for(ir=0; ir<nr; ir++) sum += Out_Ptr->Tt_r[ir];
	Out_Ptr->Tt = sum;
}

/***********************************************************
*	Scale Rd and Tt properly.
*
*	"a" stands for angle alpha.
****
*	Scale Rd(r,a) and Tt(r,a) by
*      (area perpendicular to photon direction)
*		x(solid angle)x(No. of photons).
*	or
*		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
*	or
*		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
****
*	Scale Rd(r) and Tt(r) by
*		(area on the surface)x(No. of photons).
****
*	Scale Rd(a) and Tt(a) by
*		(solid angle)x(No. of photons).
****
**	Scale Rd_z1 and Rx_z2 by
**		1.0/(No. of photons).
****/
void ScaleRdTt(InputStruct In_Parm, OutStruct *	Out_Ptr)
{
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short nz = In_Parm.nz;
	double dr = In_Parm.dr;
	double da = In_Parm.da;
	short ir,ia,iz;
	double scale1, scale2;

	scale1 = 4.0*PI*PI*dr*sin(da/2)*dr*In_Parm.num_photons;
	/* The factor (ir+0.5)*sin(2a) to be added. */

	for(ir=0; ir<nr; ir++)  
		for(ia=0; ia<na; ia++) {
			scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);
			Out_Ptr->Rd_ra[ir][ia] *= scale2;
			Out_Ptr->Tt_ra[ir][ia] *= scale2;
		}

		scale1 = 2.0*PI*dr*dr*In_Parm.num_photons;  
		/* area is 2*PI*[(ir+0.5)*dr]*dr.*/ 
		/* ir+0.5 to be added. */

		for(ir=0; ir<nr; ir++) {
			scale2 = 1.0/((ir+0.5)*scale1);
			Out_Ptr->Rd_r[ir] *= scale2;
			Out_Ptr->Tt_r[ir] *= scale2;
		}

		scale1  = 2.0*PI*da*In_Parm.num_photons;
		/* solid angle is 2*PI*sin(a)*da. sin(a) to be added. */

		for(ia=0; ia<na; ia++) {
			scale2 = 1.0/(sin((ia+0.5)*da)*scale1);
			Out_Ptr->Rd_a[ia] *= scale2;
			Out_Ptr->Tt_a[ia] *= scale2;
		}

		scale2 = 1.0/(double)In_Parm.num_photons; 
		Out_Ptr->Rd *= scale2;
		Out_Ptr->Tt *= scale2;

		/** Sacle Rd_z **/
		scale1 = 1.0/(double)(In_Parm.num_photons); /* for OCT scaling */
		for(iz = 0; iz < nz; iz++) {
			Out_Ptr->Rd_z1[iz] *= scale1;	/** OCT class I photons A-line **/
			Out_Ptr->Rd_z2[iz] *= scale1;	/** OCT class II photons A-line **/
		}

}

/***********************************************************
*	Scale absorption arrays properly.
****/
void ScaleA(InputStruct In_Parm, OutStruct * Out_Ptr)
{

	short nz  = In_Parm.nz;
	short nr  = In_Parm.nr;
	double dz = In_Parm.dz;
	double dr = In_Parm.dr;
	short nl  = In_Parm.num_layers;
	short ix,iy,iz,ir;
	short il;
	double scale1, scale2; /* Scale 1 for rz and scale 2 for xyz. */

	/* Scale A_rz. */
	scale1 = 2.0*PI*dr*dr*dz*In_Parm.num_photons;	

	/* volume is 2*pi*(ir+0.5)*dr*dr*dz.*/ 
	/* ir+0.5 to be added. */
	for(iz=0; iz<nz; iz++) 
		for(ir=0; ir<nr; ir++) 
			Out_Ptr->A_rz[ir][iz] /= (ir+0.5)*scale1;

	/* Scale A_z. */
	scale1 = 1.0/(dz*In_Parm.num_photons);
	for(iz=0; iz<nz; iz++) { 
		Out_Ptr->A_z[iz] *= scale1;
	}

	/* Scale A_l. Avoid int/int. */
	scale1 = 1.0/(double)In_Parm.num_photons;	
	for(il=0; il<=nl+1; il++)
		Out_Ptr->A_l[il] *= scale1;

	Out_Ptr->A *=scale1;
	Out_Ptr->A_obj *=scale1; /** Scale the absorbance within the object incremented in DropObj function **/
}

/***********************************************************
*	Sum and scale results of current run.
****/
void SumScaleResult(InputStruct In_Parm, 
	OutStruct * Out_Ptr)
{
	/* Get 1D & 0D results. */
	Sum2DRd(In_Parm, Out_Ptr);
	Sum2DA(In_Parm, Out_Ptr);
	Sum2DTt(In_Parm, Out_Ptr);

	ScaleRdTt(In_Parm, Out_Ptr);
	ScaleA(In_Parm, Out_Ptr);
}

/***********************************************************
*	Write the version number as the first string in the 
*	file.
*	Use chars only so that they can be read as either 
*	ASCII or binary.
****/
void WriteVersion(FILE *file, char *Version)
{
	fprintf(file, 
		"%s \t# Version number of the file format.\n\n", 
		Version);
	fprintf(file, "####\n# Data categories include: \n");
	fprintf(file, "# InParm, RAT\n");
	fprintf(file, "# A_l, Absorbance in EO, A_z, Rd_r, Rd_a, Tt_r, Tt_a, \n");
	fprintf(file, "# A_rz, Rd_ra, Tt_ra \n");
	fprintf(file, "# Rd_z1, Rd_z2, Rd_z1_bscan, Rd_z2_bscan \n####\n\n");
}

/***********************************************************
*	Write the input parameters to the file.
****/
void WriteInParm(FILE *file, InputStruct In_Parm)
{
	short i;

	fprintf(file, 
		"InParm \t\t\t# Input parameters. cm is used.\n");

	fprintf(file, 
		"%s \tA\t\t# output file name, ASCII.\n", 
		In_Parm.out_fname);
	fprintf(file, 
		"%lf \t\t\t# No. of photons\n", In_Parm.num_photons);
	fprintf(file, 
		"%G\t%G\t\t# dz, dr [cm]\n", In_Parm.dz,In_Parm.dr);
	fprintf(file, "%hd\t%hd\t%hd\t# No. of dz, dr, da.\n\n", 
		In_Parm.nz, In_Parm.nr, In_Parm.na);

	fprintf(file, 
		"%hd\t\t\t\t\t# Number of layers\n", 
		In_Parm.num_layers);

	fprintf(file, 
		"%hd\t\t\t\t\t# Embedded object code\n", 
		In_Parm.objCode);

	fprintf(file, 
		"%hd\t\t\t\t\t# Layer of object\n", 
		In_Parm.objLayer); /**Write the layer of object. */

	fprintf(file, 
		"#n\tmua\tmus\tg\td\t# One line for each layer\n"); 
	fprintf(file, 
		"%G\t\t\t\t\t# n for medium above\n", 
		In_Parm.layerspecs[0].n); 
	for(i=1; i<=In_Parm.num_layers; i++)  {
		LayerStruct s;
		s = In_Parm.layerspecs[i];
		fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n",
			s.n, s.mua, s.mus, s.g, s.z1-s.z0, i);
	}
	fprintf(file, "%G\t\t\t\t\t# n for medium below\n\n", 
		In_Parm.layerspecs[i].n); 

	if (In_Parm.objLayer != 0 && In_Parm.objCode != 0){ /** Write the object parameters if object is present. **/
		ObjStruct obj;
		obj = In_Parm.ObjSpecs[0];
		fprintf(file, "%G\t%G\t%G\t%G\t\t\t\t\t\t# Object properties \n",
			obj.n, obj.mua, obj.mus, obj.g);

		switch (In_Parm.objCode) { /** Write the parameters specific to the object **/
		case 1: 
			fprintf(file, "%G\t\t\t\t\t\t# Sphere depth \n",
				obj.cz);
			fprintf(file, "%G\t\t\t\t\t# Sphere radius (cr) \n\n\n",
				obj.crz);
			break;
		case 2:
			fprintf(file, "%G\t\t\t\t\t\t# Cylinder depth \n",
				obj.cz);
			fprintf(file, "%G\t\t\t\t\t# Cylinder radius (cr) \n\n\n",
				obj.crz);
			break;
		case 3:
			fprintf(file, "%G\t\t\t\t\t\t# Ellipsoid depth \n",
				obj.cz);
			fprintf(file, "%G\t%G\t%G\t\t\t\t\t# Ellipsoid radius (crx, cry, crz) \n\n\n",
				obj.crx, obj.cry, obj.crz);
			break;
		case 4:
			fprintf(file, "%G\t\t\t\t\t\t# Cuboid depth \n",
				obj.cz);
			fprintf(file, "%G\t%G\t%G\t\t\t\t\t# Cuboid dimensions (Lt, Bt, Ht) \n\n\n",
				obj.crx, obj.cry, obj.crz);
			break;
		}
	}
}

/***********************************************************
*	Write reflectance, absorption, transmission. 
****/
void WriteRAT(FILE * file, OutStruct Out_Parm)
{
	fprintf(file, 
		"RAT #Reflectance, absorption, transmission. \n");
	/* flag. */

	fprintf(file, 
		"%-14.6G \t#Specular reflectance [-]\n", Out_Parm.Rsp);
	fprintf(file, 
		"%-14.6G \t#Diffuse reflectance [-]\n", Out_Parm.Rd);
	fprintf(file, 
		"%-14.6G \t#Absorbed fraction [-]\n", Out_Parm.A);
	fprintf(file, 
		"%-14.6G \t#Transmittance [-]\n", Out_Parm.Tt);

	fprintf(file, "\n");
}

/***********************************************************
**	Write absorption as a function of layer and object. 
****/
void WriteA_layer(FILE * file, 
	short Num_Layers,
	OutStruct Out_Parm)
{
	short i;

	fprintf(file, 
		"A_l #Absorption as a function of layers. [-]\n");
	/* flag. */

	for(i=1; i<=Num_Layers; i++)
		fprintf(file, "%12.4G\n", Out_Parm.A_l[i]);
	fprintf(file, "Absorption in object. [-]\n" );
	fprintf(file, "%12.4G\t\t\n ", Out_Parm.A_obj); /** Write the absorption within object. **/
}

/***********************************************************
*	5 numbers each line.
****/
void WriteRd_ra(FILE * file, 
	short Nr,
	short Na,
	OutStruct Out_Parm)
{
	short ir, ia;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Rd[r][angle]. [1/(cm2sr)].",
		"# Rd[0][0], [0][1],..[0][na-1]",
		"# Rd[1][0], [1][1],..[1][na-1]",
		"# ...",
		"# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
		"Rd_ra");

	for(ir=0;ir<Nr;ir++)
		for(ia=0;ia<Na;ia++)  {
			fprintf(file, "%12.4E ", Out_Parm.Rd_ra[ir][ia]);
			if( (ir*Na + ia + 1)%5 == 0) fprintf(file, "\n");
		}

		fprintf(file, "\n");
}

/***********************************************************
**	Write the class I photons of B-scan	
****/
void WriteRd_z1_bs(FILE * file, 
	short Nz,
	short Ns,
	OutStruct Out_Parm,
	InputStruct In_Parm)
{
	short is, iz;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Rd_z1_bs[depth][scan].",
		"# Rd[0][0], [0][1],..[0][ns-1]",
		"# Rd[1][0], [1][1],..[1][ns-1]",
		"# ...",
		"# Rd[nz-1][0], [nz-1][1],..[nz-1][ns-1]",
		"Rd_z1_bs");
	fprintf(file, "%hd \t %hd\n", In_Parm.nz, In_Parm.num_scan);
	for(iz=0;iz<Nz;iz++){
		for(is=0;is<Ns;is++)  {
			fprintf(file, "%12.4E ", Out_Parm.Rd_z1_bs[iz][is]);
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

/***********************************************************
**	Write the class II photons of B-scan	
****/
void WriteRd_z2_bs(FILE * file, 
	short Nz,
	short Ns,
	OutStruct Out_Parm,
	InputStruct In_Parm)
{
	short is, iz;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Rd_z2_bs[depth][scan].",
		"# Rd[0][0], [0][1],..[0][ns-1]",
		"# Rd[1][0], [1][1],..[1][ns-1]",
		"# ...",
		"# Rd[nz-1][0], [nz-1][1],..[nz-1][ns-1]",
		"Rd_z2_bs");
	fprintf(file, "%hd \t %hd\n", In_Parm.nz, In_Parm.num_scan);
	for(iz=0;iz<Nz;iz++){
		for(is=0;is<Ns;is++)  {
			fprintf(file, "%12.4E ", Out_Parm.Rd_z2_bs[iz][is]);
			//if( (iz*Nz + is + 1)%5 == 0) fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

/***********************************************************
*	1 number each line.
****/
void WriteRd_r(FILE * file, 
	short Nr,
	OutStruct Out_Parm)
{
	short ir;

	fprintf(file, 
		"Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");	/* flag. */

	for(ir=0;ir<Nr;ir++) {
		fprintf(file, "%12.4E\n", Out_Parm.Rd_r[ir]);
	}

	fprintf(file, "\n");
}

/***********************************************************
*	1 number each line.
****/
void WriteRd_a(FILE * file, 
	short Na,
	OutStruct Out_Parm)
{
	short ia;

	fprintf(file, 
		"Rd_a #Rd[0], [1],..Rd[na-1]. [sr-1]\n");	/* flag. */

	for(ia=0;ia<Na;ia++) {
		fprintf(file, "%12.4E\n", Out_Parm.Rd_a[ia]);
	}

	fprintf(file, "\n");
}

/***********************************************************
**	1 number each line.
****/
void WriteRd_z(FILE * file, 
	short Nz,
	OutStruct Out_Parm)
{
	short iz;

	fprintf(file, 
		"Rd_z1 #Rd[0], [1],..Rd[nz-1]. [1/(cm)]\n");	/* flag. */

	for(iz=0;iz<Nz;iz++) {
		fprintf(file, "%12.4E\n", Out_Parm.Rd_z1[iz]);
	}

	fprintf(file, "\n");

	fprintf(file, 
		"Rd_z2 #Rd[0], [1],..Rd[nz-1]. [1/(cm)]\n");	/* flag. */

	for(iz=0;iz<Nz;iz++) {
		fprintf(file, "%12.4E\n", Out_Parm.Rd_z2[iz]);
	}

	fprintf(file, "\n");
}
/***********************************************************
*	5 numbers each line.
****/
void WriteTt_ra(FILE * file, 
	short Nr,
	short Na,
	OutStruct Out_Parm)
{
	short ir, ia;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Tt[r][angle]. [1/(cm2sr)].",
		"# Tt[0][0], [0][1],..[0][na-1]",
		"# Tt[1][0], [1][1],..[1][na-1]",
		"# ...",
		"# Tt[nr-1][0], [nr-1][1],..[nr-1][na-1]",
		"Tt_ra");

	for(ir=0;ir<Nr;ir++)
		for(ia=0;ia<Na;ia++)  {
			fprintf(file, "%12.4E ", Out_Parm.Tt_ra[ir][ia]);
			if( (ir*Na + ia + 1)%5 == 0) fprintf(file, "\n");
		}

		fprintf(file, "\n");
}

/***********************************************************
*	5 numbers each line.
****/
void WriteA_rz(FILE * file, 
	short Nr,
	short Nz,
	OutStruct Out_Parm)
{
	short iz, ir;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n", /* flag. */
		"# A[r][z]. [1/cm3]",
		"# A[0][0], [0][1],..[0][nz-1]",
		"# A[1][0], [1][1],..[1][nz-1]",
		"# ...",
		"# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
		"A_rz");

	for(ir=0;ir<Nr;ir++)
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E ", Out_Parm.A_rz[ir][iz]);
			if( (ir*Nz + iz + 1)%5 == 0) fprintf(file, "\n");
		}

		fprintf(file, "\n");

}

/***********************************************************
*	1 number each line.
****/
void WriteA_z(FILE * file, 
	short Nz,
	OutStruct Out_Parm)
{
	short iz;

	fprintf(file, 
		"A_z #A[0], [1],..A[nz-1]. [1/cm]\n");	/* flag. */

	for(iz=0;iz<Nz;iz++) {
		fprintf(file, "%12.4E\n", Out_Parm.A_z[iz]);
	}

	fprintf(file, "\n");

}

/***********************************************************
*	1 number each line.
****/
void WriteTt_r(FILE * file, 
	short Nr,
	OutStruct Out_Parm)
{
	short ir;

	fprintf(file, 
		"Tt_r #Tt[0], [1],..Tt[nr-1]. [1/cm2]\n"); /* flag. */

	for(ir=0;ir<Nr;ir++) {
		fprintf(file, "%12.4E\n", Out_Parm.Tt_r[ir]);
	}

	fprintf(file, "\n");
}

/***********************************************************
*	1 number each line.
****/
void WriteTt_a(FILE * file, 
	short Na,
	OutStruct Out_Parm)
{
	short ia;

	fprintf(file, 
		"Tt_a #Tt[0], [1],..Tt[na-1]. [sr-1]\n"); /* flag. */

	for(ia=0;ia<Na;ia++) {
		fprintf(file, "%12.4E\n", Out_Parm.Tt_a[ia]);
	}

	fprintf(file, "\n");
}

/***********************************************************
****/
void WriteResult(InputStruct In_Parm, 
	OutStruct Out_Parm,
	char * TimeReport)
{
	FILE *file;

	file = fopen(In_Parm.out_fname, "w");
	if(file == NULL) nrerror("Cannot open file to write.\n");

	if(toupper(In_Parm.out_fformat) == 'A') 
		WriteVersion(file, "A1");
	else 
		WriteVersion(file, "B1");

	fprintf(file, "# %s", TimeReport);
	fprintf(file, "\n");

	WriteInParm(file, In_Parm);
	WriteRAT(file, Out_Parm);	
	/* reflectance, absorption, transmittance. */

	/* 1D arrays. */
	WriteA_layer(file, In_Parm.num_layers, Out_Parm);
	WriteA_z(file, In_Parm.nz, Out_Parm);
	WriteRd_r(file, In_Parm.nr, Out_Parm);  
	WriteRd_a(file, In_Parm.na, Out_Parm);
	WriteRd_z(file, In_Parm.nz, Out_Parm);	/**for A-line OCT photons**/
	WriteRd_z1_bs(file, In_Parm.nz,In_Parm.num_scan, Out_Parm, In_Parm);	/**for B-line class I OCT photons**/
	WriteRd_z2_bs(file, In_Parm.nz,In_Parm.num_scan, Out_Parm, In_Parm);	/**for B-line class II OCT photons**/
	WriteTt_r(file, In_Parm.nr, Out_Parm);
	WriteTt_a(file, In_Parm.na, Out_Parm);

	/* 2D arrays. */
	WriteA_rz(file, In_Parm.nr, In_Parm.nz, Out_Parm);
	WriteRd_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);
	WriteTt_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);

	fclose(file);
}
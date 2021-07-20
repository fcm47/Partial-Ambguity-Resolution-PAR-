/*------------------------------------------------------------------------------
* PAR_LAMBDA : 	Partial ambiguity fixing technique is
				developed to improve the reliability of AR, involving partial
				ambiguity decorrelation (PAD) and partial ambiguity resolution
				(PAR).
				ref: [1] Wang, J. , &  Feng, Y. . (2013). Reliability of partial ambiguity fixing with multiple gnss constellations. Journal of Geodesy, 87(1), 1-14.
*
* Contact    : Caoming Fan by email: cmfan_1992@foxmail.com

*-----------------------------------------------------------------------------*/
#include "PAR_LAMBDA.h"
#include <stdlib.h>

#define LOOPMAX     	100000           /* maximum count of search loop */
#define PAR_MIN_SAT     5
#define Psuccess        (0.95)
#define QR_VAR1	        5e-4
#define QR_VAR2	        1e-3

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define PAR

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

static double thresar;

static int record_fix_state[512] = {0};
static int fixed = 0;
static int unfixed = 0;

static double prod(const double *val, int n)
{
	int i;
	double r = val[0];
	for (i=1;i<n;i++) r*=val[i];
	return r;
}

static double normcdf(const double x)
{
	static const double cdftable1[] = {0.995338811976281,  0.995365907387537,  0.995392862160578,  0.995419676918124,  0.995446352280779,  0.995472888867033,  0.995499287293263,  0.995525548173743,  0.995551672120642,  0.995577659744029,  0.995603511651879,  0.995629228450073,  0.995654810742405,  0.995680259130583,  0.995705574214235,  0.995730756590911,  0.995755806856087,  0.995780725603171,  0.995805513423503,  0.995830170906360,  
															0.995854698638964,  0.995879097206478,  0.995903367192018,  0.995927509176649,  0.995951523739395,  0.995975411457242,  0.995999172905137,  0.996022808655997,  0.996046319280712,  0.996069705348147,  0.996092967425147,  0.996116106076542,  0.996139121865149,  0.996162015351777,  0.996184787095231,  0.996207437652315,  0.996229967577836,  0.996252377424613,  0.996274667743471,  0.996296839083254,  
															0.996318891990825,  0.996340827011071,  0.996362644686907,  0.996384345559279,  0.996405930167169,  0.996427399047600,  0.996448752735638,  0.996469991764398,  0.996491116665046,  0.996512127966805,  0.996533026196959,  0.996553811880857,  0.996574485541914,  0.996595047701621,  0.996615498879544,  0.996635839593331,  0.996656070358715,  0.996676191689519,  0.996696204097660,  0.996716108093151,  
															0.996735904184109,  0.996755592876756,  0.996775174675425,  0.996794650082564,  0.996814019598740,  0.996833283722642,  0.996852442951088,  0.996871497779027,  0.996890448699542,  0.996909296203860,  0.996928040781350,  0.996946682919528,  0.996965223104068,  0.996983661818797,  0.997001999545704,  0.997020236764945,  0.997038373954848,  0.997056411591910,  0.997074350150813,  0.997092190104418,  
															0.997109931923774,  0.997127576078123,  0.997145123034903,  0.997162573259752,  0.997179927216512,  0.997197185367235,  0.997214348172186,  0.997231416089849,  0.997248389576929,  0.997265269088358,  0.997282055077299,  0.997298747995149,  0.997315348291547,  0.997331856414375,  0.997348272809763,  0.997364597922095,  0.997380832194011,  0.997396976066412,  0.997413029978468,  0.997428994367617,  
															0.997444869669572,  0.997460656318327,  0.997476354746157,  0.997491965383627,  0.997507488659595,  0.997522925001214,  0.997538274833940,  0.997553538581533,  0.997568716666066,  0.997583809507923,  0.997598817525811,  0.997613741136757,  0.997628580756118,  0.997643336797583,  0.997658009673178,  0.997672599793269,  0.997687107566568,  0.997701533400140,  0.997715877699401,  0.997730140868127,  
															0.997744323308458,  0.997758425420902,  0.997772447604339,  0.997786390256025,  0.997800253771600,  0.997814038545087,  0.997827744968899,  0.997841373433847,  0.997854924329136,  0.997868398042379,  0.997881794959595,  0.997895115465216,  0.997908359942089,  0.997921528771486,  0.997934622333101,  0.997947641005060,  0.997960585163925,  0.997973455184694,  0.997986251440810,  0.997998974304167,  
															0.998011624145106,  0.998024201332428,  0.998036706233397,  0.998049139213739,  0.998061500637653,  0.998073790867812,  0.998086010265368,  0.998098159189955,  0.998110237999699,  0.998122247051215,  0.998134186699616,  0.998146057298516,  0.998157859200036,  0.998169592754806,  0.998181258311970,  0.998192856219194,  0.998204386822663,  0.998215850467094,  0.998227247495735,  0.998238578250370,  
															0.998249843071324,  0.998261042297469,  0.998272176266226,  0.998283245313570,  0.998294249774037,  0.998305189980723,  0.998316066265294,  0.998326878957987,  0.998337628387615,  0.998348314881574,  0.998358938765843,  0.998369500364991,  0.998380000002180,  0.998390437999174,  0.998400814676335,  0.998411130352635,  0.998421385345657,  0.998431579971599,  0.998441714545281,  0.998451789380144,  
															0.998461804788262,  0.998471761080340,  0.998481658565720,  0.998491497552388,  0.998501278346975,  0.998511001254763,  0.998520666579687,  0.998530274624346,  0.998539825689997,  0.998549320076568,  0.998558758082660,  0.998568140005548,  0.998577466141189,  0.998586736784226,  0.998595952227991,  0.998605112764508,  0.998614218684501,  0.998623270277395,  0.998632267831324,  0.998641211633129};

	static const double cdftable2[] = {0.998650101968370,  0.998693761551231,  0.998736126572328,  0.998777231306408,  0.998817109256896,  0.998855793168977,  0.998893315042591,  0.998929706145321,  0.998964997025197,  0.998999217523386,  0.999032396786782,  0.999064563280486,  0.999095744800178,  0.999125968484368,  0.999155260826541,  0.999183647687171,  0.999211154305624,  0.999237805311933,  0.999263624738446,  0.999288636031355,  
															0.999312862062084,  0.999336325138560,  0.999359047016340,  0.999381048909613,  0.999402351502066,  0.999422974957609,  0.999442938930975,  0.999462262578170,  0.999480964566793,  0.999499063086214,  0.999516575857616,  0.999533520143892,  0.999549912759408,  0.999565770079618,  0.999581108050550,  0.999595942198136,  0.999610287637418,  0.999624159081600,  0.999637570850967,  0.999650536881662,  
															0.999663070734323,  0.999675185602581,  0.999686894321419,  0.999698209375391,  0.999709142906709,  0.999719706723184,  0.999729912306036,  0.999739770817572,  0.999749293108720,  0.999758489726432,  0.999767370920964,  0.999775946653009,  0.999784226600705,  0.999792220166519,  0.999799936483993,  0.999807384424364,  0.999814572603067,  0.999821509386095,  0.999828202896254,  0.999834661019280,  
															0.999840891409842,  0.999846901497426,  0.999852698492093,  0.999858289390124,  0.999863680979554,  0.999868879845579,  0.999873892375862,  0.999878724765715,  0.999883383023185,  0.999887872974018,  0.999892200266523,  0.999896370376326,  0.999900388611024,  0.999904260114731,  0.999907989872526,  0.999911582714799,  0.999915043321502,  0.999918376226297,  0.999921585820616,  0.999924676357621,  
															0.999927651956075,  0.999930516604120,  0.999933274162970,  0.999935928370511,  0.999938482844817,  0.999940941087581,  0.999943306487466,  0.999945582323366,  0.999947771767598,  0.999949877889004,  0.999951903655982,  0.999953851939444,  0.999955725515688,  0.999957527069211,  0.999959259195441,  0.999960924403402,  0.999962525118309,  0.999964063684097,  0.999965542365885,  0.999966963352371,  
															0.999968328758167,  0.999969640626073,  0.999970900929288,  0.999972111573559,  0.999973274399281,  0.999974391183526,  0.999975463642034,  0.999976493431131,  0.999977482149611,  0.999978431340552,  0.999979342493088,  0.999980217044132,  0.999981056380049,  0.999981861838282,  0.999982634708926,  0.999983376236270,  0.999984087620281,  0.999984770018052,  0.999985424545209,  0.999986052277273,  
															0.999986654250984,  0.999987231465586,  0.999987784884075,  0.999988315434405,  0.999988824010668,  0.999989311474225,  0.999989778654816,  0.999990226351627,  0.999990655334330,  0.999991066344087,  0.999991460094529,  0.999991837272697,  0.999992198539962,  0.999992544532909,  0.999992875864198,  0.999993193123401,  0.999993496877799,  0.999993787673173,  0.999994066034554,  0.999994332466958,  
															0.999994587456092,  0.999994831469043,  0.999995064954938,  0.999995288345588,  0.999995502056111,  0.999995706485530,  0.999995902017353,  0.999996089020140,  0.999996267848040,  0.999996438841320,  0.999996602326875,  0.999996758618713,  0.999996908018431,  0.999997050815677,  0.999997187288588,  0.999997317704220,  0.999997442318961,  0.999997561378926,  0.999997675120350,  0.999997783769952,  
															0.999997887545298,  0.999997986655145,  0.999998081299780,  0.999998171671336,  0.999998257954110,  0.999998340324856,  0.999998418953081,  0.999998494001322,  0.999998565625416,  0.999998633974755,  0.999998699192546,  0.999998761416043,  0.999998820776783,  0.999998877400815,  0.999998931408906,  0.999998982916757,  0.999999032035204,  0.999999078870404,  0.999999123524027,  0.999999166093434,  
															0.999999206671848,  0.999999245348521,  0.999999282208893,  0.999999317334747,  0.999999350804357,  0.999999382692628,  0.999999413071236,  0.999999442008757,  0.999999469570797,  0.999999495820112,  0.999999520816723,  0.999999544618035,  0.999999567278938,  0.999999588851916,  0.999999609387146,  0.999999628932592,  0.999999647534102,  0.999999665235492,  0.999999682078634,  0.999999698103538};

	static double start1,t1, end1, start2, t2, end2;
	int index;
	t1 = 0.002; start1 = 2.6; end1 = 3.0-t1;		// x1 = start1:t1:end1; y1=normcdf(x1);
	t2 = 0.01;   start2 = 3.0; end2 = 5.0-t2;		// x2 = start2:t2:end2; y2=normcdf(x2);
	//
	if (x>=end2)
	{
		return 1.0;
	}
	if (x<start1)
	{
		return 0.995;
	}
	//
	if (x<=end1)
	{
		index = ROUND((x-start1)*500);
		return index==200? cdftable2[0] : cdftable1[index];
	}
	else
	{
		index = ROUND((x-start2)*100);
		return index==200? cdftable2[199] : cdftable2[index];
	}
	
}

// matrix: n*n
static void diag(const double *matrix, int n, double *val)
{
	int i;
	if (matrix==NULL)
	{
		return;
	}
	for (i = 0; i < n; i++)
	{
		val[i] = matrix[i*n+i];
	}
}
// find maximum in s[] and return its index
static int find_max(double *s, int n)
{
	int i, index = -1;
	double MIN = 0.0;
	if(s==NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (s[i]>MIN)
		{
			MIN = s[i];
			index = i;
		}
		
	}
	return index;
}

// find integer val in s[] and return its index
static int findi(const int *s, int n, int val)
{
	int i;
	if (s==NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (s[i]==val)
		{
			return i;
		}
		
	}
	return -1;	
}

// find float val in s[] and return its index
static int findf(const double *s, int n, double val)
{
	int i;
	if (s==NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (fabs(s[i]-val)<1e-8)
		{
			return i;
		}
		
	}
	return -1;	
}

// remove one row and one col at "remove_index" in Q and get new Q
static void remove_mat(const double *Q, int n, int remove_index, double *newQ)
{
	int i, j, row=0;
	int col;
	if (Q==NULL || newQ==NULL || remove_index>=n)
	{
		return;
	}
	for (i = 0; i < n; i++)
	{
		if (i!=remove_index)
		{
			col = 0;
			for (j = 0; j < n; j++)
			{
				if (j!=remove_index)
				{
					newQ[row*(n-1)+col]=Q[i*n+j];
					col++;
				}	
			}
			row++;
		}	
		
	}
	
}

// remove one row at "remove_index" in x and get new x
static void remove_vec(const double *x, int n, int remove_index, double *newx)
{
	int i, row=0;;
	if (x==NULL || newx==NULL || remove_index>=n)
	{
		return;
	}
	for (i = 0; i < n; i++)
	{
		if (i!=remove_index)
		{
			newx[row++] = x[i];
		}
	}
}

// get a block in A. C = A[row_start:(row_start+rows), col_start:(col_start+cols)]
extern void mat_block(const double *A, int n, int m, int row_start, int rows,  int col_start, int cols,  double *C)
{
	int i, j, k=0, q=0;
	if (A==NULL || C==NULL)
	{
		printf("mat_block A==NULL || C==NULL\n");
		return;
	}
	
	for (i=col_start; i<col_start+cols; i++)
	{
		for (j=row_start;j<row_start+rows; j++)
		{
			C[k*rows+q] = A[i*n+j];
			q++;
		}
		q=0;
		k++;
	}
}

//matrix type: int*float
extern void ifmatmul(const char *tr, int n, int k, int m, double alpha,
				   const int *A, const double *B, double beta, double *C)
{
	double d;
	int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

	for (i=0;i<n;i++) for (j=0;j<k;j++) {
		d=0.0;
		switch (f) {
		case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
		case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
		case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
		case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
		}
		if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
	}
}

//matrix type: float*int
extern void fimatmul(const char *tr, int n, int k, int m, double alpha,
					 const double *A, const int *B, double beta, double *C)
{
	double d;
	int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

	for (i=0;i<n;i++) for (j=0;j<k;j++) {
		d=0.0;
		switch (f) {
		case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
		case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
		case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
		case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
		}
		if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
	}
}

////matrix type: int*int
extern void imatmul(const char *tr, int n, int k, int m, int alpha,
					const int *A, const int *B, int beta, int *C)
{
	int d;
	int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

	for (i=0;i<n;i++) for (j=0;j<k;j++) {
		d=0;
		switch (f) {
		case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
		case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
		case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
		case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
		}
		if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
	}
}

static void mat_int_frac(const double *full, int *int_part, double *frac_part, int n)
{
	int i;
	for (i=0; i<n; i++)
	{
		int_part[i] = ROUND(full[i]);
		frac_part[i] = full[i]-int_part[i];
	}
}

static void mat_add(const int *int_part, const double *frac_part, double *full, int n)
{
	int i;
	for (i=0; i<n; i++)
	{
		full[i] = int_part[i] + frac_part[i];
	}
}

extern int *ieye(int n)
{
	int *p;
	int i;

	if ((p=izeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1;
	return p;
}

extern int *izeros(int n, int m)
{
	int *p;

#if NOCALLOC
	if ((p=imat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0;
#else
	if (n<=0||m<=0) return NULL;
	if (!(p=(int *)calloc(sizeof(int),n*m))) {
		fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
	}
#endif
	return p;
}

/* fatal error ---------------------------------------------------------------*/
static void fatalerr(const char *format, ...)
{
    va_list ap;
    va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
    exit(-9);
}

/* test satellite system (m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN) ---------*/
static int test_sys(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
        case SYS_QZS: return m==4;
        case SYS_IRN: return m==5;
    }
    return 0;
}

/* restore SD (single-differenced) ambiguity ---------------------------------*/
static void restamb(rtk_t *rtk, const double *bias, int nb, double *xa)
{
    int i,n,m,f,index[MAXSAT],nv=0,nf=NF(&rtk->opt);
    
    trace(3,"restamb :\n");
    
    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];
    
    for (m=0;m<5;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            if (!test_sys(rtk->ssat[i].sys,m)||rtk->ssat[i].fix[f]!=2) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
        }
        if (n<2) continue;
        
        xa[index[0]]=rtk->x[index[0]];
        
        for (i=1;i<n;i++) {
            xa[index[i]]=xa[index[0]]-bias[nv++];
        }
    }
}

/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}


/* index for SD to DD transformation matrix D --------------------------------*/
// The difference from origin is that adding records of ssat_fix. We can set their ssat_fix if some ambiguities are not fixed.
static int PAR_ddidx(rtk_t *rtk, int *ix)
{
    int i,j,k,m,f,nb=0,na=rtk->na,nf=NF(&rtk->opt),nofix;
	memset(record_fix_state, 0, sizeof(record_fix_state));
    // record_fix_state.res
    trace(3,"ddidx   :\n");
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].fix[j]=0;
    }
    for (m=0;m<6;m++) { /* m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN */
        
        nofix=(m==1&&rtk->opt.glomodear==0)||(m==3&&rtk->opt.bdsmodear==0);
        
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {
            
            for (i=k;i<k+MAXSAT;i++) {
                if (rtk->x[i]==0.0||!test_sys(rtk->ssat[i-k].sys,m)||
                    !rtk->ssat[i-k].vsat[f]||!rtk->ssat[i-k].half[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    rtk->ssat[i-k].fix[f]=2; /* fix */
                    break;
                }
                else rtk->ssat[i-k].fix[f]=1;
            }
            for (j=k;j<k+MAXSAT;j++) {
                if (i==j||rtk->x[j]==0.0||!test_sys(rtk->ssat[j-k].sys,m)||
                    !rtk->ssat[j-k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[j-k].lock[f]>0&&!(rtk->ssat[j-k].slip[f]&2)&&
                    rtk->ssat[i-k].vsat[f]&&
                    rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    ix[nb*2  ]=i; /* state index of ref bias */
                    ix[nb*2+1]=j; /* state index of target bias */
                    rtk->ssat[j-k].fix[f]=2; /* fix */
					record_fix_state[nb] = (j-k)*10+f;
					nb++;
                }
                else rtk->ssat[j-k].fix[f]=1;
            }
        }
    }
    return nb;
}

static int PAR_reset_ssat_fix(rtk_t *rtk, int *index, int unfixed, int nb)
{
	int i, sat_f, id, f, fix_sat_cnt, nf=NF(&rtk->opt);
	for (i = 0; i < unfixed; i++)
	{
		sat_f = record_fix_state[index[i]];
		id = floor(sat_f/10.0);
		f = sat_f-id*10;
		// set ambguity state that has not been fixed to 1
		// considering that restoring SD (single-differenced) ambiguity uses this flag
		rtk->ssat[id].fix[f] = 1;
	}

	// count the satelite number that ambguity fixed
	fix_sat_cnt = 0;
	for (i = 0; i < MAXSAT; i++)
	{
		for (f = 0; f < nf; f++)
		{
			if(rtk->ssat[i].fix[f] == 2)
			{
				fix_sat_cnt++;
				break;
			}
		}
	}
	return fix_sat_cnt;
}

/* The difference from origin is that recording the row transformation, 
and we can track the original ambiguity-index*/

static void PAR_reduction(int n, double *L, double *D, double *Z, int *record)
{
    int i,j,k;
    double del;
    int swap_tmp;
    j=n-2; k=n-2;
	for (i = 0; i < n; i++)
	{
		record[i] = i;
	}
	
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
#ifdef PAR
			// record indices for ahat and Qahat
			swap_tmp = record[j];
			record[j] = record[j+1];
			record[j+1] = swap_tmp;

#endif // PAR
            k=j; j=n-2;
        }
        else j--;
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}

/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}


/* PAR-lambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
*          double *fixed O  number of fixed ambiguity
*          double *record_Z_row O  a matrix, record ambguity index
*          double *delete_record O  record the ambguitiy indexs that have not been fixed
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int PAR_lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s, int *fixed, int *record_Z_row, int *delete_record)
{
    int info, i, j, k, unfixed=0, inc, ind, Qahat_index, ahat_index_old;
    double *L,*D,*Z,*z,*E, *frac_part, *Qahat, *ahat, *Qzhat, *Qvar, *P, *ZTQ;
	double delete_val, Ps;
	int *record, *int_part, *new_int_part;
	int del_num = 0;
	*fixed=0;
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1);E=mat(n,m); frac_part=mat(n, 1); 
	record=imat(n, 1); int_part=imat(n, 1); 
	Qahat = mat(n, n); ahat = mat(n, 1); Qzhat = mat(n, n); Qvar = mat(n, 1);
	ZTQ = mat(n, n);
	P = mat(n, 1);
	// depart integer part and fraction part of float ambiguities.
    mat_int_frac(a, int_part, frac_part, n);
	memcpy(Qahat, Q, n*n*sizeof(double));
	memcpy(ahat, frac_part, n*sizeof(double));
	for (k=0; k<n-PAR_MIN_SAT; k++)
	{
		if (k==0)
		{
			if (!(info=LD(n,Qahat,L,D))) {
				PAR_reduction(n,L,D,Z,record);
				matmul("TN",n,1,n,1.0,Z,ahat,0.0,z); /* z=Z'*a */
			}
			else
			{
				break;
			}
			
		}
		
		if (k>0)
		{
			matmul("TN", n-k+1, n-k+1, n-k+1, 1.0, Z, Qahat, 0.0, ZTQ); // Z'*Qahat
			matmul("NN", n-k+1, n-k+1, n-k+1, 1.0, ZTQ, Z, 0.0, Qzhat);	//Qzhat = Z'*Qahat*Z
			diag(Qzhat, n-k+1, Qvar);	// get diag convariances and store them in Qvar
			ind=find_max(Qvar, n-k+1); 	// get the index of maximum var in Qvar, the index in Qzhat
			Qahat_index = record[ind];	// according to our record,  we can find the index in Qahat, as well as in ahat
			delete_val = ahat[Qahat_index];	// remember the value will be deleted, We will use it to get its index in original ahat
			remove_mat(Qahat, n-k+1, Qahat_index, Qahat);	// remove the ambiguity with max decorrelation convariance, and its convariance
			remove_vec(ahat, n-k+1, Qahat_index, ahat);

			ahat_index_old = findf(frac_part, n, delete_val);	// get its index in original ahat
			delete_record[del_num++] = ahat_index_old;			// store all of the deleted ambiguitiy indexs
			// must initial
			memset(L, 0, n*n*sizeof(double)); //initial zero
			for(i=0; i<n-k; i++) for(j=0;j<n-k; j++) Z[i*(n-k)+j]= (i==j)? 1:0;	// initial identity
			if (!(info=LD(n-k,Qahat,L,D))) {
				PAR_reduction(n-k,L,D,Z,record);
				matmul("TN",n-k,1,n-k,1.0,Z,ahat,0.0,z); /* z=Z'*a */
			}
			else
			{
				break;
			}
			
		}	
		// calculate successful rate
		for (i=k; i<n; i++)
		{
			P[i] = 2*normcdf(0.5/sqrt(D[i]))-1;
		}
		Ps = prod(P+k, n-k);
		if (Ps<Psuccess)
		{
			continue;
		}

		search(n-k, m, L, D, z, E, s);

		
		if (s[0]<=0.0||s[1]/s[0]>=thresar) 
		{
			*fixed=n-k;
			break;
		}
		
	}
	if (*fixed)
	{
		info=solve("T",Z,E,*fixed,m,F); /* F=Z'\E */
		if ((*fixed)!=n)
		{
			inc = 0;
			unfixed = n- (*fixed);
			// generate a row transformation matrix to move the unfixed ambiguities to front, as well as corresponding convariance matrix
			memset(record_Z_row, 0, n*n*sizeof(int));
			for (i = 0; i < unfixed; i++)
			{
				// unfixed ambiguities
				record_Z_row[delete_record[i]*n+inc] = 1;
				inc++;
			}
			for (i = 0; i < n; i++)
			{
				// fixed ambiguities
				if (findi(delete_record, unfixed, i)<0)
				{
					record_Z_row[n*i+inc] = 1;
					inc++;
				}
			}
			// merge integer part and fraction part
			new_int_part = imat(n, 1);
			imatmul("NN",n,1,n,1, record_Z_row, int_part,0, new_int_part);
			for (i=0; i<m; i++)
			{
				mat_add(new_int_part+unfixed, F+i*(*fixed), F+i*(*fixed), n-unfixed);
			}
			free(new_int_part);
		}
		else
		{
			for (i=0; i<m; i++)
			{
				mat_add(int_part, F+i*n, F+i*n, n);
			}
		}
		
		
	}
	
    free(L); free(D); free(Z); free(z); free(E); free(record); free(frac_part); free(int_part);
    free(ahat); free(Qahat); free(Qzhat); free(Qvar); free(P); free(ZTQ);
	return info;
}



/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
int PAR_resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    int i,j,nb,info,nx=rtk->nx,na=rtk->na;
    double *DP,*y,*b,*db,*Qb,*Qab,*QQ,s[2] = {0};
	double *b_p,*y_p,*yb, *Qb_p,*Qab_p,*QabZT, *ZQb, *ZQbZT; 
	double qr[3];
	int *ix, fixed_num, *record_Z_row, *delete_record, unfixed_num;
	char str[64];
	prcopt_t *opt=&rtk->opt;
    thresar = opt->thresar[0];
    trace(3,"resamb_LAMBDA : nx=%d\n",nx);
    
    rtk->sol.ratio=0.0;
    
    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
	time2str(rtk->sol.time,str,3);
	
    /* single to double-difference transformation matrix (D') */

	ix = imat(nx, 2);
	if((nb=PAR_ddidx(rtk, ix))<=0)
	{
		errmsg(rtk,"no valid double-difference\n");
		free(ix);
		return 0;
	}
	y=mat(nb,1); DP=mat(nb,nx-na);
	b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);
	delete_record=imat(nb, 1); record_Z_row=imat(nb, nb);
	/* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
	/* y=D*xc, Qb=D*Qc*D', Qab=Qac*D' */
    for (i=0;i<nb;i++) {
        y[i]=rtk->x[ix[i*2]]-rtk->x[ix[i*2+1]];
    }
    for (j=0;j<nx-na;j++) for (i=0;i<nb;i++) {
        DP[i+j*nb]=rtk->P[ix[i*2]+(na+j)*nx]-rtk->P[ix[i*2+1]+(na+j)*nx];
    }
    for (j=0;j<nb;j++) for (i=0;i<nb;i++) {
        Qb[i+j*nb]=DP[i+(ix[j*2]-na)*nb]-DP[i+(ix[j*2+1]-na)*nb];
    }
    for (j=0;j<nb;j++) for (i=0;i<na;i++) {
        Qab[i+j*na]=rtk->P[i+ix[j*2]*nx]-rtk->P[i+ix[j*2+1]*nx];
    }

	tracet(0, "%s\n", str);
	tracet(0, "float x=\n");
	tracemat(0, rtk->x, 1, na, 16, 6);
	tracet(0,"float amb="); tracemat(0,y,1,nb,16,6);
	tracet(4,"Qb=\n");
	tracemat(0, Qb, nb, nb, 16,6);
	tracet(0,"Qab=\n");
	tracemat(0, Qab, na, nb, 16,6);
	tracemat(0, rtk->Pa, na, na, 16,6);
	
	/* lambda/mlambda integer least-square estimation */
	if (!(info=PAR_lambda(nb,2,y,Qb,b,s, &fixed_num, record_Z_row, delete_record))) {

		// trace(0,"N(1)="); tracemat(0,b   ,1,nb,10,3);
		// trace(0,"index"); traceimat(0,index   ,1,nb,10);
		// trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);

		rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
		if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

		/* validation by popular ratio-test */
		if (s[1]/s[0]>=opt->thresar[0]) {

			for (i=0;i<na;i++) {
				rtk->xa[i]=rtk->x[i];
				for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
			}
			if (nb==fixed_num)
			{
				for (i=0;i<nb;i++) {
					bias[i]=b[i];
					y[i]-=b[i];
				}
				if (!matinv(Qb,nb)) {
					/* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
					matmul("NN",nb,1,nb, 1.0,Qb ,y,0.0,db);
					matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,rtk->xa);

					/* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
					matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
					matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa);

					/* restore single-differenced ambiguity */
					//restamb(rtk,bias,nb,xa);
					for (i=0;i<3;i++) {
						qr[i]=(float)rtk->Pa[i+i*rtk->na];
					}
					if((qr[0]+qr[1]+qr[2])<QR_VAR1 || (rtk->sol.ns>8 && rtk->sol.ratio>10))
					{
						restamb(rtk,bias,nb,xa);
					}

				}
				else nb=0;
			}
			else
			{				
				unfixed_num = nb - fixed_num;
				if (unfixed_num>=2)
				{
					printf("%d\n", unfixed_num);
				}
				
				b_p = mat(nb, 1); y_p = mat(nb, 1); yb = mat(nb, 1);
				Qb_p = mat(fixed_num, fixed_num); Qab_p = mat(na, fixed_num);
				QabZT = mat(na, nb); ZQb = mat(nb, nb); ZQbZT = mat(nb, nb);
				// ifmatmul("NN", nb, 1, nb, 1.0, record_Z_row, b, 0.0, b_p);	// b_p=Zp*b 
				ifmatmul("NN", nb, 1, nb, 1.0, record_Z_row, y, 0.0, y_p);	//y_p=Zp*y, 
				fimatmul("NT", na, nb, nb, 1.0, Qab, record_Z_row, 0.0, QabZT);	//QabZT=Qab*Zp' = (Zp*Qab')'
				//ZQbZT=Zp*Qb*Zp';
				ifmatmul("NN", nb, nb, nb, 1.0, record_Z_row, Qb, 0.0, ZQb);
				fimatmul("NT", nb, nb, nb, 1.0, ZQb, record_Z_row, 0.0, ZQbZT);
				for (i=0;i<fixed_num;i++) {
					yb[i] = y_p[i+unfixed_num] - b[i];
				}
				mat_block(ZQbZT, nb, nb, unfixed_num, fixed_num, unfixed_num, fixed_num, Qb_p);	
				mat_block(QabZT, na, nb, 0, na, unfixed_num, fixed_num, Qab_p);	
				if (!matinv(Qb_p,fixed_num))
				{
					/* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
					matmul("NN",fixed_num,1,fixed_num, 1.0,Qb_p ,yb,0.0,db);		
					matmul("NN",na,1,fixed_num,-1.0,Qab_p,db  ,1.0,rtk->xa);

					/* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
					matmul("NN",na,fixed_num,fixed_num, 1.0,Qab_p,Qb_p ,0.0,QQ);
					matmul("NT",na,na,fixed_num,-1.0,QQ ,Qab_p,1.0,rtk->Pa);
					// 					tracemat(0,ZQbZT, nb, nb, 16, 6);
					// 					tracemat(0,Qab_p, unfixed_num, fixed_num, 16, 6);
					// 					tracemat(0,QQ, na, fixed_num, 16, 6);
					// 					tracemat(0, b_p, nb, 1, 16, 6);
					// 					tracemat(0, y_p, nb, 1, 16, 6);
					//  					tracemat(0, yb, nb, 1, 16, 6);
					// 					tracemat(0, db, fixed_num, 1, 16, 6);
					// trace(4,"N(b)="); tracemat(0, bias, 1, nb, 10,3);
					//restamb(rtk,bias,nb,xa);
					for (i=0;i<3;i++) {
						qr[i]=(float)rtk->Pa[i+i*rtk->na];
					}
					rtk->sol.ns = PAR_reset_ssat_fix(rtk, delete_record, unfixed_num, nb);
					if((qr[0]+qr[1]+qr[2])<QR_VAR1 || ((qr[0]+qr[1]+qr[2])<QR_VAR2 && rtk->sol.ns>6 && rtk->sol.ratio>5))
					{
						restamb(rtk,bias,nb,xa);
						trace(4,"%s resamb : validation ok (nb=%d fixed=%d ratio=%.2f s=%.2f/%.2f)\n",
							str, nb,fixed_num, s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
						nb = fixed_num;
					}
					else
					{
						nb = 0;
					}
					
					
				}
				else 
				{
					fixed_num=0;
					nb = 0;
				}
				free(b_p); free(y_p); free(yb); free(Qb_p); free(Qab_p); free(QabZT); free(ZQb); free(ZQbZT);
			}

		}
		else { /* validation failed */
			errmsg(rtk,"ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
				nb,s[1]/s[0],s[0],s[1]);
			nb=0;
		}
	}
	else {
		errmsg(rtk,"lambda error (info=%d)\n",info);
	}
	free(ix); free(y); free(DP); free(delete_record); free(record_Z_row);
	free(b); free(db); free(Qb); free(Qab); free(QQ);
    return nb; /* number of ambiguities */
}

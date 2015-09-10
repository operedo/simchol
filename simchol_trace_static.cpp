#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <limits.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <string.h>
#include <vector>




//---------------------------------------------------------------------------------------


using namespace std;

#define PI 3.141592654
#define DEG2RAD PI/180.0
#define MAX(x,y) (x>y?x:y)

double EPSILON = 0.000000001;


class Grid
{
public:
	Grid();

	int		nodes[3];
	double	start[3];
	double	size[3];

	int length();
	double *getPoint(int index);
	double **getPoints();

public:
	~Grid();
};

class VariographicModelStructure
{
public:
	VariographicModelStructure(void);

	static const int SphericalModel = 1;
	static const int ExponentialModel = 2;
	static const int GaussModel = 3;
	static const int PowerModel = 4;
	static const int HoleEffectModel = 5;


	int type;
	float sill;
	float angle[3];
	float range[3];
	float parameter;

	double **rotmat;

	double calculate(double hsqd, double cmax);
	void rotationMatrix();

public:
	~VariographicModelStructure(void);
};

class VariographicModel
{
public:
	VariographicModel(void);
	VariographicModel(float nuggetEffect,int size,float sill[],int vmtype[], float range[], double **rotmat,float *parameter);
	float nuggetEffect;
	VariographicModelStructure	**structures;
	int size;

	void computeCovariance(double * p1, double * p2, double *cmax, double *cova);
	void computeCovarianceModified(double * p1, double * p2, double *cmax, double *cova);
	double maxCovariance() ;
	void setup() ;

public:
	~VariographicModel(void);
};

//implementations
int split(char* line, char*** tokens){ //puntero a arreglo de *char
	std::vector<char *> vtokens;

	char* token = strtok (line," \t");
	while (token != NULL)
	{
		vtokens.push_back(token);
		token = strtok (NULL, " \t");
	}

	int numOfTokens= vtokens.size();

	*tokens = new char*[numOfTokens];

	for (int i=0;i<numOfTokens;i++) {
		(*tokens)[i] = vtokens[i];
	}
	return numOfTokens;
}

double **mallocPoints(int size) {
	double **points = new double*[size];

	for (int i = 0; i < size; i++) {
		points[i] = new double[3];
	}

	return points;
}

double sqdistance(const double p1[], const double p2[], double *rotmat[]) {
	//
	// Compute component distance vectors and the squared distance:
	//
	double sqdist = 0.0;
	double cont;

	/*
	 for (int i=0;i<3;i++) {
	 cont = rotmat[i][0] * (p1[0] - p2[0]) + rotmat[i][1] * (p1[1] - p2[1]) + rotmat[i][2] * (p1[2] - p2[2]);
	 sqdist += cont * cont;
	 }
	 */

	register double dif0 = p1[0] - p2[0];
	register double dif1 = p1[1] - p2[1];
	register double dif2 = p1[2] - p2[2];

	cont = rotmat[0][0] * dif0 + rotmat[0][1] * dif1 + rotmat[0][2] * dif2;
	sqdist += cont * cont;

	cont = rotmat[1][0] * dif0 + rotmat[1][1] * dif1 + rotmat[1][2] * dif2;
	sqdist += cont * cont;

	cont = rotmat[2][0] * dif0 + rotmat[2][1] * dif1 + rotmat[2][2] * dif2;
	sqdist += cont * cont;

	return sqdist;
}

void rotation_matrix(float ang1, float ang2, float ang3, float anis1,
		float anis2, double *rot[]) {
	double alpha;
	double beta;
	double theta;

	if (ang1 >= 0.0 && ang1 < 270.0) {
		alpha = (90.0 - ang1) * DEG2RAD;
	} else {
		alpha = (450.0 - ang1) * DEG2RAD;
	}
	beta = -ang2 * DEG2RAD;
	theta = ang3 * DEG2RAD;

	double sina = sin(alpha);
	double sinb = sin(beta);
	double sint = sin(theta);
	double cosa = cos(alpha);
	double cosb = cos(beta);
	double cost = cos(theta);

	double afac1 = 1.0 / MAX(anis1, EPSILON); //(anis1 < EPSLON ? EPSLON : anis1);
	double afac2 = 1.0 / MAX(anis2, EPSILON);

	rot[0][0] = cosb * cosa;
	rot[0][1] = cosb * sina;
	rot[0][2] = -sinb;

	rot[1][0] = afac1 * (-cost * sina + sint * sinb * cosa);
	rot[1][1] = afac1 * (cost * cosa + sint * sinb * sina);
	rot[1][2] = afac1 * (sint * cosb);

	rot[2][0] = afac2 * (sint * sina + cost * sinb * cosa);
	rot[2][1] = afac2 * (-sint * cosa + cost * sinb * sina);
	rot[2][2] = afac2 * (cost * cosb);

}

void createPoints(double ***ptr, int size) {
	*ptr = mallocPoints(size);
}

void freePoints(double **points, int size) {
	for (int i = 0; i < size; i++) {
		delete[] points[i];
	}

	delete[] points;
}

Grid::Grid()
{
}

Grid::~Grid()
{
}

int Grid::length()
{
	return nodes[0] * nodes[1] * nodes[2];
}

double *Grid::getPoint(int index) {
    int nx = this->nodes[0];
    int nxy = nx * this->nodes[1];

    int indexes[3];

    indexes[2] = (int)floor((double)(index / nxy));
    indexes[1] = (int)floor((double)((index - indexes[2] * nxy) / nx));
    indexes[0] = index - indexes[2] * nxy - indexes[1] * nx;


    double *p = new double[3];

    for (int j = 0; j < 3; j++) {
        p[j] = this->start[j] + indexes[j] * this->size[j];
    }

    return p;
}

double **Grid::getPoints() {
    int nx = nodes[0];
    int nxy = nx * nodes[1];
    int nxyz = nxy * nodes[2];
    double **p = new double*[3];

    for (int index=0;index<nxyz;index++) {
        p[index] = new double[3];
		p[index][2] = (float)floor((double)(index / nxy));
		p[index][1] = (float)floor((double)((index - p[index][2] * nxy) / nx));
		p[index][0] = index - p[index][2] * nxy - p[index][1] * nx;
		for (int j = 0; j < 3; j++) {
			p[index][j] = start[j] + p[index][j] * size[j];
		}
		//printf("Grilla indice %d: %s\n",index+1,p[index].toString().c_str());
    }
    return p;
}

VariographicModelStructure::VariographicModelStructure(void)
    :rotmat(NULL) {
	this->angle[0] = 0;
	this->angle[1] = 0;
	this->angle[2] = 0;
	this->range[0] = this->range[0];
	this->range[1] = this->range[1];
	this->range[2] = this->range[2];
	this->sill = 1.0;
	this->type = VariographicModelStructure::SphericalModel;
}

VariographicModelStructure::~VariographicModelStructure(void) {
    if(rotmat)
        freePoints(rotmat,3);

}

double VariographicModelStructure::calculate(double hsqd, double h) {
	double r = 0;
	double hr = h / range[0];

	switch (type) {
	case VariographicModelStructure::SphericalModel:
		if (hr < 1) {
			r = sill * (1.0 - hr * (1.5 - 0.5 * (hr * hr)));
		}
		break;
	case VariographicModelStructure::ExponentialModel:
		r = sill * exp(-3.0 * hr);
		break;
	default:
		break;
	}
	return r;
}

void VariographicModelStructure::rotationMatrix() {
	float max = range[0]; 
	float anis1 = (float) (range[1] / max);
	float anis2 = (float) (range[2] / max);

	createPoints(&rotmat,3);

	rotation_matrix(angle[0], angle[1], angle[2], anis1, anis2, rotmat);
}

VariographicModel::VariographicModel(void)
{
	structures = NULL;
	size = 0;
}

VariographicModel::~VariographicModel(void)
{
    for (int i = 0; i < size; i++) {
    	delete structures[i];
    }
    delete[] structures;
}

void VariographicModel::computeCovariance(double * p1, double * p2, double *cmax, double *cova)
{
    *cmax = nuggetEffect;

    int i;

    for (i = 0; i < size; i++) {
        *cmax += structures[i]->sill;
    }
    double hsqd = 0;
    double h = 0;

    if (size <= 0) {
        *cova = *cmax;
    	return;
    }
    *cova = 0.0;

    for (i = 0; i < size; i++) {
        hsqd = sqdistance(p1,p2,structures[i]->rotmat);
        if (i == 0 && hsqd < EPSILON) {
            *cova = *cmax;
            return;
        }
        h = sqrt(hsqd);
	double cov = structures[i]->calculate(hsqd, h);
	*cova += cov;
    }
}


//void VariographicModel::computeCovarianceModified(double * p1, double * p2, double *cmax, double *cova)
//{
//    int i;
//    double hsqd = 0;
//    double h = 0;
//    if (size <= 0) {
//        *cova = *cmax;
//    	return;
//    }
//    *cova = 0.0;
//    for (i = 0; i < size; i++) {
//        hsqd = sqdistance(p1,p2,structures[i]->rotmat);
//        if (i == 0 && hsqd < EPSILON) {
//            *cova = *cmax;
//            return;
//        }
//        h = sqrt(hsqd);
//	double cov = structures[i]->calculate(hsqd, h);
//	*cova += cov;
//    }
//}




double VariographicModel::maxCovariance()
{
    double cmax = 0;
    cmax = nuggetEffect;
    for (int i = 0; i < size; i++) {
		cmax += structures[i]->sill;
    }
    return cmax;
}

void VariographicModel::setup()
{
    for (int i = 0; i < size; i++) {
		this->structures[i]->rotationMatrix();
    }
}

int load(const char *fileName, double ***loadPoints,double **loadValues, 
		int xColumn, int yColumn, int zColumn, int valueColumn, 
		double trimmingStart, double trimmingEnd,
		int *rows)
{
	printf("filename = %s\n",fileName);
	ifstream infile;

	infile.open(fileName, ifstream::in);

	if (!infile)
	{
		cerr << " __FILE__ load No se puede abrir archivo." << endl;
		return 1;
	}

	string line;

	//se cuenta el largo del archivo
	int length = 0;

	//leo nombre de las variables
	while (getline(infile, line))
	{
		length++;
	}

	//reabro el archivo para leer los datos
	infile.close();
	infile.open(fileName, ifstream::in);

	int numOfVars = 0;

	//leo el titulo
	if (! getline(infile, line))
	{
		cerr << "Formato Invalido" << endl;
		return 0;
	}

	//leo el numero de variables
	if (! getline(infile, line))
	{
		cerr << "Formato Invalido" << endl;
		return 0;
	}
	sscanf(line.c_str(), "%u", &numOfVars);

	printf("Variables = %d\n",numOfVars);

	//leo nombre de las variables
	int i;
	for (i=0;i<numOfVars;i++) {
		getline(infile, line);
	}

	if(i != numOfVars)
	{
		cerr << "Formato Invalido. Numero de variables inconsistente" << endl;
		return 0;
	}

	*rows = length-(numOfVars+2);
	//printf("rows = %d\n",*rows);

	createPoints(loadPoints,*rows);

	//*loadPoints = (Point *) malloc( (*rows) * sizeof(Point));
	if (loadValues != NULL) {
		*loadValues =  new double[*rows]; //(double *) malloc( (*rows) *sizeof(double ));
	}
	//printf("Puntos y datos allocados\n");

	double x, y, z, v;

	i=0;
	char** tokens;
	int numOfTokens;
	while ( getline(infile, line) )
	{
		char *cline=new char[line.length()+1];
		strcpy(cline, line.c_str());
		//printf("linea [%s]\n",cline);
		numOfTokens = split(cline, &tokens);
		//printf("linea %d, tokens %d\n",i+1,numOfTokens);

		/*
		for (int kk=0;kk<numOfTokens;kk++) {
			printf("token %d = [%s]\n",kk,tokens[kk]);
		}
		*/
		if(xColumn > -1){
			x = atof(tokens[xColumn]);
		}
		if(yColumn > -1){
			y = atof(tokens[yColumn]);
		}
		if(zColumn > -1){
			z = atof(tokens[zColumn]);
		}
		if(valueColumn > -1){
			v = atof(tokens[valueColumn]);
		}

		//printf("[%f;%f;%f] [%f,%f]\n", x,y,z,v,w);

		//ojo con esto, si no se ponen los parentesis no sirve!
		(*loadPoints)[i][0] = x;
		(*loadPoints)[i][1] = y;
		(*loadPoints)[i][2] = z;

		//printf("Punto asignado\n");
		if (loadValues != NULL && valueColumn >= 0) {
			if (v >= trimmingStart && v <= trimmingEnd) {
				(*loadValues)[i] = v;
				i++;
			}
		} else {
			i++;
		}
		delete[] cline;
		//delete [] tokens;

	}

	infile.close();

	(*rows) = i;

	return 1;
}

//double recursiveComputeCovariance(int i, int j, int k);


Grid grid;
VariographicModel vm;
VariographicModelStructure *st;
double **inputPoints;
double *inputValues;
int inputSize;
int gridSize;
double cmax;
//double *lastColumnNew;
//double *lastColumnOld;

//-----------------------------------------------------------------------------------------

int howmany;
int whoAmI;

//#define N 2000   // 1500
//#define N 4800   // 1500
//#define N 31176
//#define N 117576

//this are ok
//#define N 9527
#define N 9576
//#define N 31127
//#define N 117527


//#define N 9600   // 1500
//#define N 36000   // 1500
//#define N 9600   // 1500
//#define N 9630   // 1500
//#define N 4830   // 1500
//#define B 500
//#define B 1000
#define B 32
//#define B 512

#define SIMS 1
//#define SIMS 5000
//#define SIMS 10000
//#define SIMS 20000
//#define SIMS 40000
//#define SIMS 80000
//#define SIMS 160000
//#define SIMS 320000

//#define N 14   // 1500
//#define B 2

#define MIN(x,y)  ((x) < (y) ? (x) : (y))

/* double A[N][N]; */
double **A;
double **Atmp;
double **simval;
double *randval;
//int N;

//--------------------------------------------------------------
void init_geostat_stuff(char *argv[]){
/*
 *
 *
Grilla 1 grande:

Nx,ny,nz: 20 x 30 x 12
Xmn,ymn,zmn: -180, -280, 6
Xsiz,ysiz,zsiz: 40 x 40 x 12

 Grilla 2 grande:

Nx,ny,nz: 40 x 60 x 12
Xmn,ymn,zmn: -190, -290, 6
Xsiz,ysiz,zsiz: 20 x 20 x 12

 Grilla 3 grande:

Nx,ny,nz: 80 x 120 x 12
Xmn,ymn,zmn: -195, -295, 6
Xsiz,ysiz,zsiz: 10 x 10 x 12

 Grilla 1chica:

Nx,ny,nz: 20 x 30 x 12
Xmn,ymn,zmn: 152.5, 227.5, 6
Xsiz,ysiz,zsiz: 5 x 5 x 12

 Grilla 2 chica:

Nx,ny,nz: 40 x 60 x 12
Xmn,ymn,zmn: 151.25, 226.25, 6
Xsiz,ysiz,zsiz: 2.5 x 2.5 x 12

 Grilla 3 chica:

Nx,ny,nz: 80 x 120 x 12
Xmn,ymn,zmn: 150.625, 225.625, 6
Xsiz,ysiz,zsiz: 1.25 x 1.25 x 12
 *
 *
 */


	//Grid grid = Grid();
	grid = Grid();

	// Grilla 1 grande
	grid.nodes[0] = 20; 	// nx
	grid.nodes[1] = 30; 	// ny
	grid.nodes[2] = 12; 	// nz
	grid.start[0] = -180.0;	// xmn
	grid.start[1] = -280.0;	// ymn
	grid.start[2] = 6;	// zmn
	grid.size[0] = 40.0;	// xsiz
	grid.size[1] = 40.0;	// ysiz
	grid.size[2] = 12;	// zsiz
	////gridSize=7200, inputSize=2376, total=9576
	////gridSize=7200, inputSize=2327, total=9527  <- deleted_rows

	// Grilla 2 grande
	//grid.nodes[0] = 40; 	// nx
	//grid.nodes[1] = 60; 	// ny
	//grid.nodes[2] = 12; 	// nz
	//grid.start[0] = -190.0;	// xmn
	//grid.start[1] = -290.0;	// ymn
	//grid.start[2] = 6;	// zmn
	//grid.size[0] = 20.0;	// xsiz
	//grid.size[1] = 20.0;	// ysiz
	//grid.size[2] = 12;	// zsiz
	////gridSize=28800, inputSize=2376, total=31176

//	// Grilla 3 grande
//	grid.nodes[0] = 80; 	// nx
//	grid.nodes[1] = 120; 	// ny
//	grid.nodes[2] = 12; 	// nz
//	grid.start[0] = -195.0;	// xmn
//	grid.start[1] = -295.0;	// ymn
//	grid.start[2] = 6;	// zmn
//	grid.size[0] = 10.0;	// xsiz
//	grid.size[1] = 10.0;	// ysiz
//	grid.size[2] = 12;	// zsiz
//	//gridSize=115200, inputSize=2376, total=117576
//	//gridSize=115200, inputSize=2327, total=117527  nscore

	// Grilla 2 chica
	//grid.nodes[0] = 40; 	// nx
	//grid.nodes[1] = 60; 	// ny
	//grid.nodes[2] = 12; 	// nz
	//grid.start[0] = 151.25;	// xmn
	//grid.start[1] = 226.25;	// ymn
	//grid.start[2] = 6;	// zmn
	//grid.size[0] = 2.5;	// xsiz
	//grid.size[1] = 2.5;	// ysiz
	//grid.size[2] = 12;	// zsiz
	////gridSize=28800, inputSize=2376, total=31176
	//gridSize=28800, inputSize=2327, total=31127  nscore


	// Grilla 3 chica
	//grid.nodes[0] = 80; 	// nx
	//grid.nodes[1] = 120; 	// ny
	//grid.nodes[2] = 12; 	// nz
	//grid.start[0] = 150.625;	// xmn
	//grid.start[1] = 225.625;	// ymn
	//grid.start[2] = 6;	// zmn
	//grid.size[0] = 1.25;	// xsiz
	//grid.size[1] = 1.25;	// ysiz
	//grid.size[2] = 12;	// zsiz
	////gridSize=115200, inputSize=2376, total=117576

	//VariographicModel vm = VariographicModel();
	vm = VariographicModel();
	vm.nuggetEffect = 0.1;
	st = new VariographicModelStructure();
	st->angle[0] = 0;
	st->angle[1] = 0;
	st->angle[2] = 0;
	st->range[0] = 100;
	st->range[1] = 100;
	st->range[2] = 100;
	st->sill = 0.9;
	st->type = VariographicModelStructure::SphericalModel;
	vm.size = 1;
	vm.structures = new VariographicModelStructure*[1];
	vm.structures[0] = st;
	vm.setup();
	//0.1 Pepa + 0.9 Esférico (100,100,100)

	//double **inputPoints;
	//double *inputValues;

	//int gridSize = grid.length();
	gridSize = grid.length();
	//int inputSize;
	char *inputFile = argv[1] ;//"/home/exequiel/workspace/PGSL-lib/data/muestras.dat";
	//char *outputFile = argv[2] ;//"/home/exequiel/Projects/alges/oscar_peredo/covamatrix-dense.txt";

	load(inputFile,&inputPoints,&inputValues,0,1,2,3,0,1e100, &inputSize);

	//FILE *fout = fopen(outputFile,"w");
	//FILE *fout = stdout;


	printf("gridSize=%d, inputSize=%d, total=%d\n",gridSize,inputSize,gridSize+inputSize);

	//double * pi;
	//double * pj;
	//double cmax;
	//double cova;
    	cmax = vm.nuggetEffect;
    	for (int i = 0; i < vm.size; i++) {
        	cmax += vm.structures[i]->sill;
    	}

}

void delete_geostat_stuff(){
	freePoints(inputPoints,inputSize);
    	delete [] inputValues;
}
//--------------------------------------------------------------
// normal random generators

/******************************************************************************/
//	Standard version with trigonometric calls
//#define PI 3.14159265358979323846

/******************************************************************************/
//	"Polar" version without trigonometric calls

//double randn_notrig(double mu, double sigma);
//double randn_trig(double mu, double sigma);

/*
int main(int argc, char** argv){
	double a,b;
	int i=0;
	
	srand(time(0));
	
	while(i<atoi(argv[1])){
		a=randn_notrig(0.0,1.0);
		b=randn_trig(0.0,1.0);
		printf("%f %f\n",a,b);
		i++;
	}
}
*/

//double randn_notrig(double mu, double sigma) {
double randn_notrig() {
	static int deviateAvailable=0;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double polar, rsquared, var1, var2;
	

	//	If no deviate has been stored, the polar Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose pairs of uniformly distributed deviates, discarding those 
		//	that don't fall within the unit circle
		do {
			var1=2.0*( (double)(rand())/(double)(RAND_MAX) ) - 1.0;
			var2=2.0*( (double)(rand())/(double)(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);
		
		//	calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);
		
		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=1;
		
		//	return second deviate
		//return var2*polar*sigma + mu;
		return var2*polar;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=0;
		//return storedDeviate*sigma + mu;
		return storedDeviate;
	}
}



//double randn_trig(double mu, double sigma) {
double randn_trig() { // mu=0.0, sigma=1.0
	static int deviateAvailable=0;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double dist, angle;
	

	//	If no deviate has been stored, the standard Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	//if (!deviateAvailable) {
		
		//	choose a pair of uniformly distributed deviates, one for the
		//	distance and one for the angle, and perform transformations
		dist=sqrt( -2.0 * log( (double)(rand()) / (double)(RAND_MAX) ));
		angle=2.0 * PI * ((double)(rand()) / (double)(RAND_MAX));
		
		//	calculate and store first deviate and set flag
		storedDeviate=dist*cos(angle);
		deviateAvailable=1;
		
		//	calcaulate return second deviate
		//return dist * sin(angle) * sigma + mu;
		return dist * sin(angle) ;
	//}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	//else {
	//	deviateAvailable=0;
	//	//return storedDeviate*sigma + mu;
	//	return storedDeviate;
	//}
}










//--------------------------------------------------------------
void genmat_seq ()
{
   int init_val, i, j, diag_dist;

   init_val = 1325;

   for (i = 0; i < N; i++) 
      for (j = 0; j < N; j++)
   	{
	init_val = (3125 * init_val) % 65536;
        diag_dist = (i<j) ? j-i: i-j;
        A[i][j] = 1000.0/(diag_dist+1)+ (init_val - 32768.0) / 16384.0;
   	}
}

void printres_seq (int volcar)
{
   int i,j;

   if (volcar == 1) {
      for (i = 0; i < N; i++) {
         for (j = 0; j < N; j++)
           printf("%f ", A[i][j]);
         printf ("\n");
	 }
   } else
      printf("A[%d][%d]=%f ", N/2, N/2, A[N/2][N/2]);
}



void genmat ()
{
   int init_val, i, j,kk, diag_dist,init_val_local,iinf,isup,NB;
   double inv;
   MPI_Status st;
   
   NB = (N+B-1)/B;
   printf("%d: total blocks %d\n",whoAmI,NB);

   inv=1.0/16384.0;

   for(kk=0;kk<NB;kk++){
      if(kk%howmany == whoAmI){
         printf("%d: generating block %d\n",whoAmI,kk);
         iinf = kk*B;
         isup = MIN((kk+1)*B, N);
         /*set or receive the value of init_val*/
         if(kk==0){
            init_val = 1325;
         }
         else{
            if(whoAmI==0)
               MPI_Recv(&init_val,1,MPI_INT,howmany-1,0,MPI_COMM_WORLD,&st);
            else
               MPI_Recv(&init_val,1,MPI_INT,whoAmI-1,0,MPI_COMM_WORLD,&st);
         }

         /*initialization of A*/
//#pragma omp parallel shared(init_val)
//{
//#pragma omp for schedule(runtime) private(i,diag_dist) nowait
         for (i = iinf; i < isup; i++){
            for (j = 0; j < N; j++){
	       init_val = (3125 * init_val) % 65536;
               diag_dist = (i<j) ? j-i: i-j;
               A[i%B+(kk/howmany)*B][j] = 1000.0/(diag_dist+1)+ (init_val - 32768.0) * inv;
   	    }
         }
//}

         /*send the value of init_val*/
         if(kk<NB-1){
            if(whoAmI==howmany-1)   
               MPI_Send(&init_val,1,MPI_INT,0,0,MPI_COMM_WORLD);
            else
               MPI_Send(&init_val,1,MPI_INT,whoAmI+1,0,MPI_COMM_WORLD);
         }

      }
   }



}


//-----------------------------------------------------

void setN(){
//   N= inputSize+gridSize;
   printf("N=%d\tinputSize=%d\tgridSize=%d\n",N,inputSize,gridSize);
}


void genrandval(){
   int i, j,kk,iinf,isup,NB;
   double inv, cova;
   MPI_Status st;

   NB = (N+B-1)/B;

   srand((unsigned)((time(NULL)*whoAmI) % INT_MAX));

   for(kk=0;kk<NB;kk++){
      if(kk%howmany == whoAmI){
         iinf = kk*B;
         isup = MIN((kk+1)*B, N);
         for (i = iinf; i < isup; i++){
            for (j = 0; j < N; j++){
               randval[i%B+(kk/howmany)*B] = ((double) rand() / (RAND_MAX));;
   	    }
         }
      }
   }	
	 
}


void genmat_geostat_stuff()
{
   int init_val, i, j,kk, diag_dist,init_val_local,iinf,isup,NB;
   double inv, cova;
   MPI_Status st;

   NB = (N+B-1)/B;

   //inv=1.0/16384.0;

   double *pi;
   double *pj;




   for(kk=0;kk<NB;kk++){
      if(kk%howmany == whoAmI){
         iinf = kk*B;
         isup = MIN((kk+1)*B, N);
         for (i = iinf; i < isup; i++){
            if (i<gridSize) {
               pi = grid.getPoint(i);
            } else {
               pi = inputPoints[i-gridSize];
            }

            for (j = 0; j < N; j++){
               if (j<gridSize) {
                  pj = grid.getPoint(j);
               } else {
                  pj = inputPoints[j-gridSize];
               }

               vm.computeCovariance(pi,pj,&cmax,&cova);
               A[i%B+(kk/howmany)*B][j] = cova;

               if (j<gridSize) {
                  delete pj;
               }
   	    }
            if (i<gridSize) {
               delete pi;
            }
         }
      }
   }
}
//-----------------------------------------------------


void printres (int volcar)
{
   int i,j,kk,NB,iinf,isup,index;

   NB = (N+B-1)/B;

   if (volcar == 1) {
      for(kk=0;kk<NB;kk++){
         if(kk%howmany == whoAmI){
            iinf = kk*B;
            isup = MIN((kk+1)*B, N);
            for (i = iinf; i < isup; i++) {
               for (j = 0; j < N; j++)
                  printf("%f ", A[i%B+(kk/howmany)*B][j]);
               printf ("  %d %d \n",i,whoAmI);
	    }
         }
         sleep(2);
         MPI_Barrier(MPI_COMM_WORLD);
      }
   }
   else{
      kk=(N/2)/B;
      if(kk % howmany == whoAmI){
         index=(N/2)%B+(kk/howmany)*B;  
         printf("A[%d][%d]=%f (thread %d)\n", N/2, N/2, A[index][N/2],whoAmI);
      }
   }
}

void lu0(int kk)
{
	int i, j, k;
	int kinf, ksup;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);

        for (k=kinf; k<ksup; k++)
           for (i=k+1; i<ksup; i++) {
              if(A[k][k]==0.0){
                 printf("divide by cero!!!!!!\n");
                 exit (-1);
              }
              A[i][k] = A[i][k] / A[k][k];
              for (j=k+1; j<ksup; j++)
            	 A[i][j] = A[i][j] - A[i][k] * A[k][j];
	      }
}


void lu0local(int kk)
{
	int i, j, k, ilocal, klocal;
	int kinf, ksup;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);

        for (k=kinf; k<ksup; k++){
           klocal= k%B+(kk/howmany)*B;
           for (i=k+1; i<ksup; i++) {
              ilocal= i%B+(kk/howmany)*B;
              if(A[klocal][k]==0.0){
                 printf("divide by cero!!!!!!\n");
                 exit (-1);
              }
              A[ilocal][k] = A[ilocal][k] / A[klocal][k];
              for (j=k+1; j<ksup; j++)
            	 A[ilocal][j] = A[ilocal][j] - A[ilocal][k] * A[klocal][j];
	      }
        }
}

void lu0localChol(int kk)
{
	int i, j, k, ilocal, klocal, p;
	int kinf, ksup;

	double A_klocal_k, A_ilocal_k, randval_klocal, A_klocal_k_inv;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);

        for (k=kinf; k<ksup; k++){
           klocal= k%B+(kk/howmany)*B;
           
           A[klocal][k] =  sqrt(A[klocal][k]);
           A_klocal_k = A[klocal][k];
           A_klocal_k_inv = 1.0/A_klocal_k;

           randval_klocal= randval[klocal]; 

           for(int ss=0;ss<SIMS;ss++)
           	simval[klocal][ss]+= A_klocal_k * randn_notrig();
           	//simval[klocal][ss]+= A_klocal_k * randval_klocal;

           printf("%d: klocal=%d simval[%d][0]=%f\n",whoAmI,klocal,k,simval[klocal][0]);

//           for (i=k+1; i<ksup; i++) {
//              ilocal= i%B+(kk/howmany)*B;
//              //A[ilocal][k] = A[ilocal][k] / A[klocal][k];
//              //A[ilocal][k] = A[ilocal][k] / A_klocal_k;
//              A[ilocal][k] = A[ilocal][k] * A_klocal_k_inv;
//              A_ilocal_k = A[ilocal][k];  
//              for (j=k+1; j<ksup; j++)
//            	 //A[ilocal][j] = A[ilocal][j] - A[ilocal][k] * A[klocal][j];
//            	 A[ilocal][j] = A[ilocal][j] - A_ilocal_k * A[klocal][j];
//
//              //simval[ilocal]+= A[ilocal][k] * randval[klocal];
//              for(int ss=0;ss<SIMS;ss++)
//                 simval[ilocal][ss]+= A_ilocal_k * randval_klocal;
//	      }

           for (i=k+1; i<ksup; i++) {
              ilocal= i%B+(kk/howmany)*B;
              //A[ilocal][k] = A[ilocal][k] / A[klocal][k];
              //A[ilocal][k] = A[ilocal][k] / A_klocal_k;
              A[ilocal][k] = A[ilocal][k] * A_klocal_k_inv;
              A_ilocal_k = A[ilocal][k];  
              for (j=k+1; j<ksup; j++)
            	 //A[ilocal][j] = A[ilocal][j] - A[ilocal][k] * A[klocal][j];
            	 A[ilocal][j] = A[ilocal][j] - A_ilocal_k * A[klocal][j];

              //simval[ilocal]+= A[ilocal][k] * randval[klocal];
              for(int ss=0;ss<SIMS;ss++)
                 simval[ilocal][ss]+= A_ilocal_k * randn_notrig();
                 //simval[ilocal][ss]+= A_ilocal_k * randval_klocal;




//              ilocal= (i+1)%B+(kk/howmany)*B;
//              //A[ilocal][k] = A[ilocal][k] / A[klocal][k];
//              //A[ilocal][k] = A[ilocal][k] / A_klocal_k;
//              A[ilocal][k] = A[ilocal][k] * A_klocal_k_inv;
//              A_ilocal_k = A[ilocal][k];  
//              for (j=k+1; j<ksup; j++)
//            	 //A[ilocal][j] = A[ilocal][j] - A[ilocal][k] * A[klocal][j];
//            	 A[ilocal][j] = A[ilocal][j] - A_ilocal_k * A[klocal][j];
//
//              //simval[ilocal]+= A[ilocal][k] * randval[klocal];
//              for(int ss=0;ss<SIMS;ss++)
//                 simval[ilocal][ss]+= A_ilocal_k * randval_klocal;

	   }



//           for (; i<ksup; i++) {
//              ilocal= i%B+(kk/howmany)*B;
//              //A[ilocal][k] = A[ilocal][k] / A[klocal][k];
//              //A[ilocal][k] = A[ilocal][k] / A_klocal_k;
//              A[ilocal][k] = A[ilocal][k] * A_klocal_k_inv;
//              A_ilocal_k = A[ilocal][k];  
//              for (j=k+1; j<ksup; j++)
//            	 //A[ilocal][j] = A[ilocal][j] - A[ilocal][k] * A[klocal][j];
//            	 A[ilocal][j] = A[ilocal][j] - A_ilocal_k * A[klocal][j];
//
//              //simval[ilocal]+= A[ilocal][k] * randval[klocal];
//              for(int ss=0;ss<SIMS;ss++)
//                 simval[ilocal][ss]+= A_ilocal_k * randval_klocal;
//	   }




        }
}

void bdiv(int ii, int kk) {
	int i, j, k;
	int kinf, ksup, iinf, isup;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        iinf = ii*B; 
        isup = MIN((ii+1)*B, N);

	for (k=kinf; k<ksup; k++) {
	   for (i=iinf; i<isup; i++) {
              A[i][k] = A[i][k] / A[k][k];
              for (j=k+1; j<ksup; j++)
                 A[i][j] = A[i][j] - A[i][k]*A[k][j];
              }
	   }
}

void bdivlocal(int ii, int kk) {
	int i, j, k,ilocal,klocal;
	int kinf, ksup, iinf, isup;
        double tmp;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        iinf = ii*B; 
        isup = MIN((ii+1)*B, N);

        if(kk%howmany != whoAmI){
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
	   for (i=iinf; i<isup; i++) {
              ilocal= i%B+(ii/howmany)*B;
	      for (k=kinf; k<ksup; k++) {
                 A[ilocal][k] = A[ilocal][k] / Atmp[k-kinf][k];
                 tmp=A[ilocal][k]; 
                 for (j=k+1; j<ksup; j++)
                    A[ilocal][j] = A[ilocal][j] - tmp*Atmp[k-kinf][j];
              }
	   }
}
        }
        else{
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
	   for (i=iinf; i<isup; i++) {
              ilocal= i%B+(ii/howmany)*B;
	      for (k=kinf; k<ksup; k++) {
                 klocal= k%B+(kk/howmany)*B;
                 A[ilocal][k] = A[ilocal][k] / A[klocal][k];
                 tmp=A[ilocal][k]; 
                 for (j=k+1; j<ksup; j++)
                    A[ilocal][j] = A[ilocal][j] - tmp*A[klocal][j];
              }
	   }
}
        }
}


void bdivlocalChol(int ii, int kk) {
	int i, j, k,ilocal,klocal, kinflocal;
	int kinf, ksup, iinf, isup;
        double tmp, A_ilocal_k, randval_klocal;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        iinf = ii*B; 
        isup = MIN((ii+1)*B, N);

        if(kk%howmany != whoAmI){
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
	   for (i=iinf; i<isup; i++) {
              ilocal= i%B+(ii/howmany)*B;
	      for (k=kinf; k<ksup; k++) {
                 kinflocal=k-kinf;
                 A[ilocal][k] = A[ilocal][k] / Atmp[kinflocal][k];
                 A_ilocal_k=A[ilocal][k];

                 for (j=k+1; j<ksup; j++)
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*Atmp[kinflocal][j];

                 for(int ss=0;ss<SIMS;ss++)                 
                    simval[ilocal][ss]+= A_ilocal_k * randn_notrig();
                    //simval[ilocal][ss]+= A_ilocal_k * randval_klocal;
              }
	   }
}
        }
        else{
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
	   for (i=iinf; i<isup; i++) {
              ilocal= i%B+(ii/howmany)*B;
	      for (k=kinf; k<ksup; k++) {
                 klocal= k%B+(kk/howmany)*B;
                 A[ilocal][k] = A[ilocal][k] / A[klocal][k];
                 A_ilocal_k=A[ilocal][k];

                 //randval_klocal =randval[klocal];  
                 for (j=k+1; j<ksup; j++)
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*A[klocal][j];

                 for(int ss=0;ss<SIMS;ss++)                 
                    simval[ilocal][ss]+= A_ilocal_k * randn_notrig();
                    //simval[ilocal][ss]+= A_ilocal_k * randval_klocal;
              }
	   }
}
        }
}


void bmod(int ii, int jj, int kk) {
        int i, j, k;
        int kinf, ksup, iinf, isup, jinf, jsup;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        iinf = ii*B;
        isup = MIN((ii+1)*B, N);

        for (k=kinf; k<ksup; k++) 
           for (i=iinf; i<isup; i++)
              for (j=jinf; j<jsup; j++)
                 A[i][j] = A[i][j] - A[i][k]*A[k][j];
}

void bmodlocal(int ii, int jj, int kk) {
        int i, j, k,ilocal,klocal;
        int kinf, ksup, iinf, isup, jinf, jsup;
        double tmp;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        iinf = ii*B;
        isup = MIN((ii+1)*B, N);

        if(kk%howmany != whoAmI){
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
           for (i=iinf; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              for (k=kinf; k<ksup; k++){
                 tmp=A[ilocal][k]; 
                 for (j=jinf; j<jsup; j++)
                    A[ilocal][j] = A[ilocal][j] - tmp*Atmp[k-kinf][j];
              }
           }
}
        }
        else{
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,klocal,tmp) nowait
           for (i=iinf; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              for (k=kinf; k<ksup; k++){
                 klocal= k%B+(kk/howmany)*B;
                 tmp=A[ilocal][k]; 
                 for (j=jinf; j<jsup; j++)
                    A[ilocal][j] = A[ilocal][j] - tmp*A[klocal][j];
              }
           }
}
        }
}


void bmodlocalChol(int ii, int jj, int kk) {
        int i, j, k,ilocal,klocal,kinflocal,kinflocal2;
        int kinf, ksup, iinf, isup, jinf, jsup;
        double tmp, A_ilocal_k, A_ilocal_k2;

        kinf = kk*B;
        ksup = MIN((kk+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        iinf = ii*B;
        isup = MIN((ii+1)*B, N);

        if(kk%howmany != whoAmI){
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,tmp) nowait
           for (i=iinf; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              for (k=kinf; k<ksup; k++){
                 kinflocal=k-kinf;
                 A_ilocal_k=A[ilocal][k]; 
                 for (j=jinf; j<jsup; j++)
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*Atmp[kinflocal][j];
              }

/*
	      for (k=kinf; k<ksup-1; k=k+2){
                 kinflocal=k-kinf;
                 A_ilocal_k=A[ilocal][k]; 
                 kinflocal2=k+1-kinf;
                 A_ilocal_k2=A[ilocal][k+1]; 

                 for (j=jinf; j<jsup-3; j=j+4){
                    A[ilocal][j]   = A[ilocal][j]   - A_ilocal_k*Atmp[kinflocal][j];
                    A[ilocal][j+1] = A[ilocal][j+1] - A_ilocal_k*Atmp[kinflocal][j+1];
                    A[ilocal][j+2] = A[ilocal][j+2] - A_ilocal_k*Atmp[kinflocal][j+2];
                    A[ilocal][j+3] = A[ilocal][j+3] - A_ilocal_k*Atmp[kinflocal][j+3];

		    A[ilocal][j]   = A[ilocal][j]   - A_ilocal_k2*Atmp[kinflocal2][j];
                    A[ilocal][j+1] = A[ilocal][j+1] - A_ilocal_k2*Atmp[kinflocal2][j+1];
                    A[ilocal][j+2] = A[ilocal][j+2] - A_ilocal_k2*Atmp[kinflocal2][j+2];
                    A[ilocal][j+3] = A[ilocal][j+3] - A_ilocal_k2*Atmp[kinflocal2][j+3];

		 }
                 for (; j<jsup; j++){
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*Atmp[kinflocal][j];
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k2*Atmp[kinflocal2][j];
                 }


              }

	      for (; k<ksup; k++){
                 kinflocal=k-kinf;
                 A_ilocal_k=A[ilocal][k]; 
                 for (j=jinf; j<jsup; j++)
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*Atmp[kinflocal][j];
              }
*/

           }
}
        }
        else{
#pragma omp parallel
{
#pragma omp for schedule(runtime) private(i,ilocal,klocal,tmp) nowait
           for (i=iinf; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              for (k=kinf; k<ksup; k++){
                 klocal= k%B+(kk/howmany)*B;
                 A_ilocal_k=A[ilocal][k]; 
                 for (j=jinf; j<jsup; j++)
                    A[ilocal][j] = A[ilocal][j] - A_ilocal_k*A[klocal][j];
              }
           }
}
        }
}




void fwd(int ii, int jj) {
        int i, j, k;
        int kinf, ksup, iinf, isup, jinf, jsup;

        kinf = ii*B;
        ksup = MIN((ii+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        isup = MIN((ii+1)*B, N);

        for (k=kinf; k<ksup; k++) 
           for (i=k+1; i<isup; i++)
              for (j=jinf; j<jsup; j++)
                 A[i][j] = A[i][j] - A[i][k]*A[k][j];
}

void fwdlocal(int ii, int jj) {
        int i, j, k, ilocal,klocal;
        int kinf, ksup, iinf, isup, jinf, jsup;

        kinf = ii*B;
        ksup = MIN((ii+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        isup = MIN((ii+1)*B, N);

        for (k=kinf; k<ksup; k++){ 
           klocal= k%B+(ii/howmany)*B;
           for (i=k+1; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              for (j=jinf; j<jsup; j++)
                 A[ilocal][j] = A[ilocal][j] - A[ilocal][k]*A[klocal][j];
           }
        }
}


void fwdlocalChol(int ii, int jj) {
        int i, j, k, ilocal,klocal;
        int kinf, ksup, iinf, isup, jinf, jsup;

	double A_ilocal_k;

        kinf = ii*B;
        ksup = MIN((ii+1)*B, N);
        jinf = jj*B;
        jsup = MIN((jj+1)*B, N);
        isup = MIN((ii+1)*B, N);

        for (k=kinf; k<ksup; k++){ 
           klocal= k%B+(ii/howmany)*B;
           for (i=k+1; i<isup; i++){
              ilocal= i%B+(ii/howmany)*B;
              A_ilocal_k = A[ilocal][k];
              for (j=jinf; j<jsup; j++)
                 A[ilocal][j] = A[ilocal][j] - A_ilocal_k*A[klocal][j];
           }
        }
}




long usecs (void)
{
  struct timeval t;

  gettimeofday(&t,NULL);
  return t.tv_sec*1000000+t.tv_usec;
}

void sendBlocks(int kk,int receiverId, MPI_Request* req){
      int iinf, isup, offset,i,tag;
      iinf = kk*B;
      isup = MIN((kk+1)*B, N);
         for (i = iinf; i < isup; i++) {
            offset= i%B+(kk/howmany)*B;
            tag=i;
            MPI_Isend(&A[offset][0],N,MPI_DOUBLE,receiverId,i,MPI_COMM_WORLD,&req[i-iinf]);
         }
}


void recvBlocks(int kk,int senderId,MPI_Request* req){
         MPI_Status st;
         int i,tag;
         int NB;
   
         NB = (N+B-1)/B;
	 
         if(kk<NB-1){
            for (i = 0; i < B; i++) {
               tag=kk*B+i;
               MPI_Irecv(&Atmp[i][0],N,MPI_DOUBLE,senderId,tag,MPI_COMM_WORLD,&req[i]);
               MPI_Wait(&req[i],&st);         
            }
         }
         else{
            for (i = 0; i < N%B ; i++) {
               tag=kk*B+i;
               MPI_Irecv(&Atmp[i][0],N,MPI_DOUBLE,senderId,tag,MPI_COMM_WORLD,&req[i]);
               MPI_Wait(&req[i],&st);         
            }
         }
}

int main(int argc, char* argv[])
{
   long t_start,t_end;
   double time;
   int NB;
   int ii, jj, kk,nrows,mod,id,i,j,count,offset,p;
   MPI_Status st,st2;
   MPI_Request req[48][B],req2[48][B];
      
   printf("%d: llego aca N=%d, B=%d\n",whoAmI,N,B);

   MPI_Init(&argc,&argv);

   MPI_Comm_size(MPI_COMM_WORLD, &howmany);
   MPI_Comm_rank(MPI_COMM_WORLD, &whoAmI);


   //-----------------------------
   init_geostat_stuff(argv);
   //setN();
   //-----------------------------


   
   NB = (N+B-1)/B;
   nrows= NB / howmany; //118 para N=117000
   mod= NB % howmany;
   
   
   for(id=0;id<mod;id++){
      if(whoAmI==id)
         nrows++;
   }

   //printf("%d: N=%d, B=%d, NB=%d, N%B=%d, nrows=%d, A(%d*%d,%d)...\n",whoAmI,N,B,NB,(N%B),nrows,nrows,B,N);
   //printf("%d: N=%d, B=%d, NB=%d, N%B=%d, nrows=%d, Atmp(%d,%d)...\n",whoAmI,N,B,NB,(N%B),nrows,B,N);
   printf("%d: N=%d, B=%d, NB=%d, nrows=%d, A(%d*%d,%d)...\n",whoAmI,N,B,NB,nrows,nrows,B,N);
   printf("%d: N=%d, B=%d, NB=%d, nrows=%d, Atmp(%d,%d)...\n",whoAmI,N,B,NB,nrows,B,N);


      printf("%d: allocating A(%d*%d,%d)...\n",whoAmI,nrows,B,N);

      //simval = (double **)calloc(SIMS,sizeof(double));
      //for (ii=0; ii<SIMS; ii++)
      //   simval[ii] = (double *)calloc(nrows*B,sizeof(double));
      simval = (double **)calloc(nrows*B,sizeof(double));
      for (ii=0; ii<nrows*B; ii++)
         simval[ii] = (double *)calloc(SIMS,sizeof(double));

      randval = (double *)calloc(nrows*B,sizeof(double));

      A = (double **)malloc(nrows*B*sizeof(double *));
      for (ii=0; ii<nrows*B; ii++)
         A[ii] = (double *)malloc(N*sizeof(double));

      Atmp = (double **)malloc(B*sizeof(double *));
      for (ii=0; ii<B; ii++)
         Atmp[ii] = (double *)malloc(N*sizeof(double));
   
      t_start=usecs();
      //genmat();
      genrandval();
      genmat_geostat_stuff();
      t_end=usecs();
      time = ((double)(t_end-t_start))/1000000;
      printf("process %d/%d, time to randval and gen_cova_mat = %f (N=%d,B=%d)\n", whoAmI,howmany, time,N,B);


      MPI_Barrier(MPI_COMM_WORLD);


      t_start=usecs();

      for (kk=0; kk<NB; kk++) {
         if(kk%howmany == whoAmI){
            lu0localChol(kk);

#pragma omp parallel
{
#pragma omp for schedule(runtime) private(jj) nowait
            for (jj=kk+1; jj<NB; jj++)
               fwdlocalChol(kk, jj);
}
  
            //sending the blocks to the other threads
            for(p=0;p<howmany;p++){
               if(p!=whoAmI){
                  sendBlocks(kk,p,req[p]);
               }
            }
            for(p=0;p<howmany;p++){
               if(p!=whoAmI){
                  for(i=0;i<B;i++)
                     MPI_Wait(&req[p][i],&st);
               }
            }
            //end sending

         }
         else{
            recvBlocks(kk,kk%howmany,req2[whoAmI]);
         }

         for (ii=kk+1; ii<NB; ii++) {
               //bdiv (ii, kk);
            if(ii%howmany == whoAmI){
               bdivlocalChol(ii, kk);
               for (jj=kk+1; jj<NB; jj++) 
                  bmodlocalChol(ii, jj, kk);
            }

         }
      }
      t_end=usecs();

      printres(0);

      time = ((double)(t_end-t_start))/1000000;
      //if(whoAmI==0){
         //printf("processes %d, time to compute = %f (N=%d,B=%d)\n", howmany, time,N,B);
         printf("process %d/%d, time to compute = %f (N=%d,B=%d)\n", whoAmI,howmany, time,N,B);
      //}
   



   MPI_Barrier(MPI_COMM_WORLD);


   //for (ii=0; ii<SIMS; ii++)
   for (ii=0; ii<nrows*B; ii++)
      free(simval[ii]);
   free(simval);

   free(randval);

   for (ii=0; ii<nrows*B; ii++)
      free(A[ii]);
   free(A);

   for (ii=0; ii<B; ii++)
      free(Atmp[ii]);
   free(Atmp);   

   delete_geostat_stuff();

   MPI_Finalize();
}


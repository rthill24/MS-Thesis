//////////////////////////////////////////////////////////////////////
//Rapid implementation of the Smith Collapse method for Python
//Based on original Python and C implementation by M. Lankowski
//Extended by Matt Collette
// (c) 2012 Regents of the University of Michigan
//This code is distributed under an open-source licence please see the 
//include file License for detail
///////////////////////////////////////////////////////////////////////
#include "C:/Program Files (x86)/GnuWin32/include/gsl/gsl_spline.h"
#include "C:/Program Files (x86)/GnuWin32/include/gsl/gsl_errno.h"
#include "C:/Program Files (x86)/GnuWin32/include/gsl/gsl_math.h"
#include "C:/Program Files (x86)/GnuWin32/include/gsl/gsl_min.h"

typedef struct compT
{
	//Represents a single element in the Smith collapse methodology
	int id;			  //Integer ID for the element
	double x_global;  //Global x-location
	double y_global;  //Global y-location
	double y_eff;     //effective 
	int epp;		  //Defines if elastic-perfectly plastic is assumed
	                  //On tension side (1) or if spline data should 
	                  //used for all points (0)
	double E;		  //Elastic tension modulus
	double Ys;        //Yield stress used if epp=1
	double Area;	  //Cross sectional area of element
	int numPts;		  //Number of points in stress-strain curve
	double maxStrain; //Maximum (tensile, positive) strain
	double minStrain; //Minimum (compressive, negative) strain
	//Following are GSL elements for spline
	gsl_interp_accel *acc_stress;
	gsl_spline *spline_stress;
	
	//Pointer for next item in cross section
	struct compT* next_comp;
} comp;


typedef struct 
{
	//Simple structure for creating cross-section
    int numPanels;
    comp* firstComp;
    comp* lastComp;
} crossSection;



//Simple adaptor struct to find the maximum point the the load-
//shortening curve
typedef struct 
{
	crossSection *section;		//Section to apply curvature to 
    double NApos;				//NA position
    double fTol;			//Tolerance for force inbalance
    double *elementStrains;	 //Array of element strains
    double actualMoment;		//Actual moment with curvature
                                //returned moment is sign-adapted
    
} momentAdaptor;



crossSection* CreateEmptyCross(void);

comp* CreateAndAddComponent( int id,
					double x_global, 
                    double y_global, 
                    int epp, 
                    double E,
                    double Ys,
                    double Area,
                    int numPts,
                    double stress[], 
                    double strain[],
                    crossSection *section);
                   
void testPanels(crossSection *section);

void cleanUp(crossSection* section);

double getCompForce(comp* element, double strain);

void getForMom(crossSection *section, double radcurvature,
                  double NApos, double* netForce, double* netMom,
                  double *NAShift, double elementStrains[]);


double momentPureBending(crossSection *section, 
                         double radcurvature,
                         double *NApos, const double fTol,
                         double elementStrains[]);


double momentGSL(double radcurvature, void* params);

void setupVerticalBending(crossSection *section);

void doZeroCurvature(crossSection *section, 
                         double *NApos,
                         double elementStrains[]);

int findMaxMom(crossSection *section, double minCurv, double maxCurv,
				double forceTol, double curvTol, 
                double* NA, double *moment, double* actualCurv,
                double elementStrains[]);

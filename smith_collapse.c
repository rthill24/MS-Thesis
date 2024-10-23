//////////////////////////////////////////////////////////////////////
//Rapid implementation of the Smith Collapse method for Python
//Based on original Python and C implementation by M. Lankowski
//Extended by Matt Collette
// (c) 2012 Regents of the University of Michigan
//This code is distributed under an open-source licence please see the 
//include file License for detail
///////////////////////////////////////////////////////////////////////



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "smith_collapse.h"

/*!
 * This collection of functions implements a simple Smith progessive
 * collapse of ship hull girder or other box-beam like structure
 * 
 * Sign convention is - distances positive up from keel (assumed origin)
 * tensile strains and stresses positive
 * hogging bending moments positive
 * Postive radius of curvature gives sagging (negative) response
 * 
 * Should work for any consistent stress, distance, and component cross
 * sectional area units.  Returned moments have units of distance*area*
 * stress 
 */

//! Creates and empty cross section in memory
crossSection* CreateEmptyCross(void)
{
	crossSection *newCross;
	newCross = (crossSection *)malloc(sizeof(crossSection));

	//Check for failed allocation
	if (newCross == 0)
	{
		printf("ERROR: Out of memory to allocate new cross section\n");
		printf("Returning void pointer and hoping for the best...\n");
		return 0;
	}
	newCross->firstComp = 0;
	newCross->lastComp = 0;
	return newCross;
}
//! Creates a new component and appends it to the provided cross section
/*!
\param id - integer ID for the element will be used in printouts
\param x_global Global x location of centriod of the component 
\param y_global Global y location of centriod of the component 
\param epp False(0) or True(1) for elastic-perfectly plastic behavior
* if true, provided stress-strain curve is only used for compression and
* elastic--perfectly plastic behavior is assumed for tensile behavior
\param E elastic modulus of base material.  Used for neutral axis calcs
\param Ys Yield stress of material, used with epp=1
\param Area cross-sectional area of component
\param numPts - number of points in stress-strain curve
\param strain - array of stress strain curve strain values  - must be
* strictly increasing from most negative (compressive) to most positive
* (tensile).  This is NOT CHECKED explicitly.
\param stress - array of stress-strain curve stress values as above
\param section - cross section to append the component to - optional
\return The component object.  Will be null pointer if error or no mem
*/
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
                    crossSection *section)
{
	comp *newComp;
	newComp = (comp *)malloc(sizeof(comp));
	if (newComp == 0)
	{
		printf("ERROR: Out of memory to allocate new comp\n");
		printf("Returning void pointer and hoping for the best...\n");
		return 0;
	}


	//Copy over standard definitions
	newComp->id = id;
	newComp->x_global = x_global;
	newComp->y_global = y_global;
	newComp->epp = epp;
	newComp->E = E;
	newComp->Ys = Ys;
	newComp->Area = Area;
	newComp->next_comp = 0;
	newComp->numPts = numPts;
	newComp->maxStrain = strain[numPts-1];
	newComp->minStrain = strain[0];
	
	//Try to set up the interpretive splines
	//accelerator object first
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    if (!acc)
 	{
		printf("ERROR: Out of memory to allocate new comp\n");
		printf("failed during spline accelerator\n");
		printf("Returning void pointer and hoping for the best...\n");
		return 0;
	}   
	else
	{
		newComp->acc_stress = acc;
	}
	gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, numPts);
    if (!spline)
 	{
		printf("ERROR: Out of memory to allocate new comp\n");
		printf("failed during spline creation\n");
		printf("Returning void pointer and hoping for the best...\n");
		return 0;
	}   
	else
	{
		int retCode = gsl_spline_init(spline, strain, stress, numPts);
		newComp->spline_stress = spline;
	}     
         
    //Insert into linked list for section if passed
    if (section)
    {
		if (section->firstComp == 0) //No comps yet registered
		{
			section->firstComp = newComp;
			section->lastComp = newComp;
		}
		else
		{
			section->lastComp->next_comp  = newComp;
			section->lastComp = newComp;
		}
	}
	return newComp;
		
}
//! Really simple testing function 
void testPanels(crossSection *section)
{
	comp *currComp;
	currComp = section->firstComp;
	while (currComp)
	{
		printf("Comp x y location is (%f %f)", currComp->x_global, 
		       currComp->y_global);
		currComp = currComp->next_comp;
	}
	return;
	
}


//! Frees all memory associated with a component
void cleanUpComp(comp* toRemove)
{
	//Remove the spline if present
	if (toRemove->acc_stress)
	{
		gsl_interp_accel_free (toRemove->acc_stress);
	}
	if (toRemove->spline_stress)
	{
		gsl_spline_free (toRemove->spline_stress);
	}
	free(toRemove);
	return;
}


//! Frees all memory associated with section AND all of its components
/*!
\param section - cross section pointer, will be invalid when complete
*/
void cleanUp(crossSection* section)
{
	comp *currComp; 
	comp *temp;
	currComp = section->firstComp;
	while (currComp)
	{
		temp = currComp;
		currComp = currComp->next_comp;
		cleanUpComp(temp);
	}
	free(section);
	return;	
}

//! Gets the force associated with the strain on a single component
/*!
\param element - component of cross-section
\param strain - local strain value
*/
double getCompForce(comp* element, double strain)
{
	double result;
	#ifdef DEBUG
		printf("In getforce with element ID %i\n", element->id);
		printf("Area is %f\n", element->Area);
	#endif
	if ((element->epp == 1) && (strain >= 0)) //Use EPP formulation
	{
		#ifdef DEBUG
			printf("Using EPP");
		#endif
		if (strain*element->E > element->Ys) //Above Ys, yield*area
		{
			result = element->Ys*element->Area;
		}
		else  //Below YS, use hooke's law to get net force
		{
			result = strain*element->E*element->Area;
		}
	}
	else  //We have to use GSL spline
	{
		#ifdef DEBUG
			printf("using GPL spline\n");
		#endif
		if (strain > element->maxStrain) //Out of bounds
		{
			printf(
			      "Warning strain out of range of element s/s curve\n");
			printf("For element id %d \n", element->id);
			printf("Provide strain was  %f\n", strain);
			printf("Maximum strain was %f\n", element->maxStrain);
			printf("Assuming stress same as last valid point....\n");
			result = gsl_spline_eval(element->spline_stress, 
			                       element->maxStrain,
			                       element->acc_stress)*element->Area;
			                      
		}
		else if (strain < element->minStrain) //Out of bounds
		{
			/*printf(
			      "Warning strain out of range of element s/s curve\n");
			printf("For element id %d \n", element->id);
			printf("Provide strain was  %f\n", strain);
			printf("Minimum strain was %f\n", element->minStrain);
			printf("Assuming stress same as last valid point....\n");*/
			result = gsl_spline_eval(element->spline_stress, 
			                       element->minStrain,
			                       element->acc_stress)*element->Area;
			                       			
		}
		else
		{	
			result = gsl_spline_eval(element->spline_stress, 
			                       strain,
			                       element->acc_stress)*element->Area;


			                       			
		}  
		
	}
	return result;
}

//! Calculates the net axial force and moment from a given neutral axis 
//! position
//! See note about sign convention
//! Also calculated the suggested NA shift by the Rutherford method
//! extending to use the secant modulus of the materials
/*! 
 * 
\param section - cross-section to compute 
\param radcurvature - radius of curvature to be applied to the section
\param NApos - NA position in y_eff coordinates 
\param netForce - pointer, will be updated with net force in section
\param netMom - pointer, will be updated with net moment in section
\param NAShift - recommend shift in the neutral axis by Rutherford
\param elementStrains - array of all the element strains at this radc
*/ 
void getForMom(crossSection *section, double radcurvature,
                   double NApos, double* netForce, double* netMom,
                  double *NAShift, double elementStrains[])
{
	//Reset variables
	double force = 0;
	double mom = 0;
	
	//Go through compoents on section and calculate strain and total
	//up force
	comp *currComp;
	double currStrain, currForce;
	//We will build up the top and bottom term's of Rutherford's shift
	//formula as we go - top term is just the force!
	double NAShiftBottomTerm = 0;
	currComp = section->firstComp;
	int kount = 0;
	while (currComp)
	{
		//Strain - positive curvatures leads to sagging 
		currStrain = -(currComp->y_eff - NApos)/radcurvature;
		elementStrains[kount] = currStrain;
		//Force
		currForce = getCompForce(currComp, currStrain);
		force = force + currForce;
		//Moment
		mom = mom + currForce*(currComp->y_eff - NApos);
		//Bottom term - secant modulus (s/str) * A  - just force/str
		//but handle zero-strain correctly 
		if (fabs(currStrain) > pow(10.,-9))
		{
			NAShiftBottomTerm += fabs(currForce/currStrain);
		}
		else
		{
			NAShiftBottomTerm += (currComp->E)*(currComp->Area);
		}
		currComp = currComp->next_comp;
		kount += 1;
	}	
	*netForce = force;
	*netMom = mom;
	*NAShift = -radcurvature*force/NAShiftBottomTerm;
	return;
}


//! Special routing for calculating the NA position and element
//! strains for zero curvature that could cause all sorts of converg.
//! issues for normal mechanism
/*! 
 * 
\param section - cross-section to compute 
\param NApos - NA position in y_eff coordinates for zero net force 
no guess required, will be set to elastic NA location
\param elementStrains - array of final element strains will be filled
* with zeros
*/
void doZeroCurvature(crossSection *section, 
                         double *NApos,
                         double elementStrains[])
{
	double areaTotal = 0;
	double firstMom = 0;
	double Eref;
	comp *currComp;
	int kount = 0;
	currComp = section->firstComp;
	if (currComp)
	{
		Eref = currComp->E;
		while (currComp)
		{
			areaTotal += currComp->Area*currComp->E/Eref;
			firstMom += currComp->Area*currComp->E/Eref*currComp->y_eff;
			elementStrains[kount] = 0;
			kount += 1;
 			currComp = currComp->next_comp;
		}
		*NApos = firstMom/areaTotal;
	}
	else
	{
		printf(
		"Warning - zero curvature properties passed empty section");
	}
	return;
}

//! Calculates moment for an applied curvature.  Also returns the 
//! Neutral axis location.  
/*! 
 * 
\param section - cross-section to compute 
\param curvature - curvature to be applied to the section - see signs
\param NApos - NA position in y_eff coordinates for zero net force 
should be filled with a guess value for the first iteration
\param fTol - force tolerance for convergence of NA position
\param elementStrains - array of final element strains

*/ 
double momentPureBending(crossSection *section, 
                         double radcurvature,
                         double *NApos, double fTol,
                         double elementStrains[])
{
	//Handle the zero curvature case as a special case
	if (radcurvature == 0.)
	{
		doZeroCurvature(section, NApos, elementStrains);
		return 0;
	}
	double unbalancedForce = 0; //remaining unbalanced force
	double moment = 0;          //Resisting moment
	double NAShift = 0;         //Neutral axis shift recommended
	double NAorig = *NApos;     //Hold for a restart
	
	//Some tracking variables
	double alpha = 1.0;   //step size 1.0=perfect rutherford correction
	int kount = 0;        //How many iterations we have done                       
	const int ITERLIMIT = 500;  //Limit on how long we will iterate
	
	do
	{
		getForMom(section,radcurvature, *NApos, &unbalancedForce,
		          &moment, &NAShift, elementStrains);
		*NApos += NAShift*alpha;
	//	printf("NA_Pos = %f\n", *NApos);
		kount += 1;
		#ifdef DEBUG
			printf("Moment and kount %f %i \n", moment, kount);
		#endif
	}
	while (
	    (fabs(unbalancedForce) > fabs(fTol)) && (kount <= ITERLIMIT));
	
	//If this does not work try one more time with a smaller step size
	if (kount > ITERLIMIT) 
	{
		alpha = 0.5;
		kount = 0;
		*NApos = NAorig;
		do
		{
			getForMom(section,radcurvature, *NApos, &unbalancedForce,
		          &moment, &NAShift, elementStrains);
			*NApos += NAShift*alpha;
			kount += 1;
		}
		while ((fabs(unbalancedForce) > fabs(fTol)) && 
		      (kount <= ITERLIMIT));
		if (kount > ITERLIMIT)
		{
			//Ideally we woudl add a different fall back here
			//Something like a Brent root finder
			printf("Smith method unable to find neutral axis location");
			printf("Supplied radius of curvature was %f", radcurvature);
			printf("Error currently unhandled, program stopping...");
			exit(1);
		}		
		
	}
	#ifdef DEBUG
		printf("Iterations taken within zero moment %i\n", kount);
	#endif 
	
	return moment;
}


//internal adaptor function to allow momentPureBending to be called
//by GSL minimization library
double momentGSL(double curvature, void* params)
{
	momentAdaptor *ad;
	
	double radcurvature = 1./curvature;
	
	ad = (momentAdaptor*)params;
	
	//Call the function
	double temp = momentPureBending(ad->section,radcurvature,
                         &(ad->NApos), ad->fTol, ad->elementStrains);
    
    //Save the correctly signed value for later use
    ad->actualMoment = temp;
    
    //Return in the form for minimization
    return -fabs(temp);
   
}

//! Prepares the crossSection for vertical bending ony by setting 
//! y_eff to y_global
//! @TODO - update for bending about any axis by taking Cx, Cy
/*! 
 *
\param section - cross-section to update
*/

void setupVerticalBending(crossSection *section)
{
	comp *currComp;
	currComp = section->firstComp;
	while (currComp)
	{
		currComp->y_eff = currComp->y_global;
		currComp = currComp->next_comp;
	}
	return;
}


//! Calculates the maximum moment between twocurvature values
//! if there are multiple peaks, which peak is found is not gaurunteed 
/*! 
 * 
\param section - cross-section to compute 
\param minCurv - minumum curvature (NOT radius) at this point
\param maxCurv - maximum curvature (NOT radiues) at this point
\param forceTol - equilibrium force tolerance in floating point for each
* curvature calculation
\param curveTol - overall tolerance on the curvature for convergence to
* the maximum bending moment
\param NA - starting guess for Neutral axis position, 
* will then have final NA position noted
\param moment - resulting maximum resisting moment
\param actualCurv - curvature (NOT radius) at which this occurs
\param elementStrains - array of final element strains
\return Error code - 0 equals success, higher indicates a problem
*/ 
int findMaxMom(crossSection *section, double minCurv, double maxCurv,
				double forceTol, double curvTol, 
                double* NA, double *moment, double* actualCurv,
                double elementStrains[])
{
	
	//Return code
	int retCode = 0;
	//Construct the parameter object for the minimizer
	momentAdaptor Adaptor;
	Adaptor.section = section;
	Adaptor.NApos = *NA;
	Adaptor.fTol = forceTol;
	Adaptor.elementStrains = &elementStrains[0];
	
	//Set up for the optimizer
	int status;    			//GSL return status code
      int check;                  //Error Code for gsl_fminimizer_set
	int iter = 0;			//Current iteration
	const int MAXIT = 500;	//Maximum times to look for min
	
	const gsl_min_fminimizer_type *T;
     gsl_set_error_handler_off ();
	gsl_min_fminimizer *s;
	gsl_function F;
	F.function = &momentGSL;
	F.params = &Adaptor;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	check = gsl_min_fminimizer_set (s, &F, (minCurv+maxCurv)/2., 
	                        minCurv, maxCurv);
	if (check ==  GSL_EINVAL)
	{
            return 4;			
	}    
     
	do
	{
		//Do the iteration and update the values
		iter++;
		status = gsl_min_fminimizer_iterate (s);
		*actualCurv = gsl_min_fminimizer_x_minimum (s);
		minCurv = gsl_min_fminimizer_x_lower (s);
		maxCurv = gsl_min_fminimizer_x_upper (s);
		
          #ifdef DEBUG
			printf("Max locator iteration %i\n", iter);
			printf("Lower bound is %14.6g with value %14.6e\n", minCurv, 
				      gsl_min_fminimizer_f_lower(s));
			printf("upper bound is %14.6g with value %14.6e\n", minCurv, 
				      gsl_min_fminimizer_f_upper(s));
			printf("minimum is %14.6g with value %14.6e\n", *actualCurv, 
				      gsl_min_fminimizer_f_minimum(s));			
          #endif 
		//Do some error checking
        	
		if (status == GSL_EBADFUNC)
		{
			printf("Was unable to calcaluate required moment for\n");
			printf("curvature min/max was %f %f\n", minCurv, maxCurv);
			printf("returning error code");
			break;
		}
		else if (status ==  GSL_FAILURE)
		{
			printf("Smith method maximizer failed on iterval\n");
			printf("curvature min/max was %f %f\n", minCurv, maxCurv);
			printf("returning error code");
			break;			
		}
		else if (status ==  GSL_EINVAL)
		{
			printf("BEAKING\n");
			break;				
		}                        
		//See if we are done
		else if (fabs(maxCurv - minCurv) < curvTol)
		{
			status = GSL_SUCCESS;
		}
		else
		{
			status = GSL_CONTINUE;
		}

	} while (status == GSL_CONTINUE && iter < MAXIT);
	
	//If all worked, just call function one final time to be sure
	//last call was at the minimum so NA, elementStrains all updated
	if (status == GSL_SUCCESS)
	{
		momentGSL(*actualCurv, &Adaptor);
		*NA = Adaptor.NApos;
		*moment= Adaptor.actualMoment;
		//*actualCurv is up-to-date from above
		
	}
	else
	{
		if (iter >= MAXIT)
		{
			printf("Smith method maximizer failed to converge\n");
			printf("over iteration limit");		
			retCode = 1;	
		}
		else
		{
			printf("Error in Smith method maximizer");
			retCode = 1;
		}
	}
	gsl_min_fminimizer_free (s);
	return retCode;
}

                         




		

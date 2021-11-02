#include <iostream>
#include <cmath>
#include <string.h>
#include<float.h>
#include<complex.h>
#include <vector>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include <globes/globes.h>
#include<fstream>

extern "C"
{
	#include "snu.h"
}

using namespace std;

char AEDLFILE[] = "lbne_1300.glb";

int main(int argc, char * argv[])
{

	glbInit(argv[0]);
	glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);

	ofstream out;
	out.open("probability_epsmue.dat");

	double theta12 = asin(sqrt(0.307));
	double theta13 = 0.5*asin(sqrt(0.09));
	double theta23 = 0.5*asin(sqrt(1.0));
	double deltacp = 0.0 * M_PI/180.0;
	double sdm = 7.54e-5;
	double ldm = 2.43e-3;

	snu_init_probability_engine_3();

	glbRegisterProbabilityEngine(6 * 9 - 3,
                               &snu_probability_matrix,
							   &snu_set_oscillation_parameters,
  							   &snu_get_oscillation_parameters,
  							   NULL);
	/* Define "true" oscillation parameter vector */
	glb_params true_values = glbAllocParams();
  
	for(unsigned int i=0; i < 51; i++)
	{
	glbSetOscParams(true_values, 0.0, i);
	}

	glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);

	//############ NSI Parameter #################################//
	glbSetOscParams(true_values, 0.0, 24);   // eps_ee 
	glbSetOscParams(true_values, 0.1, 25);  // eps_mue magnitude
        glbSetOscParams(true_values, 0.0, 26);  // eps_mue phase
        glbSetOscParams(true_values, 0.0, 27);  // eps_etau 
        glbSetOscParams(true_values, 0.0, 28);  // eps_etau phase
        glbSetOscParams(true_values, 0.0, 29);  // eps_mumu
        glbSetOscParams(true_values, 0.0, 30);  // eps_mutau
        glbSetOscParams(true_values, 0.0, 31);  // eps_mutau phase
        glbSetOscParams(true_values, 0.0, 32);  // eps_tautau
	glbSetDensityParams(true_values,1.0,GLB_ALL);

	glb_params input_errors = glbAllocParams();
  	
	glbSetDensityParams(input_errors, 0.05, GLB_ALL);
	glbSetInputErrors(input_errors);

	glbSetOscillationParameters(true_values);
	glbSetRates();

	double energy, prob_NH;
	double emin= 0.5 ; //GeV
	double emax=20 ; //GeV
	double step= 1000;
  
	for (energy=emin;energy<=emax;energy+=(emax-emin)/step)
	{
   
	glbSetOscillationParameters(true_values);
	prob_NH=glbProfileProbability(0,1,2,+1,energy);
 

	out<<energy<<"  "<<prob_NH<<endl;
	}


	out.close();
	glbFreeParams(true_values);
 	return 0;

}

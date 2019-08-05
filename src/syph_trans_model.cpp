#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


/// function for matrix multiplication (%*% in R)

// [[Rcpp::export]]
  arma::vec mv_mult(
      arma::mat& lhs, 
      arma::vec& rhs) {
          return lhs * rhs;
  }

  // [[Rcpp::export]]
  arma::mat mat_mult(
      arma::mat& a, 
      arma::mat& b) { 
        return(a % b); 
    }

/// main transmission model code 

//[[Rcpp::export]]  
List syphSim(
    List  x, //params
    double dt, // timestep
    List cm, //mixing 
    NumericVector pabx, //p.abx
    NumericVector rep_count, // when to turn counting of prior infections on 
    int nYrs, //model run time
    NumericVector initPop, //initial population
    NumericVector n_sa, //size of sexually active population
    arma::mat births, //birth rate
    arma::mat births_sa, //births p_sexually active
    arma::mat births_nsa, //births p_not sexually active
    arma::mat aging, //aging rate
    arma::mat aging_nsa //aging p.not sexually active
    
) {
  int i = 5; //number of subpopulations
  int j = 2; //number of sexes
  int k = 2; //number of sexual activity groups
  int l = 2; //number of age groups
  int dim = i*j*k*l; //total number of population subgroups in output matrix
  NumericVector b = as<NumericVector>(x["b"]); //transmission rate
  double delta = as<double>(x["delta"]); //1/dur incubation 
  double return_rate_for_resusceptible = as<double>(x["return_rate_for_resusceptible"]); // rate of return to first-time-susceptible from resusceptible
  NumericVector gamma = as<NumericVector>(x["gamma"]); //1/dur infectious
  NumericVector trt1 = as<NumericVector>(x["p.trt.1"]); //treatment rate, primary syph
  NumericVector trt2 = as<NumericVector>(x["p.trt.2"]); //treatment rate, secondary syph
  NumericVector trt3 = as<NumericVector>(x["p.trt.3"]); //treatment rate, early latent syph
  NumericVector trt4 = as<NumericVector>(x["p.trt.4"]); //treatment rate, late latent syph
  NumericVector report = as<NumericVector>(x["rep"]); //reporting probability, screened case
  NumericMatrix report_symp = as<NumericMatrix>(x["rep.symp"]); //reporting probability matrix, case seeking care
  NumericMatrix screen = as<NumericMatrix>(x["screen"]); //screening rate matrix
  NumericVector alpha(dim); //screening rate (updated annually)
	NumericMatrix screen_repeat = as<NumericMatrix>(x["screen_repeat"]); // screening rates for individuals with diagnosis history
	NumericVector alpha_repeat(dim); // screening rate (updated annually) for individuals with diagnosis history 
  NumericVector rep_s(dim); //reporting probabilities (updataed annually)
  NumericVector c_msm = as<NumericVector>(x["behav"]); //transmission rr in MSM 
  NumericVector dur_imm = as<NumericVector>(x["dur.imm"]); //duration immunity after treatment
  NumericVector cmlow1 = as<NumericVector>(cm["cm.low1"]); //contact matrix low activity, black 
  NumericVector cmlow2 = as<NumericVector>(cm["cm.low2"]); //contact matrix low activity, other 
  NumericVector cmlow3 = as<NumericVector>(cm["cm.low3"]); //contact matrix low activity, Hispanic 
  NumericVector cmlow4 = as<NumericVector>(cm["cm.low4"]); //contact matrix low activity, HIV- MSM
  NumericVector cmlow5 = as<NumericVector>(cm["cm.low5"]); //contact matrix low activity, HIV+ MSM
  NumericVector cmhigh1 = as<NumericVector>(cm["cm.high1"]); //contact matrix high activity, black 
  NumericVector cmhigh2 = as<NumericVector>(cm["cm.high2"]); //contact matrix high activity, other 
  NumericVector cmhigh3 = as<NumericVector>(cm["cm.high3"]); //contact matrix high activity, Hispanic 
  NumericVector cmhigh4 = as<NumericVector>(cm["cm.high4"]); //contact matrix high activity, HIV- MSM
  NumericVector cmhigh5 = as<NumericVector>(cm["cm.high5"]); //contact matrix high activity, HIV+ MSM
	bool output_weekly = as<LogicalVector>(x["output_weekly"])[0] == TRUE; // output every timestep or only yearly (default)
  arma::mat birth_rate = births; //birth rate
  arma::mat aging_rate = aging; //aging rate
  arma::vec S(dim);  //susceptible
  arma::vec E(dim);  //incubating
  arma::vec I1(dim); //primary syphilis
  arma::vec I2(dim); //secondary syphilis
  arma::vec L1(dim); //early latent syphilis
  arma::vec L2(dim); //late latent syphilis
  arma::vec T1(dim); //treated, primary and secondary
  arma::vec T2(dim); //treated, early latent
  arma::vec T3(dim); //treated, late latent
  arma::vec SR(dim); //susceptible, prior treated infection
  arma::vec ER(dim); //incubating, prior treated infection
  arma::vec IR1(dim); //primary, prior treated infection
  arma::vec IR2(dim); //secondary, prior treated infection
  arma::vec LR1(dim); //early latent, prior treated infection
  arma::vec LR2(dim); //late latent, prior treated infection
  arma::vec Inc(dim); //cumulative incidence
  arma::vec Incr(dim); //cumulative incidence, prior infection
  arma::vec D1(dim);  //reported cases, primary 
  arma::vec D2(dim);  //reported cases, secondary
  arma::vec D3(dim);  //reported cases, early latent
  arma::vec NSA(dim); //not sexually active 
  arma::vec D4(dim);  //reported cases, late latent
  arma::vec DR(dim);  //reported cases, primary, secondary & early latent, prior infection
  arma::vec p_trt1(dim); //trt rate, primary
  arma::vec p_trt2(dim); //trt rate, secondary
  arma::vec p_trt3(dim); //trt rate, early latent
  arma::vec p_trt4(dim); //trt rate, late latent
  arma::vec lambda(dim); //force of infection
  arma::vec Pop(dim*24); //population matrix
	
	int output_length; 
	if (output_weekly) {
		output_length = nYrs * 52;
	} else {
		output_length = nYrs;
	}
  NumericMatrix outputs(output_length, dim*24+1);
  NumericMatrix out(output_length, dim*24+1);
  
  double gamma1 = gamma[0]; //1/dur primary
  double gamma2 = gamma[1]; //1/dur secondary
  double gamma3 = gamma[2]; //1/dur latent
  double lambda1 = 365/dur_imm[0]; //1/dur imm primary and secondary
  double lambda2 = 365/dur_imm[1]; //1/dur imm early latent
  double lambda3 = 365/dur_imm[2]; //1/dur imm late latent
  
  //update treatment parameters for males
  for(int n=0; n<((i-2)*k*l); n++) {
    p_trt1[n] = trt1[0];
    p_trt2[n] = trt2[0];
    p_trt3[n] = trt3[0];
    p_trt4[n] = trt4[0];
  }
  //update treatment parameters for HIV- MSM
  for(int n=((i-2)*k*l); n<((i-1)*k*l); n++) {
    p_trt1[n] = trt1[1];
    p_trt2[n] = trt2[1];
    p_trt3[n] = trt3[1];
    p_trt4[n] = trt4[1];
  }
  //update treatment parameters for HIV+ MSM
  for(int n=((i-1)*k*l); n<i*k*l; n++) {
    p_trt1[n] = trt1[2];
    p_trt2[n] = trt2[2];
    p_trt3[n] = trt3[2];
    p_trt4[n] = trt4[2];
  }
  //update treatment parameters for females
  for(int n=((i*k*l)); n<(2*i*k*l); n++) {
    p_trt1[n] = trt1[3];
    p_trt2[n] = trt2[3];
    p_trt3[n] = trt3[3];
    p_trt4[n] = trt4[3];
  }
  
  /// year loop
  for(int y=0; y<nYrs; y++){
    for(int w=0; w<52; w++){

			/* for (int d=0; d<7; d++){ */
      
				//updata reporting (in screened), transmission rr, and background abx use
				double c_rr = c_msm[y];
				double rep = report[y];
				double abx = pabx[y];
				double rep_on = rep_count[y];
				
				if(y==0 && w==0){   ///populate the model at t=0
					for(int n=0; n<dim*24; n++){
						Pop[n] = initPop[n];
					}
				}
				//int s = y*52+w; //number of weeks since start
				else{
				for(int n=0; n<dim; n++) { //calculate population size at start of each week
					S[n] =   Pop[n];
					E[n] =   Pop[n+dim];
					I1[n] =  Pop[n+2*dim];
					I2[n] =  Pop[n+3*dim];
					L1[n] =  Pop[n+4*dim];
					L2[n] =  Pop[n+5*dim];
					T1[n] =  Pop[n+6*dim];
					T2[n] =  Pop[n+7*dim];
					T3[n] =  Pop[n+8*dim];
					SR[n] =  Pop[n+9*dim];
					ER[n] =  Pop[n+10*dim];
					IR1[n] = Pop[n+11*dim];
					IR2[n] = Pop[n+12*dim];
					LR1[n] = Pop[n+13*dim];
					LR2[n] = Pop[n+14*dim];
					Inc[n] = Pop[n+15*dim];
					Incr[n] = Pop[n+16*dim];
					D1[n] =  Pop[n+17*dim];
					D2[n] =  Pop[n+18*dim];
					D3[n] =  Pop[n+19*dim];
					NSA[n] = Pop[n+20*dim];
					D4[n] =  Pop[n+21*dim];
					DR[n] =  Pop[n+22*dim];
					alpha[n] = screen(y,n);  //update screening and reporting
					alpha_repeat[n] = screen_repeat(y,n);
					rep_s[n] = report_symp(y,n);
				}
				
				//calculate force of infection for each sex and subpopulation
				for(int n=0; n<4; n++) { //males subpop 1
					int m = n;
					lambda[n] =  (b[0] * (cmlow1[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow1[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow1[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow1[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow1[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow1[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh1[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) +  //F, high AC, age=1, subpop1
																cmhigh1[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) +  //F, high AC, age=2, subpop1
																cmhigh1[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh1[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh1[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh1[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=4; n<8; n++) {   //males subpop 2
					int m = n-4;
					lambda[n] =  (b[0] * (cmlow2[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow2[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow2[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow2[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow2[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow2[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh2[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh2[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh2[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh2[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh2[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh2[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=8; n<12; n++) {   //males subpop 3
					int m = n-8;
					lambda[n] =  (b[0] * (cmlow3[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow3[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow3[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow3[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow3[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow3[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh3[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh3[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh3[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh3[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh3[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh3[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=12; n<16; n++) {   //males subpop 4
					int m = n-12;
					lambda[n] =  ((b[0] *  (cmlow4[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																	cmlow4[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																	cmlow4[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																	cmlow4[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																	cmlow4[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																	cmlow4[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																	cmhigh4[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																	cmhigh4[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																	cmhigh4[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																	cmhigh4[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																	cmhigh4[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																	cmhigh4[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))) +//F, high AC, age=2, subpop3
								
											((b[2]*c_rr)*
																 (cmlow4[m+24]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																	cmlow4[m+28]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																	cmlow4[m+32]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																	cmlow4[m+36]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																	cmhigh4[m+24] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																	cmhigh4[m+28] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																	cmhigh4[m+32] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																	cmhigh4[m+36] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19])))); //M, high AC, age=2, subpop5
										
						}  
				
				for(int n=16; n<20; n++) {   //males subpop 5
					int m = n-16;
					lambda[n] =  ((b[0] *(cmlow5[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow5[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow5[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow5[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow5[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow5[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh5[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh5[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh5[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh5[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh5[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh5[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))) + //F, high AC, age=2, subpop3
																
													((b[2]*c_rr)*
															 (cmlow5[m+24]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow5[m+28]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow5[m+32]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow5[m+36]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																cmhigh5[m+24] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh5[m+28] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh5[m+32] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh5[m+36] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19])))); //M, high AC, age=2, subpop5
															
				}
				
				for(int n=20; n<24; n++) {   //females subpop 1
					int m = n;
					lambda[n] =  (b[1] * (cmlow1[m+20]     * (I1[0]  + I2[0]  + IR1[0]  + IR2[0])  /(n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow1[m+4+20]   * (I1[2]  + I2[2]  + IR1[2]  + IR2[2])  /(n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow1[m+8+20]   * (I1[4]  + I2[4]  + IR1[4]  + IR2[4])  /(n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow1[m+12+20]  * (I1[6]  + I2[6]  + IR1[6]  + IR2[6])  /(n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow1[m+16+20]  * (I1[8]  + I2[8]  + IR1[8]  + IR2[8])  /(n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow1[m+20+20]  * (I1[10] + I2[10] + IR1[10] + IR2[10]) /(n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow1[m+24+20]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) /(n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow1[m+28+20]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) /(n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow1[m+32+20]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) /(n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow1[m+36+20]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) /(n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh1[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh1[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh1[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh1[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh1[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh1[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh1[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh1[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh1[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh1[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				for(int n=24; n<28; n++) {   //females subpop 2
					int m = n-4;
					lambda[n] =  (b[1] * (cmlow2[m+20]     * (I1[0]  + I2[0]  + IR1[0]  + IR2[0])  / (n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow2[m+4+20]   * (I1[2]  + I2[2]  + IR1[2]  + IR2[2])  / (n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow2[m+8+20]   * (I1[4]  + I2[4]  + IR1[4]  + IR2[4])  / (n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow2[m+12+20]  * (I1[6]  + I2[6]  + IR1[6]  + IR2[6])  / (n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow2[m+16+20]  * (I1[8]  + I2[8]  + IR1[8]  + IR2[8])  / (n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow2[m+20+20]  * (I1[10] + I2[10] + IR1[10] + IR2[10]) / (n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow2[m+24+20]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow2[m+28+20]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow2[m+32+20]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow2[m+36+20]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh2[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh2[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh2[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh2[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh2[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh2[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh2[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh2[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh2[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh2[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				for(int n=28; n<32; n++) {   //females subpop 3
					int m = n-8;
					lambda[n] =  (b[1] * (cmlow3[m+20]    * (I1[0]   + I2[0]  + IR1[0]  + IR2[0])  / (n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow3[m+4+20]  * (I1[2]   + I2[2]  + IR1[2]  + IR2[2])  / (n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow3[m+8+20]  * (I1[4]   + I2[4]  + IR1[4]  + IR2[4])  / (n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow3[m+12+20] * (I1[6]   + I2[6]  + IR1[6]  + IR2[6])  / (n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow3[m+16+20] * (I1[8]   + I2[8]  + IR1[8]  + IR2[8])  / (n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow3[m+20+20] * (I1[10]  + I2[10] + IR1[10] + IR2[10]) / (n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow3[m+24+20] * (I1[12]  + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow3[m+28+20] * (I1[14]  + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow3[m+32+20] * (I1[16]  + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow3[m+36+20] * (I1[18]  + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh3[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh3[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh3[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh3[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh3[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh3[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh3[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh3[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh3[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh3[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				// population aging and movement in/out of sexually active group
				arma::vec agingS = mv_mult(aging_rate,S);
				arma::vec agingE = mv_mult(aging_rate,E);
				arma::vec agingI1 = mv_mult(aging_rate,I1);
				arma::vec agingI2 = mv_mult(aging_rate,I2);
				arma::vec agingL1 = mv_mult(aging_rate,L1);
				arma::vec agingL2 = mv_mult(aging_rate,L2);
				arma::vec agingT1 = mv_mult(aging_rate,T1);
				arma::vec agingT2 = mv_mult(aging_rate,T2);
				arma::vec agingT3 = mv_mult(aging_rate,T3);
				arma::vec agingSR = mv_mult(aging_rate,SR);
				arma::vec agingER = mv_mult(aging_rate,ER);
				arma::vec agingIR1 = mv_mult(aging_rate,IR1);
				arma::vec agingIR2 = mv_mult(aging_rate,IR2);
				arma::vec agingLR1 = mv_mult(aging_rate,LR1);
				arma::vec agingLR2 = mv_mult(aging_rate,LR2);
				arma::mat agingNSAmat = mat_mult(aging_nsa, aging_rate);
				arma::vec agingNSA = mv_mult(agingNSAmat, NSA); //aging from NSA to NSA 
				arma::mat agingSNSAmat = aging_rate - mat_mult(aging_nsa, aging_rate);
				arma::vec agingSNSA = mv_mult(agingSNSAmat, NSA); //aging from NSA to SA
				
				arma::mat birthsSmat = mat_mult(birth_rate, births_sa);
				arma::vec Ntot = S+E+I1+I2+L1+L2+T1+T2+T3+SR+ER+IR1+IR2+LR1+LR2+NSA;
				arma::vec birthsS = mv_mult(birthsSmat,Ntot); //births into S
				arma::mat birthsNSAmat = mat_mult(birth_rate, births_nsa);
				arma::vec birthsNSA = mv_mult(birthsNSAmat,Ntot); //births into NSA

				
				// Initialize the equations and integration
				for(int n=0; n<dim; n++) {   
					Pop[n] += (-lambda[n]*S[n] + return_rate_for_resusceptible*SR[n] + abx*E[n] + (1-rep_on)*(lambda1*T1[n] + lambda2*T2[n] + lambda3*T3[n]) + agingS[n] + agingSNSA[n] + birthsS[n]) * dt ;  // dS/dt
					Pop[n+dim] += (lambda[n]*S[n] - delta*E[n] - abx*E[n]  + agingE[n]) * dt; // dE/dt
					Pop[n+dim*2] += ( delta*E[n]  - gamma1*I1[n] - p_trt1[n]*I1[n] -alpha[n]*I1[n] - abx*I1[n] + agingI1[n] ) * dt; //dI1/dt
					Pop[n+dim*3] += ( gamma1*I1[n] - gamma2*I2[n] - p_trt2[n]*I2[n] -alpha[n]*I2[n] - abx*I2[n] + agingI2[n]) * dt; //dI2/dt
					Pop[n+dim*4] += (gamma2*I2[n] - gamma3*L1[n] - p_trt3[n]*L1[n] -alpha[n]*L1[n] - abx*L1[n] + agingL1[n]) * dt; //dL1/dt
					Pop[n+dim*5] += (gamma3*L1[n] - p_trt4[n]*L2[n] -alpha[n]*L2[n] - abx*L2[n] + agingL2[n] ) * dt; //dL2/dt

					//dT1/dt
					Pop[n+dim*6] += (p_trt1[n]*(I1[n]+IR1[n]) + // treatment seeking for infected primary
							             p_trt2[n]*(I2[n]+IR2[n]) + // treatment seeking for infected secondary
							alpha[n]*(I1[n]+I2[n])+ // screening for first infection
							alpha_repeat[n]*(IR1[n]+IR2[n]) + // screening for repeat infection after diagnosis
							abx*(I1[n]+I2[n]+IR1[n]+IR2[n]) - // background antibiotics usage
							lambda1*T1[n] + agingT1[n])*dt; // treatment exit rate and aging out

					Pop[n+dim*7] += (p_trt3[n]*(L1[n] + LR1[n]) + alpha[n]*L1[n]+alpha_repeat[n]*LR1[n] + abx*(L1[n]+LR1[n]) - lambda2*T2[n] + agingT2[n])*dt; //dT2/dt
					Pop[n+dim*8] += (p_trt4[n]*(L2[n] + LR2[n]) + alpha[n]*L2[n]+alpha_repeat[n]*LR2[n] + abx*(L2[n]+LR2[n]) - lambda3*T3[n] + agingT3[n])*dt ; //dT3/dt
					Pop[n+dim*9] += (-lambda[n]*SR[n] - return_rate_for_resusceptible*SR[n] + abx*ER[n] + rep_on*(lambda1*T1[n] + lambda2*T2[n] + lambda3*T3[n])  + agingSR[n])*dt;  // dSR/dt
					Pop[n+dim*10] += (lambda[n]*SR[n] - delta*ER[n] - abx*ER[n] + agingER[n])*dt; // dER/dt
					Pop[n+dim*11] += (delta*ER[n] - gamma1*IR1[n] - p_trt1[n]*IR1[n] -alpha_repeat[n]*IR1[n] - abx*IR1[n]  + agingIR1[n])*dt; //dIR1/dt
					Pop[n+dim*12] += (gamma1*IR1[n] - gamma2*IR2[n] - p_trt2[n]*IR2[n] -alpha_repeat[n]*IR2[n] - abx*IR2[n]  + agingIR2[n]) * dt; //dIR2/dt
					Pop[n+dim*13] += (gamma2*IR2[n] - gamma3*LR1[n] - p_trt3[n]*LR1[n] -alpha_repeat[n]*LR1[n] - abx*LR1[n] + agingLR1[n]) * dt; //dLR1/dt
					Pop[n+dim*14] += (gamma3*LR1[n] - p_trt4[n]*LR2[n] -alpha_repeat[n]*LR2[n] - abx*LR2[n]  + agingLR2[n]) * dt; //dLR2/dt
					Pop[n+dim*15] += (lambda[n]*(S[n]+SR[n])) * dt; //dInc/dt
					Pop[n+dim*16] += (lambda[n]*SR[n]) * dt; //dIncr/dt
					Pop[n+dim*17] += (rep_s[n]*p_trt1[n]*(I1[n]+IR1[n]) + rep*(alpha[n]*I1[n] + alpha_repeat[n]*IR1[n])) * dt; //dD1/dt
					Pop[n+dim*18] += (rep_s[n]*p_trt2[n]*(I2[n]+IR2[n]) + rep*(alpha[n]*I2[n] + alpha_repeat[n]*IR2[n])) * dt; //dD2/dt
					Pop[n+dim*19] += (rep_s[n]*p_trt3[n]*(L1[n]+LR1[n]) + rep*(alpha[n]*L1[n] + alpha_repeat[n]*LR1[n])) * dt; //dD3/dt
					Pop[n+dim*20] += (agingNSA[n] + birthsNSA[n]) * dt; //dNSA/dt
					Pop[n+dim*21] += (rep_s[n]*p_trt4[n]*(L2[n]+LR2[n]) + rep*(alpha[n]*L2[n] + alpha_repeat[n]*LR2[n])) * dt; //dD4/dt

					Pop[n+dim*22] += (rep_s[n]*p_trt1[n]*IR1[n] + rep*alpha_repeat[n]*
							IR1[n] + rep_s[n]*p_trt2[n]*IR2[n] + rep*alpha_repeat[n]*IR2[n] +
							rep_s[n]*p_trt3[n]*LR1[n] + rep*alpha_repeat[n]* LR1[n] ) * dt; //dDR/dt

					Pop[n+dim*23] += (alpha[n]*(S[n] + E[n] + I1[n] + I2[n] + L1[n] + L2[n]) + 
							alpha_repeat[n]*(SR[n] + ER[n] + IR1[n] + IR2[n] + LR1[n] + LR2[n])) * dt; //dCumulativeTesting/dt
				}
				}
					
			if (output_weekly) {
					// fill results table every week 
					outputs(y*52 + w, 0) = (float) y + ((float) w+1.0)/52.0; //year
					for(int n=0; n<dim*24; n++) {
						outputs(y*52 + w, 1+n) = Pop[n]; // 
						}
			} else {
					// fill results table with annual outputs
				if(w==0) {
					outputs(y,0) = y; //year
					for(int n=0; n<dim*24; n++) {
						outputs(y,1+n) = Pop[n]; // 
						}
					}
			}

    } /// end of week loop
  } /// end of year loop
  
  /* for(int m=0; m<nYrs; m++) { */
  /*   for(int n=0; n<dim*24+1; n++) { */
  /*     out(m,n) = outputs(m,n); */
  /*   } */
  /* } */
  
  return Rcpp::List::create(
    Rcpp::Named("initPop") = initPop,
    Rcpp::Named("out") = outputs
  );
}





// ###########################################
//   contact tracing counterfactual module
// ###########################################


//[[Rcpp::export]]  
List syphSimCTCounterfactual(
    List  x, //params
    double dt, // timestep
    List cm, //mixing 
    NumericVector pabx, //p.abx
    NumericVector rep_count, // when to turn counting of prior infections on 
    int nYrs, //model run time
    NumericVector initPop, //initial population
    NumericVector n_sa, //size of sexually active population
    arma::mat births, //birth rate
    arma::mat births_sa, //births p_sexually active
    arma::mat births_nsa, //births p_not sexually active
    arma::mat aging, //aging rate
    arma::mat aging_nsa //aging p.not sexually active
    
) {
  int i = 5; //number of subpopulations
  int j = 2; //number of sexes
  int k = 2; //number of sexual activity groups
  int l = 2; //number of age groups
  int dim = i*j*k*l; //total number of population subgroups in output matrix
  NumericVector b = as<NumericVector>(x["b"]); //transmission rate
  double delta = as<double>(x["delta"]); //1/dur incubation 
  NumericVector gamma = as<NumericVector>(x["gamma"]); //1/dur infectious
  NumericVector trt1 = as<NumericVector>(x["p.trt.1"]); //treatment rate, primary syph
  NumericVector trt2 = as<NumericVector>(x["p.trt.2"]); //treatment rate, secondary syph
  NumericVector trt3 = as<NumericVector>(x["p.trt.3"]); //treatment rate, early latent syph
  NumericVector trt4 = as<NumericVector>(x["p.trt.4"]); //treatment rate, late latent syph
  NumericVector report = as<NumericVector>(x["rep"]); //reporting probability, screened case
  NumericMatrix report_symp = as<NumericMatrix>(x["rep.symp"]); //reporting probability matrix, case seeking care
  NumericMatrix screen = as<NumericMatrix>(x["screen"]); //screening rate matrix
	NumericMatrix screen_ps = as<NumericMatrix>(x["screen_ps"]); // screening rate matrix for primary/secondary infections
	NumericMatrix screen_el = as<NumericMatrix>(x["screen_el"]); // screening rate matrix for primary/secondary infections
  NumericVector alpha(dim); //screening rate (updated annually)
	NumericVector alpha_ps(dim); //screening rate for primary/secondary
	NumericVector alpha_el(dim); //screening rate for early/late latent 
  NumericVector rep_s(dim); //reporting probabilities (updataed annually)
  NumericVector c_msm = as<NumericVector>(x["behav"]); //transmission rr in MSM 
  NumericVector dur_imm = as<NumericVector>(x["dur.imm"]); //duration immunity after treatment
  NumericVector cmlow1 = as<NumericVector>(cm["cm.low1"]); //contact matrix low activity, black 
  NumericVector cmlow2 = as<NumericVector>(cm["cm.low2"]); //contact matrix low activity, other 
  NumericVector cmlow3 = as<NumericVector>(cm["cm.low3"]); //contact matrix low activity, Hispanic 
  NumericVector cmlow4 = as<NumericVector>(cm["cm.low4"]); //contact matrix low activity, HIV- MSM
  NumericVector cmlow5 = as<NumericVector>(cm["cm.low5"]); //contact matrix low activity, HIV+ MSM
  NumericVector cmhigh1 = as<NumericVector>(cm["cm.high1"]); //contact matrix high activity, black 
  NumericVector cmhigh2 = as<NumericVector>(cm["cm.high2"]); //contact matrix high activity, other 
  NumericVector cmhigh3 = as<NumericVector>(cm["cm.high3"]); //contact matrix high activity, Hispanic 
  NumericVector cmhigh4 = as<NumericVector>(cm["cm.high4"]); //contact matrix high activity, HIV- MSM
  NumericVector cmhigh5 = as<NumericVector>(cm["cm.high5"]); //contact matrix high activity, HIV+ MSM
	bool output_weekly = as<LogicalVector>(x["output_weekly"])[0] == TRUE; // output every timestep or only yearly (default)
  arma::mat birth_rate = births; //birth rate
  arma::mat aging_rate = aging; //aging rate
  arma::vec S(dim);  //susceptible
  arma::vec E(dim);  //incubating
  arma::vec I1(dim); //primary syphilis
  arma::vec I2(dim); //secondary syphilis
  arma::vec L1(dim); //early latent syphilis
  arma::vec L2(dim); //late latent syphilis
  arma::vec T1(dim); //treated, primary and secondary
  arma::vec T2(dim); //treated, early latent
  arma::vec T3(dim); //treated, late latent
  arma::vec SR(dim); //susceptible, prior treated infection
  arma::vec ER(dim); //incubating, prior treated infection
  arma::vec IR1(dim); //primary, prior treated infection
  arma::vec IR2(dim); //secondary, prior treated infection
  arma::vec LR1(dim); //early latent, prior treated infection
  arma::vec LR2(dim); //late latent, prior treated infection
  arma::vec Inc(dim); //cumulative incidence
  arma::vec Incr(dim); //cumulative incidence, prior infection
  arma::vec D1(dim);  //reported cases, primary 
  arma::vec D2(dim);  //reported cases, secondary
  arma::vec D3(dim);  //reported cases, early latent
  arma::vec NSA(dim); //not sexually active 
  arma::vec D4(dim);  //reported cases, late latent
  arma::vec DR(dim);  //reported cases, primary, secondary & early latent, prior infection
  arma::vec p_trt1(dim); //trt rate, primary
  arma::vec p_trt2(dim); //trt rate, secondary
  arma::vec p_trt3(dim); //trt rate, early latent
  arma::vec p_trt4(dim); //trt rate, late latent
  arma::vec lambda(dim); //force of infection
  arma::vec Pop(dim*24); //population matrix
	
	int output_length; 
	if (output_weekly) {
		output_length = nYrs * 52;
	} else {
		output_length = nYrs;
	}
  NumericMatrix outputs(output_length, dim*24+1);
  NumericMatrix out(output_length, dim*24+1);
  
  double gamma1 = gamma[0]; //1/dur primary
  double gamma2 = gamma[1]; //1/dur secondary
  double gamma3 = gamma[2]; //1/dur latent
  double lambda1 = 365/dur_imm[0]; //1/dur imm primary and secondary
  double lambda2 = 365/dur_imm[1]; //1/dur imm early latent
  double lambda3 = 365/dur_imm[2]; //1/dur imm late latent
  
  //update treatment parameters for males
  for(int n=0; n<((i-2)*k*l); n++) {
    p_trt1[n] = trt1[0];
    p_trt2[n] = trt2[0];
    p_trt3[n] = trt3[0];
    p_trt4[n] = trt4[0];
  }
  //update treatment parameters for HIV- MSM
  for(int n=((i-2)*k*l); n<((i-1)*k*l); n++) {
    p_trt1[n] = trt1[1];
    p_trt2[n] = trt2[1];
    p_trt3[n] = trt3[1];
    p_trt4[n] = trt4[1];
  }
  //update treatment parameters for HIV+ MSM
  for(int n=((i-1)*k*l); n<i*k*l; n++) {
    p_trt1[n] = trt1[2];
    p_trt2[n] = trt2[2];
    p_trt3[n] = trt3[2];
    p_trt4[n] = trt4[2];
  }
  //update treatment parameters for females
  for(int n=((i*k*l)); n<(2*i*k*l); n++) {
    p_trt1[n] = trt1[3];
    p_trt2[n] = trt2[3];
    p_trt3[n] = trt3[3];
    p_trt4[n] = trt4[3];
  }
  
  /// year loop
  for(int y=0; y<nYrs; y++){
    for(int w=0; w<52; w++){

			/* for (int d=0; d<7; d++){ */
      
				//updata reporting (in screened), transmission rr, and background abx use
				double c_rr = c_msm[y];
				double rep = report[y];
				double abx = pabx[y];
				double rep_on = rep_count[y];
				
				if(y==0 && w==0){   ///populate the model at t=0
					for(int n=0; n<dim*24; n++){
						Pop[n] = initPop[n];
					}
				}
				//int s = y*52+w; //number of weeks since start
				else{
				for(int n=0; n<dim; n++) { //calculate population size at start of each week
					S[n] =   Pop[n];
					E[n] =   Pop[n+dim];
					I1[n] =  Pop[n+2*dim];
					I2[n] =  Pop[n+3*dim];
					L1[n] =  Pop[n+4*dim];
					L2[n] =  Pop[n+5*dim];
					T1[n] =  Pop[n+6*dim];
					T2[n] =  Pop[n+7*dim];
					T3[n] =  Pop[n+8*dim];
					SR[n] =  Pop[n+9*dim];
					ER[n] =  Pop[n+10*dim];
					IR1[n] = Pop[n+11*dim];
					IR2[n] = Pop[n+12*dim];
					LR1[n] = Pop[n+13*dim];
					LR2[n] = Pop[n+14*dim];
					Inc[n] = Pop[n+15*dim];
					Incr[n] = Pop[n+16*dim];
					D1[n] =  Pop[n+17*dim];
					D2[n] =  Pop[n+18*dim];
					D3[n] =  Pop[n+19*dim];
					NSA[n] = Pop[n+20*dim];
					D4[n] =  Pop[n+21*dim];
					DR[n] =  Pop[n+22*dim];
					alpha[n] = screen(y,n);  //update screening and reporting
					alpha_ps[n] = screen_ps(y,n);
					alpha_el[n] = screen_el(y,n);
					rep_s[n] = report_symp(y,n);
				}
				
				//calculate force of infection for each sex and subpopulation
				for(int n=0; n<4; n++) { //males subpop 1
					int m = n;
					lambda[n] =  (b[0] * (cmlow1[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow1[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow1[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow1[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow1[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow1[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh1[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) +  //F, high AC, age=1, subpop1
																cmhigh1[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) +  //F, high AC, age=2, subpop1
																cmhigh1[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh1[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh1[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh1[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=4; n<8; n++) {   //males subpop 2
					int m = n-4;
					lambda[n] =  (b[0] * (cmlow2[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow2[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow2[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow2[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow2[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow2[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh2[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh2[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh2[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh2[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh2[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh2[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=8; n<12; n++) {   //males subpop 3
					int m = n-8;
					lambda[n] =  (b[0] * (cmlow3[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow3[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow3[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow3[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow3[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow3[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh3[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh3[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh3[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh3[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh3[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh3[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))); //F, high AC, age=2, subpop3
				}
				
				for(int n=12; n<16; n++) {   //males subpop 4
					int m = n-12;
					lambda[n] =  ((b[0] *  (cmlow4[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																	cmlow4[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																	cmlow4[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																	cmlow4[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																	cmlow4[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																	cmlow4[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																	cmhigh4[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																	cmhigh4[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																	cmhigh4[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																	cmhigh4[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																	cmhigh4[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																	cmhigh4[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))) +//F, high AC, age=2, subpop3
								
											((b[2]*c_rr)*
																 (cmlow4[m+24]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																	cmlow4[m+28]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																	cmlow4[m+32]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																	cmlow4[m+36]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																	cmhigh4[m+24] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																	cmhigh4[m+28] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																	cmhigh4[m+32] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																	cmhigh4[m+36] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19])))); //M, high AC, age=2, subpop5
										
						}  
				
				for(int n=16; n<20; n++) {   //males subpop 5
					int m = n-16;
					lambda[n] =  ((b[0] *(cmlow5[m]     * (I1[20] + I2[20] + IR1[20] + IR2[20]) / (n_sa[20]) + //F, low AC, age=1, subpop1
																cmlow5[m+4]   * (I1[22] + I2[22] + IR1[22] + IR2[22]) / (n_sa[22]) + //F, low AC, age=2, subpop1
																cmlow5[m+8]   * (I1[24] + I2[24] + IR1[24] + IR2[24]) / (n_sa[24]) + //F, low AC, age=1, subpop2
																cmlow5[m+12]  * (I1[26] + I2[26] + IR1[26] + IR2[26]) / (n_sa[26]) + //F, low AC, age=2, subpop2
																cmlow5[m+16]  * (I1[28] + I2[28] + IR1[28] + IR2[28]) / (n_sa[28]) + //F, low AC, age=1, subpop3
																cmlow5[m+20]  * (I1[30] + I2[30] + IR1[30] + IR2[30]) / (n_sa[30]) + //F, low AC, age=2, subpop3
																cmhigh5[m]    * (I1[21] + I2[21] + IR1[21] + IR2[21]) / (n_sa[21]) + //F, high AC, age=1, subpop1
																cmhigh5[m+4]  * (I1[23] + I2[23] + IR1[23] + IR2[23]) / (n_sa[23]) + //F, high AC, age=2, subpop1
																cmhigh5[m+8]  * (I1[25] + I2[25] + IR1[25] + IR2[25]) / (n_sa[25]) + //F, high AC, age=1, subpop2
																cmhigh5[m+12] * (I1[27] + I2[27] + IR1[27] + IR2[27]) / (n_sa[27]) + //F, high AC, age=2, subpop2
																cmhigh5[m+16] * (I1[29] + I2[29] + IR1[29] + IR2[29]) / (n_sa[29]) + //F, high AC, age=1, subpop3
																cmhigh5[m+20] * (I1[31] + I2[31] + IR1[31] + IR2[31]) / (n_sa[31]))) + //F, high AC, age=2, subpop3
																
													((b[2]*c_rr)*
															 (cmlow5[m+24]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow5[m+28]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow5[m+32]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow5[m+36]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																cmhigh5[m+24] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh5[m+28] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh5[m+32] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh5[m+36] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19])))); //M, high AC, age=2, subpop5
															
				}
				
				for(int n=20; n<24; n++) {   //females subpop 1
					int m = n;
					lambda[n] =  (b[1] * (cmlow1[m+20]     * (I1[0]  + I2[0]  + IR1[0]  + IR2[0])  /(n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow1[m+4+20]   * (I1[2]  + I2[2]  + IR1[2]  + IR2[2])  /(n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow1[m+8+20]   * (I1[4]  + I2[4]  + IR1[4]  + IR2[4])  /(n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow1[m+12+20]  * (I1[6]  + I2[6]  + IR1[6]  + IR2[6])  /(n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow1[m+16+20]  * (I1[8]  + I2[8]  + IR1[8]  + IR2[8])  /(n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow1[m+20+20]  * (I1[10] + I2[10] + IR1[10] + IR2[10]) /(n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow1[m+24+20]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) /(n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow1[m+28+20]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) /(n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow1[m+32+20]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) /(n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow1[m+36+20]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) /(n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh1[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh1[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh1[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh1[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh1[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh1[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh1[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh1[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh1[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh1[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				for(int n=24; n<28; n++) {   //females subpop 2
					int m = n-4;
					lambda[n] =  (b[1] * (cmlow2[m+20]     * (I1[0]  + I2[0]  + IR1[0]  + IR2[0])  / (n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow2[m+4+20]   * (I1[2]  + I2[2]  + IR1[2]  + IR2[2])  / (n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow2[m+8+20]   * (I1[4]  + I2[4]  + IR1[4]  + IR2[4])  / (n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow2[m+12+20]  * (I1[6]  + I2[6]  + IR1[6]  + IR2[6])  / (n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow2[m+16+20]  * (I1[8]  + I2[8]  + IR1[8]  + IR2[8])  / (n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow2[m+20+20]  * (I1[10] + I2[10] + IR1[10] + IR2[10]) / (n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow2[m+24+20]  * (I1[12] + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow2[m+28+20]  * (I1[14] + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow2[m+32+20]  * (I1[16] + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow2[m+36+20]  * (I1[18] + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh2[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh2[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh2[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh2[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh2[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh2[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh2[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh2[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh2[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh2[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				for(int n=28; n<32; n++) {   //females subpop 3
					int m = n-8;
					lambda[n] =  (b[1] * (cmlow3[m+20]    * (I1[0]   + I2[0]  + IR1[0]  + IR2[0])  / (n_sa[0])  + //M, low AC, age=1, subpop1
																cmlow3[m+4+20]  * (I1[2]   + I2[2]  + IR1[2]  + IR2[2])  / (n_sa[2])  + //M, low AC, age=2, subpop1
																cmlow3[m+8+20]  * (I1[4]   + I2[4]  + IR1[4]  + IR2[4])  / (n_sa[4])  + //M, low AC, age=1, subpop2
																cmlow3[m+12+20] * (I1[6]   + I2[6]  + IR1[6]  + IR2[6])  / (n_sa[6])  + //M, low AC, age=2, subpop2
																cmlow3[m+16+20] * (I1[8]   + I2[8]  + IR1[8]  + IR2[8])  / (n_sa[8])  + //M, low AC, age=1, subpop3
																cmlow3[m+20+20] * (I1[10]  + I2[10] + IR1[10] + IR2[10]) / (n_sa[10]) + //M, low AC, age=2, subpop3
																cmlow3[m+24+20] * (I1[12]  + I2[12] + IR1[12] + IR2[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
																cmlow3[m+28+20] * (I1[14]  + I2[14] + IR1[14] + IR2[14]) / (n_sa[14]) + //M, low AC, age=2, subpop4
																cmlow3[m+32+20] * (I1[16]  + I2[16] + IR1[16] + IR2[16]) / (n_sa[16]) + //M, low AC, age=1, subpop5
																cmlow3[m+36+20] * (I1[18]  + I2[18] + IR1[18] + IR2[18]) / (n_sa[18]) + //M, low AC, age=2, subpop5
																
																cmhigh3[m+20]    * (I1[1]  + I2[1]  + IR1[1]  + IR2[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
																cmhigh3[m+4+20]  * (I1[3]  + I2[3]  + IR1[3]  + IR2[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
																cmhigh3[m+8+20]  * (I1[5]  + I2[5]  + IR1[5]  + IR2[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
																cmhigh3[m+12+20] * (I1[7]  + I2[7]  + IR1[7]  + IR2[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
																cmhigh3[m+16+20] * (I1[9]  + I2[9]  + IR1[9]  + IR2[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
																cmhigh3[m+20+20] * (I1[11] + I2[11] + IR1[11] + IR2[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
																cmhigh3[m+24+20] * (I1[13] + I2[13] + IR1[13] + IR2[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
																cmhigh3[m+28+20] * (I1[15] + I2[15] + IR1[15] + IR2[15]) / (n_sa[15]) + //M, high AC, age=2, subpop4
																cmhigh3[m+32+20] * (I1[17] + I2[17] + IR1[17] + IR2[17]) / (n_sa[17]) + //M, high AC, age=1, subpop5
																cmhigh3[m+36+20] * (I1[19] + I2[19] + IR1[19] + IR2[19]) / (n_sa[19]))); //M, high AC, age=2, subpop5
					
				}
				
				// population aging and movement in/out of sexually active group
				arma::vec agingS = mv_mult(aging_rate,S);
				arma::vec agingE = mv_mult(aging_rate,E);
				arma::vec agingI1 = mv_mult(aging_rate,I1);
				arma::vec agingI2 = mv_mult(aging_rate,I2);
				arma::vec agingL1 = mv_mult(aging_rate,L1);
				arma::vec agingL2 = mv_mult(aging_rate,L2);
				arma::vec agingT1 = mv_mult(aging_rate,T1);
				arma::vec agingT2 = mv_mult(aging_rate,T2);
				arma::vec agingT3 = mv_mult(aging_rate,T3);
				arma::vec agingSR = mv_mult(aging_rate,SR);
				arma::vec agingER = mv_mult(aging_rate,ER);
				arma::vec agingIR1 = mv_mult(aging_rate,IR1);
				arma::vec agingIR2 = mv_mult(aging_rate,IR2);
				arma::vec agingLR1 = mv_mult(aging_rate,LR1);
				arma::vec agingLR2 = mv_mult(aging_rate,LR2);
				arma::mat agingNSAmat = mat_mult(aging_nsa, aging_rate);
				arma::vec agingNSA = mv_mult(agingNSAmat, NSA); //aging from NSA to NSA 
				arma::mat agingSNSAmat = aging_rate - mat_mult(aging_nsa, aging_rate);
				arma::vec agingSNSA = mv_mult(agingSNSAmat, NSA); //aging from NSA to SA
				
				arma::mat birthsSmat = mat_mult(birth_rate, births_sa);
				arma::vec Ntot = S+E+I1+I2+L1+L2+T1+T2+T3+SR+ER+IR1+IR2+LR1+LR2+NSA;
				arma::vec birthsS = mv_mult(birthsSmat,Ntot); //births into S
				arma::mat birthsNSAmat = mat_mult(birth_rate, births_nsa);
				arma::vec birthsNSA = mv_mult(birthsNSAmat,Ntot); //births into NSA

				
				// Initialize the equations and integration
				for(int n=0; n<dim; n++) {   
					Pop[n] += (-lambda[n]*S[n] + abx*E[n] + (1-rep_on)*(lambda1*T1[n] + lambda2*T2[n] + lambda3*T3[n]) + agingS[n] + agingSNSA[n] + birthsS[n]) * dt ;  // dS/dt
					Pop[n+dim] += (lambda[n]*S[n] - delta*E[n] - abx*E[n]  + agingE[n]) * dt; // dE/dt
					Pop[n+dim*2] += ( delta*E[n]  - gamma1*I1[n] - p_trt1[n]*I1[n] -alpha_ps[n]*I1[n] - abx*I1[n] + agingI1[n] ) * dt; //dI1/dt
					Pop[n+dim*3] += ( gamma1*I1[n] - gamma2*I2[n] - p_trt2[n]*I2[n] -alpha_ps[n]*I2[n] - abx*I2[n] + agingI2[n]) * dt; //dI2/dt
					Pop[n+dim*4] += (gamma2*I2[n] - gamma3*L1[n] - p_trt3[n]*L1[n] -alpha_el[n]*L1[n] - abx*L1[n] + agingL1[n]) * dt; //dL1/dt
					Pop[n+dim*5] += (gamma3*L1[n] - p_trt4[n]*L2[n] -alpha_el[n]*L2[n] - abx*L2[n] + agingL2[n] ) * dt; //dL2/dt
					Pop[n+dim*6] += (p_trt1[n]*(I1[n]+IR1[n]) + p_trt2[n]*(I2[n] + IR2[n]) + alpha_ps[n]*(I1[n]+I2[n]+IR1[n]+IR2[n]) + abx*(I1[n]+I2[n]+IR1[n]+IR2[n]) - lambda1*T1[n] + agingT1[n])*dt ; //dT1/dt
					Pop[n+dim*7] += (p_trt3[n]*(L1[n] + LR1[n]) + alpha_el[n]*(L1[n]+LR1[n]) + abx*(L1[n]+LR1[n]) - lambda2*T2[n] + agingT2[n])*dt; //dT2/dt
					Pop[n+dim*8] += (p_trt4[n]*(L2[n] + LR2[n]) + alpha_el[n]*(L2[n]+LR2[n]) + abx*(L2[n]+LR2[n]) - lambda3*T3[n] + agingT3[n])*dt ; //dT3/dt
					Pop[n+dim*9] += (-lambda[n]*SR[n] + abx*ER[n] + rep_on*(lambda1*T1[n] + lambda2*T2[n] + lambda3*T3[n])  + agingSR[n])*dt;  // dSR/dt
					Pop[n+dim*10] += (lambda[n]*SR[n] - delta*ER[n] - abx*ER[n] + agingER[n])*dt; // dER/dt
					Pop[n+dim*11] += (delta*ER[n] - gamma1*IR1[n] - p_trt1[n]*IR1[n] -alpha_ps[n]*IR1[n] - abx*IR1[n]  + agingIR1[n])*dt; //dIR1/dt
					Pop[n+dim*12] += (gamma1*IR1[n] - gamma2*IR2[n] - p_trt2[n]*IR2[n] -alpha_ps[n]*IR2[n] - abx*IR2[n]  + agingIR2[n]) * dt; //dIR2/dt
					Pop[n+dim*13] += (gamma2*IR2[n] - gamma3*LR1[n] - p_trt3[n]*LR1[n] -alpha_el[n]*LR1[n] - abx*LR1[n] + agingLR1[n]) * dt; //dLR1/dt
					Pop[n+dim*14] += (gamma3*LR1[n] - p_trt4[n]*LR2[n] -alpha_el[n]*LR2[n] - abx*LR2[n]  + agingLR2[n]) * dt; //dLR2/dt
					Pop[n+dim*15] += (lambda[n]*(S[n]+SR[n])) * dt; //dInc/dt
					Pop[n+dim*16] += (lambda[n]*SR[n]) * dt; //dIncr/dt
					Pop[n+dim*17] += (rep_s[n]*p_trt1[n]*(I1[n]+IR1[n]) + rep*alpha_ps[n]*(I1[n] + IR1[n])) * dt; //dD1/dt
					Pop[n+dim*18] += (rep_s[n]*p_trt2[n]*(I2[n]+IR2[n]) + rep*alpha_ps[n]*(I2[n] + IR2[n])) * dt; //dD2/dt
					Pop[n+dim*19] += (rep_s[n]*p_trt3[n]*(L1[n]+LR1[n]) + rep*alpha_el[n]*(L1[n] + LR1[n])) * dt; //dD3/dt
					Pop[n+dim*20] += (agingNSA[n] + birthsNSA[n]) * dt; //dNSA/dt
					Pop[n+dim*21] += (rep_s[n]*p_trt4[n]*(L2[n]+LR2[n]) + rep*alpha_el[n]*(L2[n] + LR2[n])) * dt; //dD4/dt
					Pop[n+dim*22] += (rep_s[n]*p_trt1[n]*IR1[n] + rep*alpha_ps[n]* IR1[n] + rep_s[n]*p_trt2[n]*IR2[n] + rep*alpha_ps[n]*IR2[n] + rep_s[n]*p_trt3[n]*LR1[n] + rep*alpha_el[n]* LR1[n] ) * dt; //dDR/dt
					Pop[n+dim*23] += (alpha[n]*(S[n] + E[n] + SR[n] + ER[n]) + alpha_ps[n] * (I1[n] + IR1[n] + I2[n] + IR2[n]) + alpha_el[n] * (L1[n] + LR1[n] + L2[n] + LR2[n])) * dt; //dCumulativeTesting/dt
				}
				}
					
			if (output_weekly) {
					// fill results table every week 
					outputs(y*52 + w, 0) = (float) y + ((float) w+1.0)/52.0; //year
					for(int n=0; n<dim*24; n++) {
						outputs(y*52 + w, 1+n) = Pop[n]; // 
						}
			} else {
					// fill results table with annual outputs
				if(w==0) {
					outputs(y,0) = y; //year
					for(int n=0; n<dim*24; n++) {
						outputs(y,1+n) = Pop[n]; // 
						}
					}
			}

    } /// end of week loop
  } /// end of year loop
  
  /* for(int m=0; m<nYrs; m++) { */
  /*   for(int n=0; n<dim*24+1; n++) { */
  /*     out(m,n) = outputs(m,n); */
  /*   } */
  /* } */
  
  return Rcpp::List::create(
    Rcpp::Named("initPop") = initPop,
    Rcpp::Named("out") = outputs
  );
}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// mv_mult(aging,rep(1,40))
// mat_mult(aging.nsa,aging)
// syphSim(params,
//               ts,
//               cm.list,
//               p.abx,
//               5, ##model.end
//               yinit, 
//               n.sa,
//               births,
//               births.sa,
//               births.nsa,
//               aging,
//               aging.nsa)
// */

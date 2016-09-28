// *********************************************************************** //
// MCMC code implementing integrative Bayesian variable selection for the  //
// Dirichlet-Multinomial regression model with application to the analysis //
// of taxonomic abundances in microbiome data                              //
// *********************************************************************** //
//
// Copyright (C) 2016, Raffaele Argiento and W. Duncan Wadsworth
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//    
// *********************************************************************** //
// standard includes
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// GSL includes
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h> // has both the lngamma and lnbeta functions
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

// define the RNG globally
gsl_rng *rando;

// data dimensions required in the function prototypes below
static int n_cats;
static int n_obs;
static int n_vars;

// ----- Function declarations ----- //
// adaptive MH helper functions
double online_mean(int iteration, double last_mean, double curr_obs);
double online_var(int iteration, double last_mean, double last_var, 
                  double curr_mean, double curr_obs);
double adap_prop(double curr_var);
// log Beta-Bernoulli spike-and-slab prior evaluation
double lprior_bbsas(double betajk, int sjk, double sig_bejk, double mu_bejk, 
                    double aa_hp, double bb_hp);
// calculates the linear predictor for each category of the Dirichlet-Multinomial
double calculate_gamma(double **XX, double *alpha, double *beta, int jj, 
                       int ii, int Log);
// MH update for the regression parameters
void update_beta_jj(double **XX, double **JJ, double **loggamma, 
                    double *beta_temp, int *inclusion_indicator, 
                    double *prop_per_beta, double mu_be[n_cats][n_vars], 
                    double sig_be[n_cats][n_vars], double aa_hp, double bb_hp, 
                    int jj);
// MH update for the intercept parameters
void update_alpha_jj(double **JJ, double **loggamma, double *alpha, 
                     double *prop_per_alpha, int *accepted_alpha_flag, 
                     double mu_al[n_cats], double sig_al[n_cats], int jj);
// for the Savitsky et al. inclusion proposal step
void between_models_jj(double **XX, double **JJ, double **loggamma, 
                       double *beta_temp, int *accepted_beta_flag, 
                       int *inclusion_indicator, double mu_be[n_cats][n_vars], 
                       double sig_be[n_cats][n_vars], double aa_hp, double bb_hp, 
                       int jj);
// Principal function
void dmbvs_gibbs(double **XX, int **YY, double *alpha, double *beta, 
                 double mu_al[n_cats], double sig_al[n_cats], 
                 double mu_be[n_cats][n_vars], double sig_be[n_cats][n_vars], 
                 double aa_hp, double bb_hp, double *prop_per_alpha, 
                 double *prop_per_beta, const int GG, const int burn, 
                 const int thin, int Log, char *temp_dir);

// -------------------------------- //
int main(int argc, char *argv[]){
  
  // command line argument parsing
  if(argc != 13){
    fprintf(stderr, "Usage: %s <Iters> <Thin> <Burn> <Intercept Variance> <Slab Variance>\n <Beta-Bernoulli Alpha> <Beta-Bernoulli Beta> <N Categories>\n <N Observations> <N Variables> <External Seed> <Temp Directory>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  const int GG = atoi(argv[1]);
  const int thin = atoi(argv[2]);
  const int burn = atoi(argv[3]);
  double alpha_variance = atof(argv[4]);
  double slab_variance = atof(argv[5]);
  double aa_hp = atof(argv[6]);
  double bb_hp = atof(argv[7]);
  n_cats = atoi(argv[8]);
  n_obs = atoi(argv[9]);
  n_vars = atoi(argv[10]);
  long int external_seed = atoi(argv[11]);
  char *temp_dir = argv[12];
  
  // hyperparameters and Metropolis-Hastings proposal parameters
  const double mean_alpha = 0.0;
  const double mean_slab = 0.0;
  
  printf("***********************************************************\n");
  printf("********************* Starting dmbvs **********************\n");
  // outfile info
  printf("***********************************************************\n");
  printf("*************************** MCMC **************************\n");
  printf("Total iterations = %i, Thinning = %i, Burnin = %i\n", GG, thin, burn);
  // random number generator initialization
  const gsl_rng_type * Tipo;
  Tipo = gsl_rng_default;
  rando = gsl_rng_alloc(Tipo);
  gsl_rng_set(rando, external_seed);
  printf("External seed = %lu, GSL seed = %lu, Generator type: %s\n", 
         external_seed, gsl_rng_get(rando), gsl_rng_name(rando));
  printf("Text output deposited in directory:\n %s\n", temp_dir);
  printf("***********************************************************\n");
  printf("*************************** Data **************************\n");
  printf("# of categories = %i, # of observations = %i, # of covariates = %i\n", 
         n_cats, n_obs, n_vars);
  printf("***********************************************************\n");
  printf("********************* Hyperparameters *********************\n");
  printf("Slab variance = %0.1f, Intercept variance = %0.1f\n", slab_variance, 
         alpha_variance);
  printf("Beta-Bernoulli alpha = %0.2f, Beta-Bernoulli beta = %0.2f\n", aa_hp, 
         bb_hp);
  printf("***********************************************************\n");
  printf("***********************************************************\n");
  
  ///////////////////////////////////////////////////
  // Indices (shifted right one)
  // ii: index of the observations in 1, ..., n_obs
  // jj: index of the categories in 1, ..., n_cats
  // kk: index of the covariates in 1, ..., n_vars
  // hh: index of the beta vector in 1, ..., n_vars * n_cats, is collapsed 
  //     row-wise so that the 2x2 identity matrix becomes (1, 0, 0, 1)
  int ii, jj, kk, hh;
  
  // Larger objects such as data and initializations are read from files
  
  // single file pointer for all files
  FILE *fin;
  
  // reading in YY: the sample-by-category count matrix
  // all other files are read in in the same way
  char td_YY[200];
  strcpy(td_YY, temp_dir);
  strcat(td_YY, "/count_matrix.txt");
  fin = fopen(td_YY, "r");
  // dynamic allocation for YY for an array of pointers
  int **YY;
  YY = (int **) malloc(n_obs * sizeof(int *));
  for(ii = 0 ; ii < n_obs ; ii++){
    // allocate each row of the matrix YY
    YY[ii] = (int *) malloc(n_cats * sizeof(int));
    // fill in each row of the matrix YY
    for(jj = 0 ; jj < n_cats ; jj++){
      fscanf(fin, "%u", &YY[ii][jj]);
    }
  }
  
  // reading in XX: the sample-by-covariate matrix
  char td_XX[200];
  strcpy(td_XX, temp_dir);
  strcat(td_XX, "/covariates.txt");
  fin = fopen(td_XX, "r");
  double **XX;
  XX = (double **) malloc(n_obs * sizeof(double *));
  for(ii = 0 ; ii < n_obs ; ii++){
    XX[ii] = (double *) malloc(n_vars * sizeof(double));
    for(jj = 0 ; jj < n_vars ; jj++){
      fscanf(fin, "%lf", &XX[ii][jj]); 
    }
  }
  
  ///////// Initializations from text files
  
  // initialization for beta
  char td_ib[200];
  strcpy(td_ib, temp_dir);
  strcat(td_ib, "/init_beta.txt");
  fin = fopen(td_ib, "r");
  double *beta;
  beta = malloc((n_vars * n_cats) * sizeof(double));
  for(hh = 0 ; hh < (n_vars * n_cats) ; hh++){
    fscanf(fin, "%lf", &beta[hh]);
  }
  
  // initialization for alpha
  char td_ia[200];
  strcpy(td_ia, temp_dir);
  strcat(td_ia, "/init_alpha.txt");
  fin = fopen(td_ia, "r");
  double *alpha;
  alpha = malloc(n_cats * sizeof(double)); 
  for(jj = 0 ; jj < n_cats ; jj++){
    fscanf(fin, "%lf", &alpha[jj]);
  }
  
  // for the beta proposal variance
  char td_pb[200];
  strcpy(td_pb, temp_dir);
  strcat(td_pb, "/proposal_beta.txt");
  fin = fopen(td_pb, "r");
  double *prop_per_beta;
  prop_per_beta = malloc((n_vars * n_cats) * sizeof(double));
  for(hh = 0 ; hh < (n_vars * n_cats) ; hh++){
    fscanf(fin, "%lf", &prop_per_beta[hh]);
  }
  
  // for the alpha proposal variance
  char td_pa[200];
  strcpy(td_pa, temp_dir);
  strcat(td_pa, "/proposal_alpha.txt");
  fin = fopen(td_pa, "r");
  double *prop_per_alpha;
  prop_per_alpha = malloc(n_cats * sizeof(double));
  for(jj = 0 ; jj < n_cats ; jj++){
    fscanf(fin, "%lf", &prop_per_alpha[jj]); 
  }
  
  // mean and standard deviation of the independent normal priors on alpha and beta
  double mu_al[n_cats], sig_al[n_cats]; 
  double mu_be[n_cats][n_vars], sig_be[n_cats][n_vars];
  for(jj = 0 ; jj < n_cats ; jj++){
    mu_al[jj] = mean_alpha;
    sig_al[jj] = alpha_variance;
    for(kk = 0 ; kk < n_vars ; kk++){
      mu_be[jj][kk] = mean_slab;
      sig_be[jj][kk] = slab_variance;
    }
  }
  
  // a flag for functions that can calculate on the log scale
  int Log = 1;
  
  ///////// Initialization finished: commence the MCMC
  
  dmbvs_gibbs(XX, YY, alpha, beta, mu_al, sig_al, mu_be, sig_be, aa_hp, bb_hp, 
              prop_per_alpha, prop_per_beta, GG, burn, thin, Log, temp_dir);
  
  // cleanup (kind of pedantic)
  fclose(fin);
  free(YY);
  free(XX);
  free(beta);
  free(alpha);
  free(prop_per_alpha);
  free(prop_per_beta);
  gsl_rng_free(rando);
  
  return(0);
  
} // end of main()

//////////////////////////////////////
// ----- Function definitions ----- //
//////////////////////////////////////

// Gibbs Sampler function
void dmbvs_gibbs(double **XX, int **YY, double *alpha, double *beta, 
                 double mu_al[n_cats], double sig_al[n_cats], 
                 double mu_be[n_cats][n_vars], double sig_be[n_cats][n_vars], 
                 double aa_hp, double bb_hp, double *prop_per_alpha, 
                 double *prop_per_beta, const int GG, const int burn, 
                 const int thin, int Log, char *temp_dir){
  
  // loop indices
  int hh, ii, jj, kk, gg;
  
  // for timing iterations
  time_t mytime;
  
  // adaptive proposal values
  double last_mean;
  double last_var;
  
  ///////// Initialize variables that don't require input or can be obtained 
  ///////// from inputs
  
  // for MH acceptance/rejection ratios, initialized with zeros
  int *accepted_beta_flag;
  accepted_beta_flag = malloc((n_vars * n_cats) * sizeof(int));
  int *accepted_beta;
  accepted_beta = malloc((n_vars * n_cats) * sizeof(int));
  for(hh = 0 ; hh < (n_vars * n_cats) ; hh++){
    accepted_beta_flag[hh] = 0;
    accepted_beta[hh] = 0;
  }
  int *accepted_alpha_flag;
  accepted_alpha_flag = malloc(n_cats * sizeof(int));
  int *accepted_alpha;
  accepted_alpha = malloc(n_cats * sizeof(int));
  for(hh = 0 ; hh < n_cats ; hh++){
    accepted_alpha_flag[hh] = 0;
    accepted_alpha[hh] = 0;
  }
  // for adaptive MH need to keep track of target density mean and variance
  double *curr_mean;
  curr_mean = malloc((n_vars * n_cats) * sizeof(double));
  double *curr_var;
  curr_var = malloc((n_vars * n_cats) * sizeof(double));
  for(hh = 0 ; hh < (n_vars * n_cats) ; hh++){
    curr_mean[hh] = 0.0;
    curr_var[hh] = 0.5;
  }
  
  // variable inclusion indicator
  int *inclusion_indicator;
  inclusion_indicator = malloc((n_vars * n_cats) * sizeof(int));
  for(hh = 0 ; hh < (n_vars * n_cats) ; hh++){
    if(beta[hh] == 0){
      inclusion_indicator[hh] = 0;
    }else{
      inclusion_indicator[hh] = 1; 
    }
  }
  
  // row sums of YY
  int *Ypiu;
  Ypiu = malloc(n_obs * sizeof(int));
  for(ii = 0 ; ii < n_obs ; ii++){
    for(jj = 0 ; jj < n_cats ; jj++){
      Ypiu[ii] = Ypiu[ii] + YY[ii][jj];
    }
  }
  
  // the beta_temp vector is an n_vars long vector used for updating each of the
  // jj-th rows in Beta corresponding the jj-th category in YY
  double *beta_temp;
  beta_temp = malloc(n_vars * sizeof(double));
  
  // the linear predictor matrix, represented as a vector, is in the log scale
  double **loggamma;
  loggamma = malloc(n_obs * sizeof(double *));
  for(ii = 0 ; ii < n_obs ; ii++){
    loggamma[ii] = malloc(n_cats * sizeof(double));
    for(jj = 0 ; jj < n_cats ; jj++){
      loggamma[ii][jj] = calculate_gamma(XX, alpha, beta, jj, ii, Log);
    }
  }
  
  // initialize latent variable: JJ ~ gamma()
  double **JJ;
  JJ = malloc(n_obs * sizeof(double *));
  // TT is the vector of normalization constants
  double *TT;
  TT = malloc(n_obs * sizeof(double *));
  for(ii = 0 ; ii < n_obs ; ii++){
    JJ[ii] = malloc(n_cats * sizeof(double));
    TT[ii] = 0.0;
    for(jj = 0 ; jj < n_cats ; jj++){
      JJ[ii][jj] = (double)YY[ii][jj];
      // a very small number
      if(JJ[ii][jj] < pow(10.0, -100.0)){
        JJ[ii][jj] = pow(10.0, -100.0);
      }
      TT[ii] = TT[ii] + JJ[ii][jj];
    }
  }
  
  // initialize the other clever latent variable: uu
  // note that this uses the data to initialize
  double *uu;
  uu = malloc(n_obs * sizeof(double));
  for(ii = 0 ; ii < n_obs ; ii++){
    uu[ii] = gsl_ran_gamma(rando, Ypiu[ii], 1.0/TT[ii]);
  }
  
  // output file pointers
  FILE *fout_alpha, *fout_beta, *fout_alpha_acceptance, *fout_beta_acceptance, *fout_beta_proposal;
  
  // this is stupid, I'm sure there's a better way
  char tf_alpha[200];
  strcpy(tf_alpha, temp_dir);
  strcat(tf_alpha, "/alpha.out");
  fout_alpha = fopen(tf_alpha, "w");
  char tf_alpha_acceptance[200];
  strcpy(tf_alpha_acceptance, temp_dir);
  strcat(tf_alpha_acceptance, "/alpha_acceptance.out");
  fout_alpha_acceptance = fopen(tf_alpha_acceptance, "w");
  char tf_beta[200];
  strcpy(tf_beta, temp_dir);
  strcat(tf_beta, "/beta.out");
  fout_beta = fopen(tf_beta, "w");
  char tf_beta_acceptance[200];
  strcpy(tf_beta_acceptance, temp_dir);
  strcat(tf_beta_acceptance, "/beta_acceptance.out");
  fout_beta_acceptance = fopen(tf_beta_acceptance, "w");
  char tf_beta_proposal[200];
  strcpy(tf_beta_proposal, temp_dir);
  strcat(tf_beta_proposal, "/beta_proposal.out");
  fout_beta_proposal = fopen(tf_beta_proposal, "w");
  
  // the magic starts here
  printf("Starting the sampler\n");	
  
  for(gg = 0 ; gg < GG ; gg++){
    
    // timing info to the printed output
    mytime = time(NULL);
    if(gg % 1000 == 0){
      printf("Starting the %uth iteration out of %u at %s", gg, GG, ctime(&mytime));
    }
    
    // first a round of the between-model step for every covariate within every taxa
    for(jj = 0 ; jj < n_cats ; jj++){
      
      // fill in beta_temp
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta_temp[kk] = beta[hh];
      }
      
      between_models_jj(XX, JJ, loggamma, beta_temp, accepted_beta_flag, 
                        inclusion_indicator, mu_be, sig_be, aa_hp, bb_hp, jj);
      
      // update beta with beta_temp
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta[hh] = beta_temp[kk];
      }
    }
    
    // now an accelerating within-model step 
    for(jj = 0 ; jj < n_cats ; jj++){
      
      // fill in beta_temp and update the proposal variance
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta_temp[kk] = beta[hh];
        // wait until each parameter has had a few iterations
        if(gg > 2 * n_vars * n_cats){
          // calculate the online variance
          last_mean = curr_mean[hh];
          last_var = curr_var[hh];
          curr_mean[hh] = online_mean(gg, last_mean, beta[hh]);
          curr_var[hh] = online_var(gg, last_mean, last_var, curr_mean[hh], beta[hh]);
          // update proposal variance
          prop_per_beta[hh] = curr_var[hh];
        }
      }
      
      update_beta_jj(XX, JJ, loggamma, beta_temp, inclusion_indicator, 
                     prop_per_beta, mu_be, sig_be, aa_hp, bb_hp, jj);
      
      // update beta with beta_temp and write to file
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta[hh] = beta_temp[kk];
        if((gg >= burn) & (gg % thin == 0)){
          accepted_beta[hh] = accepted_beta[hh] + accepted_beta_flag[hh];
          accepted_beta_flag[hh] = 0;
          fprintf(fout_beta, "%e ", beta[hh]); // space delimited output files
        }
      }
    }
    if((gg >= burn) & (gg % thin == 0)){fprintf(fout_beta, "\n");}
    
    // update alpha and write to file
    for(jj = 0 ; jj < n_cats ; jj++){
      update_alpha_jj(JJ, loggamma, alpha, prop_per_alpha, accepted_alpha_flag, 
                      mu_al, sig_al, jj);
      if((gg >= burn) & (gg % thin == 0)){
        accepted_alpha[jj] = accepted_alpha[jj] + accepted_alpha_flag[jj];
        accepted_alpha_flag[jj] = 0;
        fprintf(fout_alpha,"%e ", alpha[jj]);
      }
    }
    if((gg >= burn) & (gg % thin == 0)){fprintf(fout_alpha, "\n ");}
    
    // update JJ and consequently TT 
    for(ii = 0 ; ii < n_obs ; ii++){
      TT[ii] = 0.0;
      for(jj = 0 ; jj < n_cats ; jj++){
        JJ[ii][jj] = gsl_ran_gamma(rando, YY[ii][jj] + exp(loggamma[ii][jj]), 1.0/(uu[ii] + 1.0));
        if(JJ[ii][jj] < pow(10.0, -100.0)){
          JJ[ii][jj] = pow(10.0, -100.0);
        }
        TT[ii] = TT[ii] + JJ[ii][jj];
      }
    }
    
    // update latent variables uu
    for(ii = 0 ; ii < n_obs ; ii++){
      uu[ii] = gsl_ran_gamma(rando, Ypiu[ii], 1.0/TT[ii]);
    }
    
  } // end of iterations
  
  // print acceptance ratio files
  for(jj = 0 ; jj < n_cats ; jj++){
    fprintf(fout_alpha_acceptance, "%f\n", (float)accepted_alpha[jj]/((GG - burn)/thin));
    for(kk = 0 ; kk < n_vars ; kk++){
      fprintf(fout_beta_acceptance, "%f\n", (float)accepted_beta[kk + jj * n_vars]/((GG - burn)/thin));
      fprintf(fout_beta_proposal, "%f\n", prop_per_beta[kk + jj * n_vars]);
    }
  }
  
  // cleanup (kind of pedantic)
  fclose(fout_alpha);
  fclose(fout_beta);
  fclose(fout_alpha_acceptance);
  fclose(fout_beta_acceptance);
  free(accepted_alpha);
  free(accepted_beta);
  free(accepted_alpha_flag);
  free(accepted_beta_flag);
  free(inclusion_indicator);
  free(beta_temp);
  free(loggamma);
  free(JJ);
  free(TT);
  free(uu);
  free(Ypiu);
    
  printf("Sampling finished!\n");  
  
  return;
  
} // close Gibbs sampler function

// 
double calculate_gamma(double **XX, double *alpha, double *beta, int jj, int ii, int Log){
  
  // loop indices and out
  int kk, hh;
  double out;
  
  out = alpha[jj];
  for(kk = 0 ; kk < n_vars ; kk++){
    
    hh = kk + jj * n_vars;
    out = out + beta[hh] * XX[ii][kk];
    
  }
  
  if(Log == 1){
    return(out);
  }
  
  return(exp(out));
  
}

void update_beta_jj(double **XX, double **JJ, double **loggamma, 
                    double *beta_temp, int *inclusion_indicator, 
                    double *prop_per_beta, double mu_be[n_cats][n_vars], 
                    double sig_be[n_cats][n_vars], double aa_hp, 
                    double bb_hp, int jj){
  
  // this function loops through n_vars so it is called for one taxa at a time
  
  double sig_prop;
  double loggamma_p[n_obs];
  int hh, ii, kk;
  
  // define proposal variable and acceptance ratio variables
  double beta_p;
  double lnu, ln_acp;
  
  double log_full_beta, log_full_beta_p;
  
  for(kk = 0 ; kk < n_vars ; kk++){ 
    
    // Stride for full (n_vars * n_cats) vector
    hh = kk + jj * n_vars;
    
    if(inclusion_indicator[hh] == 1){
      
      log_full_beta = 0;
      
      for(ii = 0 ; ii < n_obs ; ii++){
        
        log_full_beta = log_full_beta - gsl_sf_lngamma(exp(loggamma[ii][jj]));
        log_full_beta = log_full_beta + exp(loggamma[ii][jj]) * log(JJ[ii][jj]);
        
      }
      
      log_full_beta = log_full_beta + lprior_bbsas(beta_temp[kk], inclusion_indicator[hh], 
                                                   sig_be[jj][kk], mu_be[jj][kk], aa_hp, bb_hp);
      sig_prop = prop_per_beta[hh];
      beta_p = beta_temp[kk] + adap_prop(sig_prop);
      
    
      for(ii = 0 ; ii < n_obs ; ii++){
        loggamma_p[ii] = loggamma[ii][jj] - beta_temp[kk] * XX[ii][kk] + beta_p * XX[ii][kk];
      }
      
      // calculate proposal probability
      log_full_beta_p = 0;
      for(ii = 0 ; ii < n_obs ; ii++){
        log_full_beta_p = log_full_beta_p - gsl_sf_lngamma(exp(loggamma_p[ii]));
        log_full_beta_p = log_full_beta_p + exp(loggamma_p[ii]) * log(JJ[ii][jj]);
      }
      log_full_beta_p = log_full_beta_p + lprior_bbsas(beta_p, inclusion_indicator[hh], sig_be[jj][kk], mu_be[jj][kk], aa_hp, bb_hp);
      
      ln_acp = log_full_beta_p - log_full_beta;
      lnu = log(gsl_rng_uniform(rando));
      
      if(lnu < ln_acp){
        
        beta_temp[kk] = beta_p;
        
        for(ii = 0 ; ii < n_obs ; ii++){
          
          loggamma[ii][jj] = loggamma_p[ii];
          
        }
        
      } // Close MH accept
    } // Close inclusion_indicator[hh] == 1
  } // Close for kk
} // Close function

void between_models_jj(double **XX, double **JJ, double **loggamma, 
                       double *beta_temp, int *accepted_beta_flag, 
                       int *inclusion_indicator, double mu_be[n_cats][n_vars],
                       double sig_be[n_cats][n_vars], double aa_hp, 
                       double bb_hp, int jj){
  
  // Metropolis-Hastings proposals
  double sig_prop;
  double mu_prop;
  
  int hh, kk, ii;
  
  // proposed value vectors
  int inclusion_indicator_p;
  double loggamma_p[n_obs];
  double beta_p;
  
  // acceptance ratio variables
  double lnu, ln_acp;
  
  // must update beta_jk for kk = 1, ..., n_vars
  double log_full_beta, log_full_beta_p;
  
  for(kk = 0 ; kk < n_vars ; kk++){
    
    hh = kk + jj * n_vars;
    
    // calculate the current log full conditional
    log_full_beta = 0;
    for(ii = 0 ; ii < n_obs ; ii++){
      log_full_beta = log_full_beta - gsl_sf_lngamma(exp(loggamma[ii][jj]));
      log_full_beta = log_full_beta + exp(loggamma[ii][jj]) * log(JJ[ii][jj]);
    }
    log_full_beta = log_full_beta + lprior_bbsas(beta_temp[kk], inclusion_indicator[hh], sig_be[jj][kk], mu_be[jj][kk], aa_hp, bb_hp);
    
    // proposing a new value for beta[jj][kk] using the prior mean and standard deviation
    sig_prop = pow(sig_be[jj][kk], 0.5);
    mu_prop = mu_be[jj][kk];
    
    // swap and sample beta
    if(inclusion_indicator[hh] == 0){
      beta_p = mu_prop + gsl_ran_gaussian_ziggurat(rando, sig_prop);
      inclusion_indicator_p = 1;
    }
    else{
      beta_p = 0.0;
      inclusion_indicator_p = 0;
    }
    
    // update proposal gamma (linear predictor)
    for(ii = 0 ; ii < n_obs ; ii++){
      loggamma_p[ii] = loggamma[ii][jj] - beta_temp[kk] * XX[ii][kk] + beta_p * XX[ii][kk];
    }
    
    // calculate the proposed log full conditional
    log_full_beta_p = 0;
    for(ii = 0 ; ii < n_obs ; ii++){
      log_full_beta_p = log_full_beta_p - gsl_sf_lngamma(exp(loggamma_p[ii]));
      log_full_beta_p = log_full_beta_p + exp(loggamma_p[ii]) * log(JJ[ii][jj]);
    }
    log_full_beta_p = log_full_beta_p + lprior_bbsas(beta_p, inclusion_indicator_p, sig_be[jj][kk], mu_be[jj][kk], aa_hp, bb_hp);
    
    ln_acp = log_full_beta_p - log_full_beta;
    lnu = log(gsl_rng_uniform(rando));
    
    // Metropolis-Hastings ratio
    if(lnu < ln_acp){
  
      accepted_beta_flag[hh] = 1;
      
      // update accepted beta into beta_temp
      beta_temp[kk] = beta_p;
      inclusion_indicator[hh] = inclusion_indicator_p;
      
      for(ii = 0 ; ii < n_obs ; ii++){
        // also update the gamma (linear predictor)
        loggamma[ii][jj] = loggamma_p[ii];  
      }
    } // Close MH if 
  } // Close for kk
} // Close function

void update_alpha_jj(double **JJ, double **loggamma, double *alpha, 
                     double *prop_per_alpha, int *accepted_alpha_flag, 
                     double mu_al[n_cats], double sig_al[n_cats], int jj){
  
  int ii;
  
  double sig_prop = prop_per_alpha[jj];
  double loggamma_p[n_obs];
  
  double alpha_p;
  double lnu, ln_acp;
  
  // Prepare the current and proposed full conditional values
  double log_full_alpha, log_full_alpha_p;
  
  // Calculate the full conditional for the current value
  log_full_alpha=0;
  for(ii = 0 ; ii < n_obs ; ii++){
  
    log_full_alpha = log_full_alpha - gsl_sf_lngamma(exp(loggamma[ii][jj]));
    log_full_alpha = log_full_alpha + exp(loggamma[ii][jj]) * log(JJ[ii][jj]);
  
  }
  log_full_alpha = log_full_alpha - 1.0/(2.0 * sig_al[jj]) * pow(alpha[jj] - mu_al[jj], 2.0);
  
  // Propose a new value for alpha[jj] using a random walk proposal centered on 
  // the current value of alpha[jj]
  alpha_p = alpha[jj] + gsl_ran_gaussian_ziggurat(rando, sig_prop);
  
  // Gamma must be updated too
  for(ii = 0 ; ii < n_obs ; ii++){
    loggamma_p[ii] = loggamma[ii][jj] - alpha[jj] + alpha_p;
  }
  
  // Calculate the full conditional for the proposed value
  log_full_alpha_p = 0;
  for(ii = 0 ; ii < n_obs ; ii++){
    log_full_alpha_p = log_full_alpha_p - gsl_sf_lngamma(exp(loggamma_p[ii]));
    log_full_alpha_p = log_full_alpha_p + exp(loggamma_p[ii]) * log(JJ[ii][jj]);
  }
  
  log_full_alpha_p = log_full_alpha_p - 1.0/(2.0 * sig_al[jj]) * pow(alpha_p - mu_al[jj], 2.0);
  ln_acp = log_full_alpha_p - log_full_alpha;
  lnu = log(gsl_rng_uniform(rando));
  
  if(lnu < ln_acp){
  
    // If accepted, update both alpha[jj] and loggamma[ii][jj], and keep 
    // track of acceptances
    accepted_alpha_flag[jj] = 1;
    alpha[jj] = alpha_p;
  
    for(ii = 0 ; ii < n_obs ; ii++){
      loggamma[ii][jj] = loggamma_p[ii];
    }
  }
  
  return;
  
} // Close function

// the probability of the Bernoulli trial is integrated out
double lprior_bbsas(double betajk, int sjk, double sig_bejk, double mu_bejk, 
                    double aa_hp, double bb_hp){
  
  // calculate additional beta factor
  double aa = sjk + aa_hp;
  double bb = 1 - sjk + bb_hp; 
  double lbeta_factor = gsl_sf_lnbeta(aa, bb) - gsl_sf_lnbeta(aa_hp, bb_hp);
  
  // piece-wise function
  if(sjk == 0){
    return(log(1.0 - exp(lbeta_factor)));
  }
  else{
    return(lbeta_factor - 0.5 * log(2.0 * M_PI * sig_bejk) - 1.0/(2.0 * sig_bejk) * pow(betajk - mu_bejk, 2.0));
  }
  
}

// simple adaptive MH proposal
// this uses equation (3) from Roberts & Rosenthal (2009) but without the full
// covariance of the target which makes it similar to the component-wise update
// of Haario et al. (2005)
double adap_prop(double curr_var){
  
  // need dimension (???)
  double dd = (double)n_cats * (double)n_vars;
  
  // ensures bounded convergence  
  int safeguard = gsl_ran_bernoulli(rando, 0.05);
  double usually = pow(2.38, 2.0) * curr_var/dd;
  double unusually = pow(0.1, 2.0)/dd;
  
  double prop_var = (1 - safeguard) * gsl_ran_gaussian_ziggurat(rando, pow(usually, 0.5)) + safeguard * gsl_ran_gaussian_ziggurat(rando, pow(unusually, 0.5));
  
  return(prop_var);
}

// sort of following http://www.johndcook.com/blog/standard_deviation/
double online_mean(int iteration, double last_mean, double curr_obs){
  
  // don't fall off by one
  int n_vals = iteration++;
  
  // assume that iteration > 0 since there's a burnin proposal
  double curr_mean = last_mean + (curr_obs - last_mean)/(double)n_vals;
  
  return(curr_mean);
}

double online_var(int iteration, double last_mean, double last_var, double curr_mean, double curr_obs){

  // note that iteration == n - 1 and that (n - 1) * last_var == last_ss
  // assume that iteration > 0 since there's a burnin proposal
  double curr_ss = (double)iteration * last_var + (curr_obs - last_mean) * (curr_obs - curr_mean);
  
  return(curr_ss/(double)iteration);
}
 
// ///////////////////////////////////////////////////
// Roy Burstein and Aaron Osgood-Zimmerman
// August 2017
// Template file for space-time-Z GPR model.
// Used for fitting IHME Geospatial MBG models
// ///////////////////////////////////////////////////

// ///////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. In our nomenclature, Z is a third interaction (ie age) which defaults to AR1
// 3. Requires same space mesh for all time-Z points
// 4. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 5. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// ///////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;


// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// Robust Inverse Logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // cooresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}


// AR funtion from neal m
template<class Type>
SparseMatrix<Type> ar_Q(int N, Type rho, Type sigma) {
  SparseMatrix<Type> Q(N,N);
  Q.insert(0,0) = (1.) / pow(sigma, 2.);
  for (size_t n = 1; n < N; n++) {
    Q.insert(n,n) = (1. + pow(rho, 2.)) / pow(sigma, 2.);
    Q.insert(n-1,n) = (-1. * rho) / pow(sigma, 2.);
    Q.insert(n,n-1) = (-1. * rho) / pow(sigma, 2.);
  }
  Q.coeffRef(N-1,N-1) = (1.) / pow(sigma, 2.);
  return Q;
}

// objective function (ie the likelihood function for the model), returns the evaluated negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ////////////////////////////////////////////////////////////////////////////
  // INPUTS
  // ////////////////////////////////////////////////////////////////////////////
  DATA_INTEGER(flag); // flag=0 => only prior

  // Indices
  DATA_INTEGER(num_i);       // number of datapts in space-time-Z (aka STZ)
  DATA_INTEGER(num_s);       // number of mesh pts in space mesh
  DATA_INTEGER(num_t);       // number of time periods
  DATA_INTEGER(num_z);       // number of Z groups

  // Data (each, excpect for X_ij is a vector of length num_i)
  DATA_VECTOR(y_i);          // obs successes per binomial experiment at point i (aka cluster)
  DATA_VECTOR(n_i);          // trials per cluster
  DATA_IVECTOR(t_i);         // time period of the data point
  DATA_IVECTOR(w_i);         // weights for observations
  DATA_MATRIX(X_ij);         // covariate design matrix (num_i by number of fixed effects matrix)

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Options
  DATA_VECTOR(options)       // boolean vector of options to be used to select different models/modelling options:
                             // 0: If 1, Include priors. All are default settings right now
                             // 1: If 1, ADREPORT is on. Used for testing for now
                             // 2: If 1, use nugget

  // Parameters
  PARAMETER_VECTOR(alpha_j);   // fixed effect coefs, including intercept as first index
  PARAMETER(logtau);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho_trans);             // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho_trans);             // Z autocorrelation parameter for AR1, natural scale
  PARAMETER(log_nugget_sigma); // log of the standard deviation of the normal error nugget term


  // Random effects
  PARAMETER_ARRAY(Epsilon_stz); // Random effects for each STZ mesh location. Should be 3D array of dimensions num_s by num_t by num_z
  PARAMETER_VECTOR(nug_i);       // Random effects of the nugget
  
  printf("Epsilon_stz size: %d \n", Epsilon_stz.size());

  // ////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  // ////////////////////////////////////////////////////////////////////////////

  // Define the joint-negative log-likelihood as a parallel_accumulator	   // Define the joint-negative log-likelihood as a parallel_accumulator
  // this allows us to add or subtract numbers to the object in parallel	   // this allows us to add or subtract numbers to the object in parallel
  // parallel_accumulator<Type> jnll(this);	   // parallel_accumulator<Type> jnll(this);
  vector<Type> jnll_comp(3);
  jnll_comp[0] = Type(0); // priors contribution	   jnll_comp[0] = Type(0); // priors contribution
  jnll_comp[1] = Type(0); // latent field contrib	   jnll_comp[1] = Type(0); // latent field contrib
  jnll_comp[2] = Type(0); // data contrib

  // print parallel info
  // max_parallel_regions = omp_get_max_threads();
  // printf("This is thread %d\n", max_parallel_regions);
  max_parallel_regions = 5;

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(logkappa, logtau, M0, M1, M2);
  // printf("Q_ss size: %d \n", Q_ss.size());

  // Make transformations of some of our parameters
  Type range     = sqrt(8.0) / exp(logkappa);
  Type sigma     = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  Type trho      = (exp(trho_trans) - 1) / (exp(trho_trans) + 1); // TRANSOFRM from -inf, inf to -1, 1.. //log((1.1 + trho) / (1.1 - trho));
  Type zrho      = (exp(zrho_trans) - 1) / (exp(zrho_trans) + 1); //TRANSOFRM from -inf, inf to -1, 1.. // log((1.1 + zrho) / (1.1 - zrho));
  Type nugget_sigma  = exp(log_nugget_sigma);
  
  // Define objects for derived values
  vector<Type> fe_i(num_i);                         // main effect X_ij %*% t(alpha_j)
  vector<Type> epsilon_stz(num_s * num_t * num_z);  // Epsilon_stz unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);                // value of gmrf at data points
  vector<Type> prob_i(num_i);                       // Logit estimated prob for each point i

  // Prior contribution to likelihood. Values are defaulted (for now). Only run if options[0]==1
  if(options[0] == 1) {
   PARALLEL_REGION jnll_comp[0] -= dnorm(logtau,    Type(0.0), Type(1.0),   true);  // N(0,1) prior for logtau
   PARALLEL_REGION jnll_comp[0] -= dnorm(logkappa,  Type(0.0), Type(1.0),   true);  // N(0,1) prior for logkappa
   if(num_t > 1) {
     PARALLEL_REGION jnll_comp[0] -= dnorm(trho_trans, Type(0.0), Type(2.582), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   if(num_z > 1) {
     PARALLEL_REGION jnll_comp[0] -= dnorm(zrho_trans, Type(0.0), Type(2.582), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   for( int j = 0; j < alpha_j.size(); j++){
     PARALLEL_REGION jnll_comp[0] -= dnorm(alpha_j(j), Type(0.0), Type(3), true); // N(0, sqrt(1/.001)) prior for fixed effects.
   }
  }

  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: S, ST, SZ, and STZ
  if (num_t == 1 & num_z == 1)  {
    printf("GP FOR SPACE  ONLY \n");
    PARALLEL_REGION jnll_comp[1] += GMRF(Q_ss,false)(epsilon_stz);
  } else if(num_t > 1 & num_z == 1) {
    printf("GP FOR SPACE-TIME \n");
    PARALLEL_REGION jnll_comp[1] += SEPARABLE(AR1(trho),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t == 1 & num_z > 1) {
    printf("GP FOR SPACE-Z \n");
    PARALLEL_REGION jnll_comp[1] += SEPARABLE(AR1(zrho),GMRF(Q_ss,false))(Epsilon_stz);
  } else if (num_t > 1 & num_z > 1) {
    printf("GP FOR SPACE-TIME-Z \n");
    PARALLEL_REGION jnll_comp[1] += SEPARABLE(AR1(zrho),SEPARABLE(AR1(trho),GMRF(Q_ss,false)))(Epsilon_stz);
  }

  // nugget contribution to the likelihood
  if(options[2] == 1 ){
    printf("adding in Nugget \n");
    for (int i = 0; i < num_i; i++){
      PARALLEL_REGION jnll_comp[1] -= dnorm(nug_i(i), Type(0.0), nugget_sigma, true);
    }
  }

  // Transform GMRFs and make vector form
  for(int s = 0; s < num_s; s++){ // space

    if(num_t == 1) { // single time
      epsilon_stz[(s)] = Epsilon_stz(s);
    } else{
      for(int t = 0; t < num_t; t++){ // more than one time
	if(num_z == 1) { // single z 
	  epsilon_stz[(s + num_s * t )] = Epsilon_stz(s,t);
	} else {
	  for(int z = 0; z < num_z; z++){ // more than one z
	    epsilon_stz[(s + num_s * t + num_s * num_t * z)] = Epsilon_stz(s,t,z);
	  } // z-dim
	}
      } // time
    }
  } // space

  // Project from mesh points to data points in order to eval likelihood at each data point
  projepsilon_i = Aproj * epsilon_stz.matrix();

  // evaluate fixed effects for alpha_j values
  fe_i = X_ij * alpha_j.matrix();

  // Return un-normalized density on request
  if (flag == 1){
    printf("Returning before data likelihood b/c flag == 1\n");
    printf("jnll of priors:       %f\n", asDouble(jnll_comp(0)));
    printf("jnll of gmrf & nug:   %f\n", asDouble(jnll_comp(1)));
    printf("jnll of data:         %f\n", asDouble(jnll_comp(2)));
    Type jnll = jnll_comp.sum();
    printf("Combined jnll is: %f\n", asDouble(jnll));

    return jnll;
  } 

  // Likelihood contribution from each datapoint i
  printf("Data likelihood \n")
  for (int i = 0; i < num_i; i++){

    prob_i(i) = fe_i(i) + projepsilon_i(i) + nug_i(i);

    if(!isNA(y_i(i))){
      PARALLEL_REGION jnll_comp[2] -= dbinom( y_i(i), n_i(i), invlogit_robust(prob_i(i)), true ) * w_i(i);
    }
    
  }

  // to help with debug, print each loglik component
  printf("jnll of priors: %f\n", asDouble(jnll_comp(0)));
  printf("jnll of gmrf:   %f\n", asDouble(jnll_comp(1)));
  printf("jnll of data:   %f\n", asDouble(jnll_comp(2)));

  // sum logliks
  Type jnll = jnll_comp.sum();
  printf("Combined jnll is: %f\n", asDouble(jnll));

  
  // Report estimates
  if(options[1] == 1){
    ADREPORT(alpha_j);
    ADREPORT(Epsilon_stz);
  }

  return jnll;
}

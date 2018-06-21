// /////////////////////////////////////////////////////////////////////////////
// Roy Burstein, Aaron Osgood-Zimmerman, Nat Henry
// March 2018
// Template file for fitting combined Birth History - VR datasets for U5M
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. In our nomenclature, Z is a third interaction (ie age) which defaults to AR1
// 3. Requires same space mesh for all time-Z points
// 4. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 5. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// /////////////////////////////////////////////////////////////////////////////

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



template<class Type>
Type objective_function<Type>::operator() ()
{
  // Indices
  DATA_INTEGER( num_i );   // Number of data points in space
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh

  // Demographic information to convert between 5q0 and 5m0
  // DATA_SCALAR( nAx );          // nAx
  // DATA_SCALAR( bin_width );    // Number of years in the age bin (aka 5 yrs for 5q0)

  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_i );    // Num occurrences (deaths) per binomial experiment at point i (cluster)
  DATA_VECTOR( Exp_i );   // Trials per cluster
  // DATA_VECTOR( is_vr );   // Whether or not the data is VR
  DATA_MATRIX( X_ij );    // Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(Aproj);   // Used to project spatial mesh to data locations

  // Options
  DATA_VECTOR( options );

  // Fixed effects
  PARAMETER_VECTOR(alpha_j); // Fixed effect coefficients, including intercept as first index
  PARAMETER(log_tau);        // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER(log_kappa);      // Log of INLA kappa (related to spatial correlation and range)
  // PARAMETER(logit_pi_vr);    // Completeness of VR data, currently fixed across space

  // Random effects
  PARAMETER_ARRAY(Epsilon_s);  // Random effect for each spatial mesh location. Currently a 1d array of num_s

  // objective function -- joint negative log-likelihood
  Type jnll = 0;

  // print parallel info
  max_parallel_regions = omp_get_max_threads();
  printf("This is thread %d\n", max_parallel_regions);

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);

  // Transform some of our parameters
  Type range = sqrt(8.0) / exp(log_kappa);
  Type sigma = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * log_tau) * exp(2.0 * log_kappa));
  // Type pi_vr = invlogit(logit_pi_vr);

  // Define objects for derived values
  vector<Type> epsilon_s(num_s);         // Epsilon_s (array) unlisted into a vector 
                                         //   for easier matrix multiplication
  vector<Type> fe_i(num_i);              // main effect X_ij %*% t(alpha_j)
  vector<Type> projepsilon_i(num_i);     // value of gmrf at data points
  vector<Type> prob_i(num_i);            // Logit estimated prob for each point i
  // Vectors used for convenience when calculating jnll contributions of data points
  //vector<Type> q_val(num_i);
  //vector<Type> m_val(num_i);

  // Priors
  // Prior contribution to likelihood. Values are defaulted (for now). Only run if options[0]==1
  if(options[0] == 1) {
   PARALLEL_REGION jnll -= dnorm(log_tau,   Type(0.0), Type(1.0), true); // N(0,1) prior for logtau
   PARALLEL_REGION jnll -= dnorm(log_kappa, Type(0.0), Type(1.0), true); // N(0,1) prior for logkappa
   // PARALLEL_REGION jnll -= dbeta(pi_vr,     Type(2.0), Type(2.0), true); // Beta distribution centered at .5 for VR completeness
   // Priors for each fixed effect
   for( int j = 0; j < alpha_j.size(); j++){
     PARALLEL_REGION jnll -= dnorm(alpha_j(j), Type(0.0), Type(3), true); // N(0, sqrt(1/.001)) prior for fixed effects.
   }
  }

  // Transform GMRFs and make vector form
  for(int s = 0; s < num_s; s++){
    epsilon_s[s] = Epsilon_s(s); // This is pretty easy since there's only 1 time period and 1 age bin
  }

  // Probability of Gaussian-Markov random fields (GMRFs)
  PARALLEL_REGION jnll += GMRF(Q_ss, false)(epsilon_s);

  // Project from mesh points to data points in order to eval likelihood at each data point
  projepsilon_i = Aproj * epsilon_s.matrix();

  // evaluate fixed effects for alpha_j values
  fe_i = X_ij * alpha_j.matrix(); 

  // Likelihood contribution from each datapoint i
  for (int i = 0; i < num_i; i++){
    prob_i(i) = fe_i(i) + projepsilon_i(i);
    if(!isNA(y_i(i))){
      // Transformed values to use in probability calculations
      //q_val(i)  = invlogit(prob_i(i));
      // m_val(i)  = q_val(i) / (bin_width + nAx * q_val(i) - bin_width * q_val(i));
      // Likelihood contribution from birth history datas
      //if (is_vr(i) == 1){
        // Likelihood contribution from VR data
        //PARALLEL_REGION jnll -= dpois( y_i(i), (Exp_i(i) * pi_vr * m_val(i)), true);
      //} else {
        // Likelihood contribution from non-VR data
        // Uses the dbinom_robust function, which takes the logit probability
        PARALLEL_REGION jnll -= dbinom_robust( y_i(i), Exp_i(i), prob_i(i), true );
	//}
    }
  }

  // Report estimates
  if(options[1] == 0){
    ADREPORT(alpha_j);
    // ADREPORT(logit_pi_vr);
    ADREPORT(Epsilon_s);
  }

  return jnll;
}





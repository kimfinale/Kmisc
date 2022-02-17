#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector reulermultinom(double size, NumericVector rate, double dt) {
  int ncol = rate.size();
  NumericVector trans(ncol); // transition events
  double p = sum(rate); //total event rate
  double tmpp = p;
  double tmpsize = R::rbinom(size, (1 - exp(-tmpp*dt))); // total number of events
  for (int k = 0; k < (ncol-1); k++) {
    double tr = R::rbinom(tmpsize, rate(k)/p);
    trans(k) = tr;
    tmpsize -= tr;
    tmpp -= rate(k);
  }
  trans(ncol-1) = tmpsize;
  return(trans);
}

// [[Rcpp::export]]
List sepiar_tauleap(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S(nsteps);
  NumericVector E(nsteps);
  NumericVector P(nsteps);
  NumericVector I(nsteps);
  NumericVector A(nsteps); //asymptomatic
  NumericVector X(nsteps); //isolated
  NumericVector R(nsteps); // recovered
  NumericVector CI(nsteps); // cumulative incidence
  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  List init = params["init"];
  S(0) = init["S"];
  E(0) = init["E"];
  P(0) = init["P"];
  I(0) = init["I"];
  A(0) = init["A"]; //asymptomatic
  X(0) = init["X"]; //isolated
  R(0) = init["R"]; // recovered
  CI(0) = init["CI"]; // cumulative incidence

  double delta = params["delta"]; // 1 / incubation period
  double gamma = params["gamma"]; // 1 / recovery period
  double epsilon = params["epsilon"]; // 1 / latent period
  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double rate_P_I = 1 / (1/delta - 1/epsilon);
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A compared to I
  double bP = params["bP"]; // relative infectivity of P compared to I
  double xA = params["xA"]; // relative isolation rate of A compared to I
  double xP = params["xP"]; // relative isolation rate of P compared to I
  double R0 = params["R0"];

  double rate_x = params["rate_isol"];
  double beta = R0 / (bP/(rate_P_I + rate_x*xP) +
                      bA*fA/(gamma + rate_x*xA) +
                      (1-fA)/(gamma + rate_x));

  // Rprintf("the value of beta : %f \n", beta);
  NumericVector P_rates = {rate_P_I*(1-fA), rate_P_I*fA};
  if (xP > 0) { // additional potential outcome (P -> X) if xP > 0
    P_rates.push_back(rate_x*xP);
  }
  // Rprintf("the value of P_rates(2): %f \n", P_rates(2));
  NumericVector A_rates = {gamma};
  if (xA > 0) {
    A_rates.push_back(rate_x*xA);
  }
  NumericVector I_rates = {gamma, rate_x};

  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < nsteps - 1; istep++) {

    double iS = S[istep];
    double iE = E[istep];
    double iP = P[istep];
    double iI = I[istep];
    double iA = A[istep];
    double iX = X[istep];
    double iR = R[istep];
    double iCI = CI[istep];

    // State Equations
    double N = iS + iE + iP + iI + iA + iX + iR;
    double inf = bP*iP + bA*iA + iI;
    double foi = beta * inf / N;

    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection = R::rbinom(iS, 1 - exp(-foi * tau));
    double EtoP = R::rbinom(iE, 1 - exp(- epsilon * tau)); //

    NumericVector from_P = reulermultinom(iP, P_rates, tau);
    double PtoI = from_P(0); //
    double PtoA = from_P(1); //
    double PtoX = 0; //
    if(xP > 0){
      PtoX = from_P(2);
    }

    NumericVector from_I = reulermultinom(iI, I_rates, tau);
    double ItoR = from_I(0); //
    double ItoX = from_I(1); //

    double AtoR = R::rbinom(iA, 1 - exp(- gamma * tau)); //
    double AtoX = 0;
    if (xA > 0) {
      NumericVector from_A = reulermultinom(iA, A_rates, tau);
      AtoR = from_A(0); //
      AtoX = from_A(1); //
    }

    // Calculate the change in each state variable
    double dS = - new_infection;
    double dE = new_infection - EtoP;
    double dP = EtoP - PtoI - PtoA - PtoX;
    double dI = PtoI - ItoR - ItoX;
    double dA = PtoA - AtoR - AtoX;
    double dR = ItoR + AtoR;
    double dX = PtoX + ItoX + AtoX;

    // Update next timestep
    S[istep + 1] = iS + dS;
    E[istep + 1] = iE + dE;
    P[istep + 1] = iP + dP;
    I[istep + 1] = iI + dI;
    A[istep + 1] = iA + dA;
    R[istep + 1] = iR + dR;
    X[istep + 1] = iX + dX;
    CI[istep + 1] = iCI + new_infection;// cumulative incidence
    time[istep + 1] = (istep + 1) * tau;// time in fractional years

  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E") = E,
    Named("P") = P,
    Named("I") = I,
    Named("A") = A,
    Named("R") = R,
    Named("X") = X,
    Named("CI") = CI);

  return sim;
};

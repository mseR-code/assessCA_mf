//********************************************************************
// Module: assessCA.cpp
// Authors: S.P. Cox, A.R. Kronlund, K.R. Holt, S. D. N. Johnson
// Procedure: Statistical catch-age model
//
// References: 
//
// Revised from single fleet to multifleet by S. D. N. Johnson
// Date Last Revised: Sept 27, 2018
// 
// Quick multifleet statistical catch at age model. Uses
// Popes approximation to approximate F for each fleet, with 
// a given fleet timing. Assumes no uncertainty in growth 
// model (all fish lie perfectly along vonB growth curve). 
// Modified from ADMB single fleet model assessCA.tpl
// 
// Features:
//  1. Pope's approx of F_tg values
//  2. Robust normal approx to multinomial likelihood
// 
// ToDo:
//  1. Check that the current prop'n at age calcs are correct:
//      - Look at difference bn vuln numbers normalised,
//        and vuln biomass scaled by wt then normalised
//  2. Consider changing age obs likelihood, or add an
//      effective sample size parameter
//  3. Ask Sean if he remembers why age obs nll was 
//      scaled by -0.5... Check Multifan-CL paper,
//      appears that likelihood had not been negged
//  4. Check beta prior moment matching scaling
//  5. Update SimulFleets code.
//  6. Add IG priors for rec, mort and obs err vars
//  7. Add normal prior for lnq
// 
//*******************************************************************/

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// invLogit
template<class Type>
Type invLogit(Type x, Type scale, Type trans)
{
  return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
}

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*\/\/\/\/\ DATA SECTION /\/\/\/\*/
  // Data Structures
  DATA_ARRAY(I_tg);             // CPUE data by gear, time (-1 == missing)
  DATA_ARRAY(C_tg);             // Catch data by gear, time
  DATA_ARRAY(A_atg);            // Age observations by age, gear, time ( -1 == missing )

  // Model dimensions
  int nG = Igt.dim(1);          // No. of surveys
  int nT = Igt.dim(0);          // No of time steps
  int nA = A_atg.dim(0);        // No of age classes

  // Model switches
  DATA_IVECTOR(survType_g);     // Type of index (0 = vuln bio, 1 = vuln numbers)
  DATA_IVECTOR(indexType_g);    // Type of survey (0 = relative, 1 = absolute)
  DATA_IVECTOR(calcIndex_g);    // Type of fleet (0 = survey, 1 = commercial) 
  DATA_IVECTOR(selType_g);      // Type of selectivity (0 = asymptotic, 1 = domed (normal))
  DATA_VECTOR(fleetTiming);     // Fraction of year before catch/survey observation
  DATA_IVECTOR(initCode);       // initialise at 0 => unfished, 1=> fished
  DATA_SCALAR(posPenFactor);    // Positive-penalty multiplication factor


  /*\/\/\/\/\ PARAMETER SECTION /\/\/\/\*/
  // Leading biological parameters //
  PARAMETER(lnB0);                      // Unfished spawning biomass
  PARAMETER(logit_ySteepness);          // SR steepness on (0,1)
  PARAMETER(lnM);                       // Natural mortality rates
  PARAMETER(log_initN_mult);            // Non-eq initial numbers multiplier (nA-vector)

  // Observation models //
  // Selectivity
  PARAMETER_VECTOR(lnSelAlpha_g);       // log scale selectivity alpha par (ageSel50 or gamma shape par)
  PARAMETER_VECTOR(lnSelBeta_g);        // log scale selectivity beta par (ageSel95 or gamma scale par)
  PARAMETER_VECTOR(lntauAge_g);         // Age observation sampling error log-SD
  PARAMETER_VECTOR(effSampleSize_g);    // Effetive sample size for age obs likelihood

  // Survey obs error and catchabilities are
  // concentrated as cond. MLEs

  // Process errors //
  PARAMETER_VECTOR(omegaR_t);           // Recruitment log-deviations
  PARAMETER(lnsigmaR);                  // log scale recruitment SD          
  PARAMETER_VECTOR(omegaM_t);           // Natural mortality log-deviations
  PARAMETER(lnsigmaM);                  // log scale natural mortality SD

  // Priors //
  PARAMETER_VECTOR(agetau2IGa);         // Inverse Gamma Prior alpha for age sampling error var prior
  PARAMETER_VECTOR(agetau2IGb);         // Inverse Gamma Prior beta for age sampling error var prior
  PARAMETER_VECTOR(obstau2IGa);         // Inverse Gamma Prior alpha for stock index obs error var prior
  PARAMETER_VECTOR(obstau2IGb);         // Inverse Gamma Prior beta for stock index obs error var prior
  PARAMETER_VECTOR(sig2RPrior);         // Hyperparameters for sig2R prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sig2MPrior);         // Hyperparameters for sig2M prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(rSteepBetaPrior);    // pars for Beta dist steepness prior
  PARAMETER_VECTOR(initMPrior);         // pars for initial M prior (mean,sd)

  // Fixed LH parameters
  PARAMETER_VECTOR( aMat );             // Maturity ogive parameters (ages at 50% and 95% mature)
  PARAMETER( Linf );                    // Asymptotic length
  PARAMETER( L1 );                      // Length at age 1
  PARAMETER( vonK );                    // vonB growth rate
  PARAMETER_VECTOR( lenWt );            // Allometric length/weight c1, c2 parameters
  
  

  // Derived Variables //
  // Transformed scalar parameters
  Type B0         = exp(lnB0);
  Type M          = exp(lnM);
  Type ySteepness = invLogit(logit_ySteepness,Type(1), Type(0.));
  Type rSteepness = (ySteepness + 0.3)/1.3;
  Type sigmaR     = exp(lnsigmaR);
  Type sigmaM     = exp(lnsigmaM);

  // Transformed parameter vectors
  vector<Type> SelAlpha_g   = exp(lnSelAlpha_g);
  vector<Type> SelBeta_g    = exp(lnSelBeta_g);
  vector<Type> tauAge_g     = exp(lntauAge_g);
  vector<Type> initN_mult   = exp(log_initN_mult);

  // Stock recruit model
  Type rec_a;             // BH a par
  Type rec_b;             // BH b par
  Type phi;               // ssbpr
  Type R0;                // unfished eqbm recruitment


  // Optimisation and modeling quantities
  Type objFun = 0.;             // Objective function
  Type posPen = 0.;             // posFun penalty
  Type obsIdxNLL = 0.;          // observation indices NLL
  Type ageObsNLL = 0.;          // Age observation NLL
  Type recNLL = 0.;             // recruitment deviations NLL
  Type mortNLL = 0.;            // Mt devations NLL


  // Life schedules
  vector<Type>  age(nA);
  vector<Type>  wt(nA);
  vector<Type>  len(nA);
  vector<Type>  mat(nA);
  array<Type>   surv(nA);

  // Selectivity values
  array<Type>   sel_ag(nA,nG);   

  // State variables
  array<Type>   N_at(nA,nT+1);          // Numbers at age
  array<Type>   B_at(nA,nT+1);          // Total biomass at age
  array<Type>   vulnB_atg(nA,nT,nG);    // Vulnerable biomass at age by gear
  array<Type>   vulnB_tg(nT,nG);        // Vulnerable biomass by gear
  array<Type>   vulnN_atg(nA,nT,nG);    // Vulnerable numbers at age by gear
  vector<Type>  SB_t(nT+1);             // Spawning biomass
  vector<Type>  R_t(nT+1);              // Recruitments
  vector<Type>  M_t(nT);                // Natural mortality vector
  array<Type>   F_tg(nT,nG);            // Fleet fishing mortality
  array<Type>   uAge_atg(nA,nT,nG);     // Vuln proportion at age in each gear at each time
  array<Type>   catAge_atg(nA,nT,nG);   // Catch at age in each gear at each time step
  array<Type>   predPA_atg(nA,nT,nG);   // Predicted proportion at age in each gear at each time step

  // Management quantities
  Type          termDep;        // Terminal depletion
  Type          projSpawnBio;   // Projected spawning biomass
  vector<Type>  projExpBio;     // Projected exploitable biomass

  // Observation model quantities
  Type          zSum;           // sum of log(It/Bt)
  Type          validObs;       // Number of observations
  array<Type>   z_tg(nT,nG);    // array to hold residuals for NLL calc

  // Conditional MLEs for survey observations
  vector<Type> lnqhat_g(nG);    // log scale q_g
  vector<Type> qhat_g(nG);      // natural scale
  vector<Type> lntauObs_g(nG);  // log scale obs err SD
  vector<Type> tauObs_g(nG);    // obs err SD


  /*\/\/\/\/\ PROCEDURE SECTION /\/\/\/\*/
  // Calculate life schedules - growth, maturity
  // Fill ages vector
  for( int a = 0; a < nA; a++)
    age(a) = a+1;
  // Calculate length curve
  len = Linf + (L1-Linf)*exp(-vonK*(age-1.));
  // Weight
  wt  = c1*pow(len,c2);
  // Maturity
  tmp = log(19.)/( aMat(1) - aMat(0) );
  mat = 1./( 1. + exp(-tmp*( age - aMat(0) ) ) );    

  // Calculate Mortality time series
  // First and last year uses the
  // average/initial M
  Mt(0) = M;
  Mt(nT) = M;
  for( int t = 1; t < nT; t++ )
    Mt(t) = Mt(t-1)*exp(omegaM_t(t));


  // Calculate selectivity
  for( int g = 0; g < nG; g++ )
  {
    // Check seltype switch (asymptotic vs dome)
    if( selType_g(g) == 0)
    {
      // asymptotic
      Type aSel50   = SelAlpha_g(g);
      Type aSel95   = SelAlpha_g(g) + SelBeta_g(g);
      Type tmp      = log(19.)/( aSel95 - aSel50 );
      sel_ag.col(g) = 1./( 1. + exp(-tmp*( age - aSel50 ) ) );
    }

    if( selType_g(g) == 1)
    {
      // domed (un-normalised gamma - or maybe normal ??)
      Type dSelAlpha  = exp ( SelAlpha_g(g));
      Type dSelBeta   = exp ( SelBeta_g(g) );
      sel_ag.col(g)   = pow(age * dSelBeta / (dSelAlpha - 1),dSelAlpha-1)*
                        exp ( dSelAlpha - 1 - age * dSelBeta );
    }    
  }

  // Loop through all fleets  in 2 nested loops, and build a 
  // new vector of indices in chronological order
  Type          minTime = 0;
  Type          prevTime = 0;
  vector<Type>  usedFleet(nG);
  vector<int>   chronIdx(nG);

  Type          lowG = 0;

  for( int pIdx = 0; p < nG; p++ )
  {
    // Loop over fleet timing vector, pull
    // get gear idx of minimum year fraction
    for( gIdx = 0; g < nG; g++ )
    {
      if( ( fleetTiming(gIdx) >= minTime) & 
          ( fleetTiming(gIdx) <= fleetTiming(lowG) ) &
          ( usedFleet(gIdx)  == 0 ) )
      {
        // Record index of new lowest time that hasn't been used
        lowG      = gIdx;
      }
    }
    chronIdx(pIdx)  = lowG;
    usedFleet(lowG) = int(1);
    prevTime        = minTime;
    minTime         = fleetTiming(lowG);
  }


  // Calculate SSBpr, recruitment parameters
  // Calculate equilibrium unfished survivorship.
  surv(0) = 1;
  for(int a = 1; a < nA; a++ )
    surv(a) = exp(-M*a);

  surv(nA-1) /= (1.0 - exp(-M));

  //phi calc
  phi = sum( mat * wt * surv );

  // Now compute R0 from B0
  R0 = B0/phi;

  // Beverton-Holt a parameter.
  rec_a = 4.*rSteepness*R0/(B0*(1.-rSteepness));
  // Beverton-Holt b parameter.
  rec_b = (5.*rSteepness-1.)/(B0*(1.-rSteepness));
    
  // Initialise population
  N_at.col(0) = R0 * surv;
  if( initCode == 1 )
    Nat.col(0) *= initN_mult; 

  // Calc biomass and spawning biomass in first year
  Bat.col(0) = Nat.col(0) * wt;
  SBt(0) = sum(Bat.col(0) * mat);


  // Loop over time steps, run pop dynamics
  for( int t = 1; t <= nT; t++ )
  {
    // First, get the previous year's beginning of year numbers
    vector<Type>  tmpN_at = N_at.col(t-1);

    // Now compute loop over fleets and take
    // catch as necessary
    Type prevTime = 0.;
    for( cIdx = 0; cIdx < nG; cIdx ++ )
    {
      // Get actual fleet number
      gIdx = chronIdx(cIdx);
      // Get fraction of M being used to reduce Nat
      Type fracM = fleetTiming(gIdx) - prevTime;

      // First, vulnerable numbers is found by reducing
      // numbers by fracM
      tmpN_at = tmpN_at * exp( - fracM * Mt(t-1) );
      vulnN_atg.col(gIdx).col(t-1) = tmpN_at * sel_ag
      vulnB_atg.col(gIdx).col(t-1) = vulnN_atg.col(gIdx).col(t-1)*wt;
      vulnB_tg.col(gIdx).col(t-1) = vulnB_atg.col(gIdx).col(t-1).sum();

      // Caclulate proportion of each age in each fleet's 
      // vuln biomass
      uAge_atg.col(gIdx).col(t-1) = vunB_atg.col(gIdx).col(t-1) / vulnB_tg(t-1,gIdx);

      // Get numbers caught at age
      catAge_atg.col(gIdx).col(t-1) = uAge_atg.col(gIdx).col(t-1) * C_tg(t-1,gIdx) / wt;

      // Calculate F - there might be a better way here...
      // Read MacCall's paper on alternatives to Pope's approx
      // Question is: what do we use for estimating F? Just U? Crank this
      // on paper first.
      // Pope's approx works better within a single age class, or
      // with DD models (F = log(Nt/Nt-1)-M)
      F_tg(t-1,gIdx) = - log( 1 - C_tg(t-1, gIdx) / vulnB_tg(t-1,gIdx) );

      // Remove numbers from the total population
      // SimulFleets //
      // Here is where we'll need to take care of
      // the vulnerable biomass for simultaneous fleet timings
      // 1. Check if the next fleet's timing is the same as this
      // fleet's
      // 2. If so, use the previous numbers to
      // calculate vuln biomass at the same timing
      for( int a = 0; a < nA; a++)
      {
        tmpN_at(a) -= catAge_atg(a,t-1,gIdx);
        tmpN_at(a) =  posfun(tmpN_at(a), 1e-3, posPen );
      }
    }

    // Finally, advance numbers at age
    // to the following time step by reducing
    // remaining numbers by the remaining mortality,
    // and adding recruitment
    Type lastFrac = 1 - prevTime;
    for(int a = 0; a < nA; a++ )
    {
      // Recruitment
      if( a == 0 )
      {
        N_at(a,t) = rec_a*SB_t(t-1)/( 1. + rec_b*SB_t(t-1) );
        if( t < nT )
          N_at(a,t) *= exp(omegaR_t(t));
      } else {
        // Depletion by lastFrac*Mt
        N_at(a,t) = tmpN_at(a-1) * exp( - lastFrac * Mt(t-1))  
        if( a == nA - 1)
          N_at(nA-1,t) += tmpN_at(nA-1) * exp( - lastFrac * Mt(t-1))
      }
      // Convert to biomass
      B_at(a,t) = N_at(a,t) * wt;
    }

    // Caclulate spawning biomass for next time step
    SB_t(t) = sum(B_at.col(t) * mat);

  } // End of population dynamics
  
  // Calculate likelihood functions
  // Observation models //
  
  // Initialise lnqhat and lntauObs at 0
  lnqhat_g.fill(0.0);
  lntauObs_g.fill(0.0);
  z_tg.fill(0.0);
  int t = 0; int gIdx = 0;
  // Loop over gear types
  for( gIdx = 0; gIdx < nG; g++ )
  {
    validObs = 0.;
    zSum = 0.;
    // Stock indices (concentrate obs error var and lnq)
    // Check if this is a survey
    if( calcIndex_g(gIdx) == 1)
    {
      // Loop over time steps
      for( t = 0; t < nT; t++ )
      {
        if( I_tg(t,gIdx) > 0.0 )
        {
          // Recover numbers or biomass
          if( survType_g(gIdx) == 1 )
            idxState = vulnN_atg.col(gIdx).col(t).sum();
          if( survType_g(gIdx) == 0 )
            idxState = vulnB_tg(t,gIdx);

          // Calculate residual
          z_tg(t,gIdx) = log(I_tg(t,gIdx)) - log(idxState);
          // Add to sum of residuals
          zSum += z_tg(t,gIdx);
          validObs += 1;
        }
      }
      // Calculate conditional MLE of q
      // if a relative index
      if( indexType_g(gIdx) == 0 )
      {
        // mean residual is lnq (intercept)
        lnqhat_g(gIdx) = zSum / validObs;

        // Subtract mean from residuals for
        // inclusion in likelihood
        for(int t = 0; t < nT; t++)
          if( I_tg(t,gIdx) > 0.0)
            z_tg(t,gIdx) -= lnqhat_g(gIdx);
      }
      // Calculate conditional MLE of observation
      // index variance
      Type SSR = pow(z_tg.col(gIdx),2).sum();
      tauObs_g(gIdx) = sqrt(SSR/validObs);

      // Add concentrated nll value using cond MLE of tauObs
      obsIdxNLL += 0.5*(validObs * log( SSR / validObs ) + validObs);
    }

    // Age observations //

    // assessCA.tpl uses vulnerable
    // biomass for the predicted proportion at age
    // in the fleet and survey. My question
    // to myself at this point is: would 
    // some measure of catch at age be better? 
    // Given it's calculated based on the katch * uAge / wt,
    // which is no longer a proportion, and would 
    // change when renormalised by the sum, I think it's
    // better. It's also consistent with the modeling
    // choice of using a version of Pope's Approx

    vector<Type> predPropAge(nA);
    vector<Type> obsPropAge(nA);

    // Loop over time steps
    for( t = 0; t < nT; t++ )
    { 
      // Now estimate predicted catch-at-age
      // This was done already if catch > 0, but 
      // not for all fleets (fishery indep. surveys
      // would be missed)
      // First, calculate prop at age in each fleet
      predPropAge = vulnB_atg.col(gIdx).col(t) / vulnB_tg(t,gIdx);
      // Convert to numbers
      predPropAge /= wt;
      // Now renormalise to 
      predPropAge /= predPropAge.sum();
      // Save to array
      predPA_atg.col(gIdx).col(t) = predPropAge;

      // Check that age observations
      // exist by checking that the plus
      // group has observations
      if( A_atg(nA,t,gIdx) > 0 )
        for( int a=1; a<=nAges;a++ )
          if(A_atg(a,t,gIdx) >= 0)
          {
            // Robust normal approx. to multinomial likelihood from MULTIFAN-CL
            ageObsNLL += log(2.*3.14*(predPropAge(a)*(1.-predPropAge(a))+0.1/nA));
            ageObsNLL += log( exp( -effSampSize_g(gIdx)*pow( A_atg(a,t,gIdx)-predPropAge(a), 2. )/(2.*(1.-predPropAge(a))*predPropAge(a) +0.1/nA) )  +.01 );
          }
    }
    ageObsNLL *= (-0.5);
  }

  // Process error priors
  // Recruitment priors
  // Add recruitment deviations to rec NLL; sigmR is estimated
  recNLL += 0.5*((nT - 1) * lnsigmaR + pow(omegaR_t,2).sum()/sigmaR/sigmaR);

  // Now do a beta prior on steepness
  // Steepness prior pars passed in as a 2ple
  // of mean and sd, use moment matching to convert
  // to beta parameters. There's some scaling here
  // I don't fully understand, requires some thought
  Type muB            = 1.3*rSteepBetaPrior(0)-0.3;
  Type tauB           = muB*(1.-muB)/( 1.5625*pow(rSteepBetaPrior(1),2) )-1.; 
  Type aB             = tauB*muB; 
  Type bB             = tauB*(1.-muB); 
  Type steepness_nlp  = (-1.)*( (aB-1.)*log(ySteepness) + (bB-1.)*log(1.-ySteepness) );    
  recNLL              += steepness_nlp;

  // Add a recruitment var IG prior here?

  // Natural mortality prior
  mortNLL += 0.5*pow(M-initMprior(0),2)/pow(initMPrior(1),2);
  // Mt deviations; sigmaM is estimated
  mortNLL += 0.5*((nT - 1) * lnsigmaM + pow(omegaM_t,2).sum() / sigmaM / sigmaM);

  // Add mortality var IG prior here?
  
  // Add positive function penalty to objFun
  objFun += posPenFactor * posPen;

  // Add NLL contributions to objFun
  objFun += mortNLL + recNLL + obsIdxNLL + ageObsNLL;

  // Convert some values to log scale
  // for sd reporting
  vector<Type> lnSB_t       = log(SB_t);
  vector<Type> lnD_t        = log(SB_t) - lnB0;
  vector<Type> lnR_t        = log(R_t);
  vector<Type> lnM_t        = log(M_t);
  array<Type>  lnF_tg       = log(F_tg);


  /*\/\/\/\/\ REPORTING SECTION /\/\/\/\*/
  // Variables we want SEs for
  ADREPORT(lnSB_t);
  ADREPORT(lnR_t);
  ADREPORT(lnM_t);
  ADREPORT(lnM);
  ADREPORT(lnqhat_g);
  ADREPORT(lnD_t);
  ADREPORT(logit_ySteepness);
  ADREPORT(log_initN_mult);

  
  // Everything else //
  // State variables
  REPORT(B_at);
  REPORT(N_at);
  REPORT(vulnB_tg);
  REPORT(SB_t);
  REPORT(M_t);
  REPORT(F_tg);
  
  // Data
  REPORT(C_tg);
  REPORT(I_tg);
  REPORT(A_atg);

  // Random effects
  REPORT(omegaM_t);
  REPORT(omegaR_t);
  REPORT(initN_mult);
  
  // Observation model quantities
  REPORT(qhat_os);
  REPORT(tauObs_g);
  REPORT(tauAge_g);
  REPORT(predPA_atg);

  // Likelihood quantities
  REPORT(objFun);
  REPORT(ageObsNLL);
  REPORT(obsIdxNLL);
  REPORT(recNLL);
  REPORT(mortNLL);
  REPORT(pospen);
  
  return objFun;
}
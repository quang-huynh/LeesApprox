



//template <class Type>
//Type objective_function<Type>::operator() () {

  using namespace SCA;
  
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_VECTOR(M);         // Natural mortality at age
  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  
  DATA_INTEGER(use_LeesEffect);
  
  // Only use if use_LeesEffect = TRUE
  DATA_INTEGER(Nbins);
  DATA_VECTOR(LenBins);   // Length bins
  
  DATA_INTEGER(ngtg);     // Number of growth type groups
  DATA_VECTOR(distGTG);
  DATA_VECTOR(rdist);
  
  DATA_MATRIX(WAA);       // Weight-at-age-and-GTG at the beginning of the year
  DATA_MATRIX(LAA);       // Length-at-age-and-GTG at the beginning of the year
  DATA_MATRIX(xout);      // Length bins and Length-at-age sorted by row
  
  DATA_IMATRIX(interp_check);
  DATA_IMATRIX(interp_check2);
  DATA_IMATRIX(integ_check);
  DATA_IVECTOR(integ_fac);
  DATA_IVECTOR(integ_ind);
  
  DATA_MATRIX(SAA);
  DATA_VECTOR(Select_at_length);
  
  // Only use if use_LeesEffect = FALSE
  DATA_VECTOR(mean_WAA);
  DATA_VECTOR(Select_at_age2);

  PARAMETER(log_F);
  PARAMETER(Arec);
  PARAMETER(Brec);

  Type F = exp(log_F);

  vector<Type> NPR(max_age);
  vector<Type> Select_at_age(max_age);
  vector<Type> Weight_at_age(max_age);
  Select_at_age.setZero();
  Weight_at_age.setZero();

  if(use_LeesEffect) { // Get NPR, Weight_at_age, and Select_at_age
    // Split integration indices
    vector<vector<int> > integ_index = split(integ_ind, integ_fac);
    
    // Update Weight_at_age and NPR
    matrix<Type> NPR_GTG = LeesApp_fn(F, rdist, M, SAA, WAA, Weight_at_age, NPR, max_age, ngtg);
    matrix<Type> probLA = calcprob(LAA, NPR_GTG, xout, LenBins, max_age, ngtg, Nbins, interp_check, interp_check2,
                                   integ_check, integ_index);
    
    for(int a=0; a<max_age; a++) { // Get Select_at_age
      for(int len=0; len<Nbins; len++) Select_at_age(a) += probLA(a, len) * Select_at_length(len);
    }
  } else {
    Select_at_age = Select_at_age2;
    Weight_at_age = mean_WAA;
    NPR = calc_NPR(F, Select_at_age, M, max_age);
  }

  // Equilibrium reference points and per-recruit quantities
  vector<Type> CPR(max_age);
  Type EPR = 0;
  Type B = 0;
  Type VB = 0;
  Type Yield = 0;
  for(int a=0; a<max_age; a++) {
    Type Z_a = F * Select_at_age(a) + M(a);
    Type meanN = NPR(a) * (1 - exp(-Z_a))/Z_a;
    CPR(a) = Select_at_age(a) * F * meanN;

    EPR += NPR(a) * Weight_at_age(a) * mat(a);
    B += NPR(a) * Weight_at_age(a);
    VB += NPR(a) * Weight_at_age(a) * Select_at_age(a);
    Yield += CPR(a) * Weight_at_age(a);
  }

  Type R;
  if(SR_type == "BH") {
    R = Arec * EPR - 1;
  } else {
    R = log(Arec * EPR);
  }
  R /= Brec * EPR;

  B *= R;
  Yield *= R;
  Type N = R * NPR.sum();
  Type E = R * EPR;

  //REPORT(probLA);
  REPORT(F);
  REPORT(R);
  REPORT(B);
  REPORT(E);
  REPORT(N);
  REPORT(VB);
  REPORT(EPR);
  REPORT(Yield);
  REPORT(NPR);
  REPORT(CPR);
  REPORT(Select_at_age);

  return -1 * Yield;

//}

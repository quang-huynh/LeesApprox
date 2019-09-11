



//template <class Type>
//Type objective_function<Type>::operator() () {

  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_VECTOR(M);         // Natural mortality at age

  DATA_MATRIX(WAA);       // Weight-at-age-and-GTG at the beginning of the year
  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used

  DATA_INTEGER(ngtg);     // Number of growth type groups
  DATA_INTEGER(Nbins);
  DATA_VECTOR(LenBins);   // Length bins

  DATA_MATRIX(SAA);
  DATA_VECTOR(Select_at_length);

  DATA_MATRIX(LAA);       // Length-at-age-and-GTG at the beginning of the year
  DATA_MATRIX(xout);      // Length bins and Length-at-age sorted by row

  DATA_VECTOR(distGTG);
  DATA_VECTOR(rdist);

  DATA_IMATRIX(interp_check);
  DATA_IMATRIX(interp_check2);
  DATA_IMATRIX(integ_check);
  DATA_IVECTOR(integ_fac);
  DATA_IVECTOR(integ_ind);

  DATA_INTEGER(use_LeesEffect);

  PARAMETER(log_F);
  PARAMETER(Arec);
  PARAMETER(Brec);

  Type F = exp(log_F);

  // Split integration indices
  vector<vector<int> > integ_index = split(integ_ind, integ_fac);

  ////// GTG calcs
  matrix<Type> NPR(max_age, ngtg);
  matrix<Type> probGTGA(max_age, ngtg);
  vector<Type> NPR_age(max_age);
  vector<Type> CPR_age(max_age);

  NPR = LeesApp_fn(F, rdist, M, SAA, max_age, ngtg);

  // Calculate prob L|A
  matrix<Type> probLA(max_age, Nbins);
  vector<Type> Select_at_age(max_age);
  Select_at_age.setZero();

  if(use_LeesEffect) {
    probLA = calcprob(LAA, NPR, xout, LenBins, max_age, ngtg, Nbins, interp_check, interp_check2,
                      integ_check, integ_index);
  } else {
    matrix<Type> NPR_F0(max_age, ngtg);
    NPR_F0 = LeesApp_fn(Type(0), rdist, M, SAA, max_age, ngtg);
    probLA = calcprob(LAA, NPR_F0, xout, LenBins, max_age, ngtg, Nbins, interp_check, interp_check2,
                      integ_check, integ_index);
  }

  // Equilibrium reference points and per-recruit quantities
  Type EPR = 0;
  Type B = 0;
  Type VB = 0;
  Type Yield = 0;
  for(int a=0; a<max_age; a++) {
    NPR_age(a) = NPR.row(a).sum();
    probGTGA.row(a) = NPR.row(a)/NPR_age(a);

    for(int len=0; len<Nbins; len++) Select_at_age(a) += probLA(a, len) * Select_at_length(len);
    Type meanN = NPR_age(a) * (1 - exp(-Select_at_age(a) * F - M(a))) / (Select_at_age(a) * F + M(a));
    CPR_age(a) = Select_at_age(a) * F * meanN;

    for(int g=0;g<ngtg;g++) {
      EPR += NPR(a,g) * WAA(a,g) * mat(a);
      B += NPR(a,g) * WAA(a,g);
      VB += NPR(a,g) * WAA(a,g) * Select_at_age(a);
      Yield += CPR_age(a) * probGTGA(a,g) * WAA(a,g);
    }
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

  REPORT(probLA);
  REPORT(F);
  REPORT(R);
  REPORT(B);
  REPORT(E);
  REPORT(N);
  REPORT(VB);
  REPORT(EPR);
  REPORT(Yield);
  REPORT(NPR);
  REPORT(CPR_age);
  REPORT(Select_at_age);

  return -1 * Yield;

//}

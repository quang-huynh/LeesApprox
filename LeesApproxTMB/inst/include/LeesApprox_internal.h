
//template <class Type>
//Type objective_function<Type>::operator() () {

  DATA_VECTOR(FVec);
  DATA_INTEGER(ngtg);
  DATA_VECTOR(LenBins);
  DATA_VECTOR(LenMids);
  DATA_VECTOR(ages);
  DATA_SCALAR(M);
  DATA_SCALAR(Linf);
  DATA_MATRIX(LAA);
  DATA_MATRIX(xout); // Length bins and Length-at-age sorted by row
  DATA_SCALAR(LFS);
  DATA_SCALAR(L5);
  DATA_SCALAR(Vmaxlen);
  DATA_INTEGER(maxage);

  DATA_VECTOR(distGTG);
  DATA_VECTOR(rdist);

  PARAMETER(p);

  // calculate selectivity-at-age for each GTG
  Type sl = (LFS - L5) / pow(-log2(0.05), 0.5);
  Type sr = (Linf - LFS) / pow(-log2(Vmaxlen), 0.5);

  matrix<Type> SAA(maxage, ngtg); // selectivity-at-age by GTG
  SAA = s_dnormal(LAA, LFS, sl, sr);

  // calculate selectivity-at-length
  vector<Type> Select_at_length(LenMids.size());
  Select_at_length = s_dnormal(LenMids, LFS, sl, sr);

  // loop over ages and calculate N per recruit for each age class
  matrix<Type> Ns(maxage, ngtg);
  Ns.row(0) = rdist;

  for (int age=1; age<maxage; age++) {
    int yr_st = maxage-age-1;
    vector<Type> Zs(ngtg);
    Zs.setZero();
    for(int g=0; g<ngtg; g++) {
      for (int age2=0; age2<age; age2++) Zs(g) += M + FVec(yr_st + age2) * SAA(age2,g);
      Ns(age,g) = Ns(0,g) * exp(-Zs(g));
    }
  }

  // Calculate prob L|A
  int Nbins = LenBins.size() - 1;
  matrix<Type> probLA(maxage, Nbins);
  probLA = calcprob_wrapper(LAA, Ns, xout, LenBins, maxage, ngtg, Nbins);

  // Calculate mean selectivity-at-age and mean length-at-age and fished length composition
  vector<Type> Select_at_age(maxage);
  vector<Type> Len_at_age(maxage);
  vector<Type> NAA(maxage);
  vector<Type> LenComp(Nbins);

  Select_at_age.setZero();
  LenComp.setZero();
  for(int age=0; age<maxage; age++) {
    Len_at_age(age) = (LAA.row(age) * Ns.row(age)).sum();

    NAA(age) = Ns.row(age).sum();
    Len_at_age(age) /= NAA(age);

    for(int len=0; len<Nbins; len++) {
      Select_at_age(age) += probLA(age, len) * Select_at_length(len);
      LenComp(len) += probLA(age, len) * Select_at_length(len) * NAA(age);
    }
  }

  LenComp /= NAA.sum();
  LenComp /= LenComp.sum();

  REPORT(probLA);
  REPORT(LenComp);
  REPORT(Ns);
  REPORT(SAA);
  REPORT(Len_at_age);
  REPORT(Select_at_age);
  REPORT(Select_at_length);
  REPORT(LAA);
  REPORT(LenMids);
  REPORT(LenBins);
  REPORT(NAA);

  return p;

//}


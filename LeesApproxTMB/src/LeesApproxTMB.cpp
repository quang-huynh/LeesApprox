#include <TMB.hpp>
#include "../inst/include/LeesApprox_TMB_fn.h"

template <class Type>
Type objective_function<Type>::operator() () {

  DATA_VECTOR(FVec);
  DATA_INTEGER(ngtg);
  //DATA_SCALAR(maxsd);
  DATA_VECTOR(LenBins);
  DATA_VECTOR(LenMids);
  DATA_VECTOR(ages);
  //DATA_SCALAR(binwidth);
  DATA_SCALAR(M);
  DATA_SCALAR(Linf);
  //DATA_SCALAR(K);
  //DATA_SCALAR(t0);
  //DATA_SCALAR(LinfCV);
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

  // create vector of Fs for last maxage years
  // F assumed 0 for years where no F are provided
  //NumericVector FVec2(maxage);
  //int nFs = FVec.length();
  //int ind = nFs-1;
  //for (int i=maxage-1; i>=0; i--) {
  //  if (ind >= 0) {
  //    FVec2(i) = FVec(ind);
  //  }
  //  ind -= 1;
  //}
  //vector<Type> FVec2(maxage);
  //FVec2.setZero();
  //int nFs = FVec.size();
  //int ind = nFs - 1;
  //for(int i=maxage-1; i>=0; i--) {
  //  FVec2(i) == CppAD::CondExpGe(Type(ind), Type(0), FVec(ind), Type(0));
  //  ind -= 1;
  //}

  // loop over ages and calculate N per recruit for each age class
  matrix<Type> Ns(maxage, ngtg);
  for(int g=0; g<ngtg; g++) Ns(0,g) = rdist(g);

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
  Len_at_age.setZero();
  NAA.setZero();
  LenComp.setZero();
  for(int age=0; age<maxage; age++) {
    for(int g=0; g<ngtg; g++) {
      Len_at_age(age) += LAA(age,g) * Ns(age,g);
      NAA(age) += Ns(age,g);
    }
    Len_at_age(age) /= NAA(age);

    for(int len=0; len<Nbins; len++) {
      Select_at_age(age) += probLA(age, len) * Select_at_length(len);
      LenComp(len) += probLA(age, len) * Select_at_length(len) * NAA(age);
    }
  }

  for(int len=0; len<Nbins; len++) LenComp(len) /= NAA.sum();
  Type LenCompSum = LenComp.sum();
  for(int len=0; len<Nbins; len++) LenComp(len) /= LenCompSum;

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
}


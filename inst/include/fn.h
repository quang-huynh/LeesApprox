
template <class Type>
Type dmultinom_robust(vector<Type> x, vector<Type> p, int give_log=0)
{
  vector<Type> xp1 = x+Type(1);
  Type logres = lgamma(x.sum() + Type(1)) - lgamma(xp1).sum();
  for(int i=0;i<p.size();i++) logres += CppAD::CondExpGt(p(i), Type(0), x(i)*log(p(i)), Type(0));
  if(give_log) return logres;
  else return exp(logres);
}


template <class Type>
Type log2(Type x) {
  return log(x)/log(2);
}


template <class Type>
vector<vector<int> > split(vector<int> x, vector<int> fac) {
  if (x.size() != fac.size()) Rf_error("x and fac must have equal length.");
  int nlevels = 0;
  for (int i = 0; i < fac.size(); i++)
    if (fac[i] >= nlevels) nlevels = fac[i] + 1;
  vector<vector<int> > ans(nlevels);
  vector<int> lngt(nlevels);
  lngt.setZero();
  for (int i = 0; i < fac.size(); i++) lngt[fac[i]]++;
  for (int i = 0; i < nlevels; i++) ans[i].resize(lngt[i]);
  lngt.setZero();
  for (int i = 0; i < fac.size(); i++) {
    ans[fac[i]][lngt[fac[i]]] = x[i];
    lngt[fac[i]]++;
  }
  return ans;
}

// Interpolated abundance at each xout (ordered concatenation of LenBin and LAA for each GTG)
// x = LAA, y = N
template <class Type>
vector<Type> linear_int(matrix<Type> x, matrix<Type> y, matrix<Type> xout, int a,
                        matrix<int> interp_check, matrix<int> interp_check2) {
  int nout = xout.cols();
  vector<Type> yout(nout);
  for(int i=0;i<nout;i++) {
    if(interp_check(a, i) == -2) { // v is within x
      Type x1 = x(a, interp_check2(a, i));
      Type x2 = x(a, interp_check2(a, i)+1);
      Type y1 = y(a, interp_check2(a, i));
      Type y2 = y(a, interp_check2(a, i)+1);
      Type m = y2 - y1;
      m /= x2 - x1;
      yout(i) = m * (xout(a, i) - x1) + y1;
    } else if(interp_check(a, i) >= 0) { // v == one of x
      yout(i) = y(a, interp_check(a, i));
    } else yout(i) = 1e-8; // v outside range of x
  }
  return yout;
}

//calcprob_internal
template <class Type>
Type calcprob_internal2(matrix<Type> x, matrix<Type> y, int i, int j, int a) {
  Type h1 = abs(y(a, j) - y(a, i));
  Type by = x(a, j) - x(a, i);
  Type h2 = CppAD::CondExpLt(y(a, j), y(a, i), y(a, j), y(a, i));
  return h2 * by + 0.5 * by * h1;
}

template <class Type>
Type calcprob_internal(matrix<Type> x, matrix<Type> y, vector<vector<int> > integ_index, int a, int i,
                       int nclasses) {
  Type area = 0;
  int ii = a * nclasses + i;
  for(int j=0;j<integ_index(ii).size()-1;j++) {
    area += calcprob_internal2(x, y, integ_index(ii)(j), integ_index(ii)(j+1), a);
  }
  return area;
}


// Determine if the j-th xout entry (sorted concatenation of length bin and LAA for each GTG)
// is between the i-th and i-th + 1 length bin
// If yes for only one j, return the interpolated value. Else, sum across j.
// x = LAA, y = N
template <class Type>
matrix<Type> calcprob(matrix<Type> LAA, matrix<Type> Ns, matrix<Type> xout, vector<Type> LenBins,
                      int maxage, int ngtg, int Nbins, matrix<int> interp_check, matrix<int> interp_check2,
                      matrix<int> integ_check, vector<vector<int> > integ_index) {
  matrix<Type> probLA(maxage, Nbins);
  matrix<Type> yout(xout.rows(), xout.cols());

  for(int a=0; a<maxage; a++) {
    yout.row(a) = linear_int(LAA, Ns, xout, a, interp_check, interp_check2);

    for(int i=0;i<Nbins;i++) {
      if(integ_check(a, i)) {
        probLA(a, i) = yout(a, i);
      } else {
        probLA(a, i) = calcprob_internal(xout, yout, integ_index, a, i, Nbins);
      }
    }
    probLA.row(a) /= probLA.row(a).sum();
  }

  return probLA;
}



//s_dnormal
template <class Type>
vector<Type> s_dnormal(vector<Type> Lengths, Type LFS, Type sl, Type sr, Type Vmaxlen) {
  vector<Type> sel(Lengths.size());
  for(int i=0;i<Lengths.size();i++) {
    Type lo = pow(2,-((Lengths(i) - LFS)/sl*(Lengths(i)-LFS)/sl));
    Type hi = pow(2,-((Lengths(i) - LFS)/sr*(Lengths(i)-LFS)/sr));
    Type hi2 = CppAD::CondExpGe(Vmaxlen, Type(0.99), Type(1), hi);
    sel(i) = CppAD::CondExpLe(Lengths(i), LFS, lo, hi2);
  }
  Type sel_max = max(sel);
  sel /= sel_max;
  return sel;
}

template <class Type>
matrix<Type> s_dnormal(matrix<Type> Lengths, Type LFS, Type sl, Type sr, Type Vmaxlen) {
  matrix<Type> sel(Lengths.rows(), Lengths.cols());
  for(int i=0;i<Lengths.rows();i++) {
    vector<Type> L2 = Lengths.row(i);
    sel.row(i) = s_dnormal(L2, LFS, sl, sr, Vmaxlen);
  }
  return sel;
}

// F     - vector of length n_y
// rdist - vector of length ngtg
// R     - vector of length n_y
// M     - vector of length maxage
// SAA   - matrix of maxage rows, ngtg cols
template <class Type>
matrix<Type> LeesApp_fn(vector<Type> F, Type F_eq, vector<Type> rdist, vector<Type> M, matrix<Type> SAA,
                        vector<Type> LenBins, matrix<Type> LAA, matrix<Type> xout, vector<Type> Select_at_length,
                        matrix<Type> &Select_at_age, vector<matrix<Type> > &NPR, vector<matrix<Type> > &probGTGA,
                        int Nbins, int maxage, int ngtg, int y, matrix<int> interp_check, matrix<int> interp_check2,
                        matrix<int> integ_check, vector<vector<int> > integ_index) {
  matrix<Type> Ns(maxage, ngtg);
  matrix<Type> probGTG(maxage, ngtg);
  Ns.row(0) = rdist;
  probGTG.row(0) = rdist;

  for (int a=1; a<maxage; a++) {
    int yr_st = y-a;
    for(int g=0; g<ngtg; g++) {
      Type Zs = 0;
      for (int a2=0; a2<a; a2++) {
        if(yr_st + a2 < 0) {
          Zs += M(a2) + F_eq * SAA(a2,g);
        } else {
          Zs += M(a2) + F(yr_st + a2) * SAA(a2,g);
        }
      }
      Ns(a,g) = Ns(0,g) * exp(-Zs);
    }
    probGTG.row(a) = Ns.row(a);
    probGTG.row(a) /= Ns.row(a).sum();
  }
  NPR(y) = Ns;
  probGTGA(y) = probGTG;

  // Calculate prob L|A
  matrix<Type> probLA = calcprob(LAA, Ns, xout, LenBins, maxage, ngtg, Nbins, interp_check, interp_check2,
                                 integ_check, integ_index);

  for(int a=0; a<maxage; a++) {
    for(int len=0; len<Nbins; len++) Select_at_age(y, a) += probLA(a, len) * Select_at_length(len);
  }
  return probLA;
}


template <class Type>
matrix<Type> LeesApp_fn(Type F, vector<Type> rdist, vector<Type> M, matrix<Type> SAA, int maxage, int ngtg) {
  matrix<Type> Ns(maxage, ngtg);
  Ns.row(0) = rdist;

  for (int a=1; a<maxage; a++) {
    for(int g=0; g<ngtg; g++) {
      Type Zs = 0;
      for (int a2=0; a2<a; a2++) Zs += M(a2) + F * SAA(a2,g);
      Ns(a,g) = Ns(0,g) * exp(-Zs);
    }
  }

  return Ns;
}


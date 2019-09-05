


template <class Type>
Type log2(Type x) {
  return log(x)/log(2);
}


// Find the conditional maximum and conditional minimum of x
// x1 is the max of x where x < v
// x2 is the min of x where x > v
// Then get the y1 and y2 that correspond to x1 and x2 and interpolate
template <class Type>
Type int_point_internal(Type v, vector<Type> x, vector<Type> y) {
  vector<Type> x_low(x.size());
  vector<Type> x_high(x.size());

  Type vl = x(x.size()-1)+1;
  for(int i=0;i<x.size();i++) {
    x_low(i) = CppAD::CondExpLt(x(i), v, x(i), Type(0));
    x_high(i) = CppAD::CondExpGt(x(i), v, x(i), vl);
  }

  Type x1 = max(x_low);
  Type x2 = min(x_high);
  Type y1 = 0;
  Type y2 = 0;
  for(int i=0;i<x.size();i++) {
    y1 = CppAD::CondExpEq(x1, x(i), y(i), y1);
    y2 = CppAD::CondExpEq(x2, x(i), y(i), y2);
  }

  Type m = y2 - y1;
  m /= x2 - x1;
  return m * (v - x1) + y1;
}


// First, check if v (the i-th xout entry) == any of x (LAA), indices of x < v (x_low), or x > v (x_high)
// The function returns: (1) the correspoding y if v = any of x, else
// (2) 0 if all x < v, else (3) 0 if all x > v, else (4) interpolate otherwise
template <class Type>
Type int_point(Type v, vector<Type> x, vector<Type> y) {
  Type y_eq = -1;
  for(int i=0;i<x.size();i++) y_eq = CppAD::CondExpEq(x(i), v, y(i), y_eq);
  Type ans = CppAD::CondExpGt(y_eq, Type(-1), y_eq,
                              CppAD::CondExpLt(v, min(x), Type(0),
                                               CppAD::CondExpGt(v, max(x), Type(0),
                                                                int_point_internal(v, x, y))));
  return ans;
}


// Interpolated abundance at each xout (ordered concatenation of LenBin and LAA for each GTG)
template <class Type>
vector<Type> linear_int(vector<Type> x, vector<Type> y, vector<Type> xout) {
  int nout = xout.size();
  vector<Type> yout(nout);
  for(int i=0;i<nout;i++) yout(i) = int_point(xout(i), x, y);
  return yout;
}


//calcprob_internal
template <class Type>
Type calcprob_internal2(vector<Type> x, vector<Type> y, int i, int j) {
  Type h1 = abs(y(j) - y(i));
  Type by = x(j) - x(i);
  Type h2 = CppAD::CondExpLt(y(j), y(i), y(j), y(i));
  return h2 * by + 0.5 * by * h1;
}

template <class Type>
Type calcprob_internal(vector<Type> x, vector<Type> y, vector<Type> ind) {
  int ind_size = ind.size();

  matrix<Type> grid(ind_size, ind_size);
  grid.setZero();
  Type area = 0;
  for(int i=0;i<ind_size-1;i++) {
    grid(i,i) = ind(i);
    for(int j=i+1;j<ind_size;j++) {
      grid(i,j) = CppAD::CondExpLt(grid(i,j-1), Type(2), ind(i) * (grid(i,j-1) + ind(j)), Type(3));
      area += CppAD::CondExpEq(grid(i,j), Type(2), calcprob_internal2(x, y, i, j), Type(0));
    }
  }
  return area;
}


// Determine if the j-th xout entry (sorted concatenation of length bin and LAA for each GTG)
// is between the i-th and i-th + 1 length bin
// If yes for only one j, return the interpolated value. Else, sum across j.
template <class Type>
vector<Type> calcprob(vector<Type> x, vector<Type> y, vector<Type> xout, vector<Type> LenBins, int nclasses) {
  int xout_size = xout.size();
  vector<Type> yout(xout_size);
  yout = linear_int(x, y, xout);

  vector<Type> Prob(nclasses);
  for(int i=0;i<nclasses;i++) {
    vector<Type> ind(xout_size);
    for(int j=0;j<xout.size();j++) {
      ind(j) = CppAD::CondExpGe(xout(j), LenBins(i),
          CppAD::CondExpLe(xout(j), LenBins(i+1), Type(1), Type(0)), Type(0));
    }
    Prob(i) = CppAD::CondExpEq(ind.sum(), Type(1), yout(i), calcprob_internal(xout, yout, ind));
  }
  return Prob;
}

template <class Type>
matrix<Type> calcprob_wrapper(matrix<Type> LAA, matrix<Type> Ns, matrix<Type> xout, vector<Type> LenBins,
                              int maxage, int ngtg, int Nbins){
  matrix<Type> probLA(maxage, Nbins);

  vector<vector<Type> > tempprobLA(maxage);
  for(int age=0; age<maxage; age++) {
    vector<Type> tempLAA(ngtg);
    tempLAA = LAA.row(age);
    vector<Type> tempNs(ngtg);
    tempNs = Ns.row(age);
    vector<Type> tempxout(xout.cols());
    tempxout = xout.row(age);

    tempprobLA(age) = calcprob(tempLAA, tempNs, tempxout, LenBins, Nbins);
    probLA.row(age) = tempprobLA(age)/tempprobLA(age).sum();
  }
  return probLA;
}

//s_dnormal
template <class Type>
matrix<Type> s_dnormal(matrix<Type> Lengths, Type LFS, Type sl, Type sr) {
  matrix<Type> sel(Lengths.rows(), Lengths.cols());
  for(int i=0;i<Lengths.rows();i++) {
    for(int g=0;g<Lengths.cols();g++) {
      Type lo = pow(2,-((Lengths(i,g) - LFS)/sl*(Lengths(i,g)-LFS)/sl));
      Type hi = pow(2,-((Lengths(i,g) - LFS)/sr*(Lengths(i,g)-LFS)/sr));
      sel(i,g) = CppAD::CondExpLe(Lengths(i,g), LFS, lo, hi);
    }
  }
  return sel;
}

template <class Type>
vector<Type> s_dnormal(vector<Type> Lengths, Type LFS, Type sl, Type sr) {
  vector<Type> sel(Lengths.size());
  for(int i=0;i<Lengths.size();i++) {
    Type lo = pow(2,-((Lengths(i) - LFS)/sl*(Lengths(i)-LFS)/sl));
    Type hi = pow(2,-((Lengths(i) - LFS)/sr*(Lengths(i)-LFS)/sr));
    sel(i) = CppAD::CondExpLe(Lengths(i), LFS, lo, hi);
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
                        matrix<Type> &Select_at_age, matrix<Type> &Len_at_age, matrix<Type> &NAA,
                        int Nbins, int maxage, int ngtg, int y) {

  matrix<Type> Ns(maxage, ngtg);
  Ns.row(0) = rdist;

  for (int age=1; age<maxage; age++) {
    int yr_st = y-age-1;
    vector<Type> Zs(ngtg);
    Zs.setZero();
    for(int g=0; g<ngtg; g++) {
      for (int age2=0; age2<age; age2++) {
        if(yr_st + age2 < 0) {
          Zs(g) += M(age2) + F_eq * SAA(age2,g);
        } else {
          Zs(g) += M(age2) + F(yr_st + age2) * SAA(age2,g);
        }
      }
      Ns(age,g) = Ns(0,g) * exp(-Zs(g));
    }
  }

  // Calculate prob L|A
  matrix<Type> probLA(maxage, Nbins);
  probLA = calcprob_wrapper(LAA, Ns, xout, LenBins, maxage, ngtg, Nbins);

  for(int age=0; age<maxage; age++) {
    for(int g=0; g<ngtg; g++) Len_at_age(y, age) += LAA(age, g) * Ns(age, g);

    NAA(y, age) = Ns.row(age).sum();
    Len_at_age(y, age) /= NAA(y, age);

    for(int len=0; len<Nbins; len++) Select_at_age(y, age) += probLA(age, len) * Select_at_length(len);

  }

  return probLA;
}


template <class Type>
matrix<Type> LeesApp_fn(Type F, vector<Type> rdist, vector<Type> M, matrix<Type> SAA, int maxage, int ngtg) {

  matrix<Type> Ns(maxage, ngtg);
  Ns.row(0) = rdist;

  for (int age=1; age<maxage; age++) {
    vector<Type> Zs(ngtg);
    Zs.setZero();
    for(int g=0; g<ngtg; g++) {
      for (int age2=0; age2<age; age2++) {
        Zs(g) += M(age2) + F * SAA(age2,g);
      }
      Ns(age,g) = Ns(0,g) * exp(-Zs(g));
    }
  }

  vector<Type> NAA(maxage);
  for(int age=0; age<maxage; age++) NAA(age) = Ns.row(age).sum();

  return NAA;
}


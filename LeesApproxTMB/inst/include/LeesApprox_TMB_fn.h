


template <class Type>
Type log2(Type x) {
  return log(x)/log(2);
}


template <class Type>
Type int_point_internal(Type v, vector<Type> x, vector<Type> x_low, vector<Type> x_high, vector<Type> y) {

  // Find the conditional maximum and conditional minimum of x
  // x1 is the max of x among indices where x < v
  // x2 is the min of x among indices where x > v
  Type x1 = 0;
  Type x2 = max(x);
  for(int i=0;i<x.size();i++) {
    x1 = CppAD::CondExpGt(x_low(i), Type(0), CppAD::CondExpGt(x1, x(i), x1, x(i)), x1);
    x2 = CppAD::CondExpGt(x_high(i), Type(0), CppAD::CondExpLt(x2, x(i), x2, x(i)), x2);
  }

  // Get the index of x1 and x2 and get corresponding y1 and y2 from y
  Type ind1 = 1;
  Type ind2 = 1;
  for(int i=0;i<x.size();i++) {
    ind1 *= CppAD::CondExpEq(x1, x(i), Type(i), Type(1));
    ind2 *= CppAD::CondExpEq(x2, x(i), Type(i), Type(1));
  }
  Type y1 = y(CppAD::Integer(ind1));
  Type y2 = y(CppAD::Integer(ind2));

  // Interpolate
  Type m = (y2-y1)/(x2-x1);

  return m * (v - x1) + y1;
}

template <class Type>
Type int_point_internal2(Type v, vector<Type> x, vector<Type> y) {
  Type ind = 1;
  for(int i=0; i<x.size(); i++) ind *= CppAD::CondExpEq(v, x(i), Type(i), Type(1));
  return y(CppAD::Integer(ind));
}


template <class Type>
Type int_point(Type v, vector<Type> x, vector<Type> y) {

  vector<Type> eq(x.size());
  vector<Type> x_low(x.size());
  vector<Type> x_high(x.size());

  // Check if v (the i-th xout entry) == any of x (LAA), indices of x < v (x_low), or x > v (x_high)
  for(int i=0;i<x.size();i++) {
    eq(i) = CppAD::CondExpEq(x(i), v, Type(1), Type(0));
    x_low(i) = CppAD::CondExpLt(x(i), v, Type(1), Type(0)); // sums to zero if all x > v
    x_high(i) = CppAD::CondExpGt(x(i), v, Type(1), Type(0)); // sums to zero if all x < v
  }

  // Return: (1) the correspoding y if v = any of x, else
  // (2) 0 if all x < v, else (3) 0 if all x > v, else (4) interpolate otherwise
  Type ans = CppAD::CondExpGt(eq.sum(), Type(0), int_point_internal2(v, x, y),
                              CppAD::CondExpLt(x_low.sum(), Type(1), Type(0),
                                               CppAD::CondExpLt(x_high.sum(), Type(1), Type(0),
                                                                int_point_internal(v, x, x_low, x_high, y))));
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
  for(int i=0;i<ind_size-1;i++) {
    grid(i,i) = ind(i);
    for(int j=i+1;j<ind_size;j++) {
      grid(i,j) = CppAD::CondExpLt(grid(i,j-1), Type(2), ind(i) * (grid(i,j-1) + ind(j)), Type(3));
    }
  }

  vector<Type> area(ind_size-1);
  area.setZero();
  for(int i=0; i<ind_size-1; i++) {
    for(int j=i+1; j<ind_size; j++) {
      area(i) += CppAD::CondExpEq(grid(i,j), Type(2), calcprob_internal2(x, y, i, j), Type(0));
    }
  }
  return area.sum();
}


template <class Type>
vector<Type> calcprob(vector<Type> x, vector<Type> y, vector<Type> xout, vector<Type> LenBins, int nclasses) {
  // interpolate at xout
  int xout_size = xout.size();
  vector<Type> yout(xout_size);
  yout = linear_int(x, y, xout);

  // Determine if the j-th xout entry (sorted concatenation of length bin and LAA for each GTG)
  // is between the i-th and i-th + 1 length bin
  // If yes for only one j, return the interpolated value. Else, sum across j.
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
    Type Probsum = tempprobLA(age).sum();
    for(int len=0; len<Nbins; len++) probLA(age, len) = tempprobLA(age)(len)/Probsum;
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

//myseq_len
//template <class Type>
//vector<Type> myseq_len(Type start, Type end, int length) {
//  Type incr = end - start;
//  incr /= asDouble(length);
//
//  vector<Type> seq_out(length);
//  seq_out(0) = start;
//  for(int i=1;i<length;i++) seq_out(i) = seq_out(i-1) + incr;
//
//  return seq_out;
//}

//joinVec
//template <class Type>
//vector<Type> joinVec(vector<Type> V1, vector<Type> V2) {
//  int Length = V1.size() + V2.cols();
//  vector<Type> Vout(Length);
//  for(int i=0;i<V1.size();i++) Vout(i) = V1(i);
//  for(int i=0;i<V2.size();i++) Vout(V1.size()+i) = V2(i);
//  return Vout;
//}

//joinVec
//template <class Type>
//vector<Type> joinVec(vector<Type> V1, matrix<Type> V2, int a) {
//  int Length = V1.size() + V2.cols();
//  vector<Type> Vout(Length);
//  for(int i=0;i<V1.size();i++) Vout(i) = V1(i);
//  for(int i=0;i<V2.size();i++) Vout(V1.size()+i) = V2(a,i);
//  return Vout;
//}

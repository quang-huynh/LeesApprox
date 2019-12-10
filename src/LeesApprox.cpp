#include <Rcpp.h>
using namespace Rcpp;

double int_point(double v, NumericVector x, NumericVector y, int opt=1) {
  LogicalVector eq = v == x;

  if (sum(eq)>0) { // return y if v %in% x
    NumericVector temp = y[eq];
    double yout = as<double>(temp);
    return(yout);
  }
  LogicalVector x_low = x < v;
  LogicalVector x_high = x > v;

  if (opt == 1) {
    if (sum(x_low) < 1) return(0);
    if (sum(x_high) < 1) return(0);
  } else {
    if (sum(x_low) < 1) throw std::range_error("v < min(x)");
    if (sum(x_high) < 1) throw std::range_error("v > max(x)");
  }

  NumericVector temp1 = x[x_low];
  double x1 = max(temp1);

  LogicalVector temp1a = x == x1;
  NumericVector y1vec = y[temp1a];
  double y1 = as<double>(y1vec);

  NumericVector temp2 = x[x_high];
  double x2 = min(temp2);
  LogicalVector temp2a = x == x2;
  NumericVector y2vec = y[temp2a];
  double y2 = as<double>(y2vec);

  double m = (y2 - y1)/(x2 - x1);
  double yout = m * (v-x1) + y1;
  return(yout);

}

//' Title
//'
//' Description
//'
//' @param x An integer vector
//' @useDynLib LeesApprox
//' @export
//[[Rcpp::export]]
NumericVector linear_int(NumericVector x, NumericVector y, NumericVector xout) {
  int nout = xout.size();
  NumericVector yout(nout);
  for (int i=0; i<nout; i++) {
    yout[i] = int_point(xout[i], x, y, 1);
  }
  return(yout);
}

//[[Rcpp::export]]
NumericVector calcprob(NumericVector x, NumericVector y, NumericVector xout,
                       NumericVector LenBins) {
  // interpolate at xout
  NumericVector yout = linear_int(x, y, xout);

  int nclasses = LenBins.size() - 1;
  NumericVector Prob(nclasses);

  for (int i=0; i<nclasses; i++) {
    double L1 = LenBins(i);
    double L2 = LenBins(i+1);
    LogicalVector ind = (xout >=L1) & (xout <=L2);
    if (sum(ind)==1) {
      Prob[i] = yout(i);
    } else {
      NumericVector Ys = yout[ind];
      NumericVector Xs = xout[ind];
      int ns = sum(ind)-1;
      NumericVector area(ns);
      for (int j=0; j<ns; j++) {
        double p1 = Ys(j);
        double p2 = Ys(j+1);
        double h1 = std::abs(p1-p2);
        double by = Xs(j+1) - Xs(j);
        NumericVector tvec = NumericVector::create(p1, p2);
        double h2 = min(tvec);
        area[j] = h2 * by + 0.5 * by * h1;
      }
      Prob[i] = sum(area);
    }
  }
  NumericVector Probout = Prob/sum(Prob);
  return(Probout);
}

struct add_multiple {
  int incr;
  int count;
  add_multiple(int incr)
    : incr(incr), count(0)
  {}
  inline int operator()(int d) {
    return d + incr * count++;
  }
};

//[[Rcpp::export]]
Rcpp::NumericVector myseq_by(double x, double y, double by) {
  NumericVector anOut(1);
  // compute sequence
  double min_by = 1.e-8;
  if (by < min_by) min_by = by/100;
  double i = x + by;
  anOut(0) = x;
  while(i/min_by < y/min_by + 1) {
    anOut.push_back(i);
    i += by;
  }
  return anOut;
}


//' Title
//'
//' Description
//'
//' @param x An integer vector
//' @useDynLib LeesApprox
//' @export
//[[Rcpp::export]]
NumericVector myseq_len(double start, double end, double length) {
  double by = end/(length-1) * 2;
  NumericVector seq_out = myseq_by(start, end, by);
  return(seq_out);
}

NumericVector s_dnormal(NumericVector Lengths, double LFS, double sl, double sr) {
  NumericVector sel(Lengths.length());
  for (int i=0; i<Lengths.length(); i++) {
    if(Lengths(i)<=LFS) {
      sel(i) = pow(2,-((Lengths(i) - LFS)/sl*(Lengths(i)-LFS)/sl));
    } else {
      sel(i) = pow(2,-((Lengths(i) - LFS)/sr*(Lengths(i)-LFS)/sr));
    }
  }
  return(sel);
}

NumericVector joinVec(NumericVector V1, NumericVector V2) {
  int Length = V1.length() + V2.length();
  NumericVector Vout(Length);
  for (int i=0; i<V1.length(); i++) {
    Vout(i) = V1(i);
  }
  for (int i=0; i<V2.length(); i++) {
    Vout(i+V1.length()) = V2(i);
  }
  return(Vout);
}

//' Title
//'
//' Description
//'
//' @param x An integer vector
//' @export
//[[Rcpp::export]]
List LeesApprox(NumericVector FVec, int ngtg, double maxsd, double binwidth,
            double M, double Linf, double K, double t0,
            double LFS, double L5, double Vmaxlen,
            double LinfCV, int maxage) {

  NumericVector distGTG = myseq_len(-maxsd, maxsd, ngtg); // sd distribution of Linf for GTGs
  NumericVector rdist = dnorm(distGTG, 0, 1.0, 0); // dist of recruits across GTGs
  rdist = rdist/sum(rdist); // sum to 1
  NumericVector Linfgtg = Linf + LinfCV*Linf*distGTG; // dist of Linf across GTGs

  NumericVector LenBins = myseq_by(0, max(Linfgtg), binwidth);
  NumericVector LenMids = LenBins - 0.5 * binwidth;
  LenMids = LenMids[LenMids>0];

  // calculate length-at-age for each GTG
  NumericMatrix LAA(maxage, ngtg);
  NumericVector ages = myseq_by(1.0, maxage, 1);

  for (int g=0; g<ngtg; g++) {
    LAA(_,g) = Linfgtg(g) * (1-exp(-K * (ages-t0)));
  }

  // calculate selectivity-at-age for each GTG
  double sl = (LFS - L5) / pow(-log2(0.05),0.5);
  double sr = (Linf - LFS) / pow(-log2(Vmaxlen),0.5);

  NumericMatrix SAA(maxage, ngtg); // selectivity-at-age by GTG
  for (int g=0; g<ngtg; g++) {
    SAA(_,g) = s_dnormal(LAA(_,g), LFS, sl, sr);
  }

  // calculate selectivity-at-length
  NumericVector Select_at_length = s_dnormal(LenMids, LFS, sl, sr);

  // create vector of Fs for last maxage years
  // F assumed 0 for years where no F are provided
  NumericVector FVec2(maxage);
  int nFs = FVec.length();
  int ind = nFs-1;
  for (int i=maxage-1; i>=0; i--) {
    if (ind >= 0) {
      FVec2(i) = FVec(ind);
    }
    ind -= 1;
  }

  // loop over ages and calculate N per recruit for each age class
  NumericMatrix Ns(maxage, ngtg);
  Ns(0,_) = rdist;
  for (int age=1; age<maxage; age++) {
    int yr_st = maxage-age-1;
    int age_end = age;
    NumericVector Zs(ngtg);
    for (int age2=0; age2<age_end; age2++) {
      Zs += M + FVec2(yr_st + age2) * SAA(age2,_);
    }
    Ns(age,_) = Ns(0,_) * Rcpp::exp(-Zs);
  }

  // Calculate prob L|A
  int Nbins = LenBins.length() - 1;
  NumericMatrix probLA(maxage, Nbins);
  for (int age=0; age<maxage; age++) {
    NumericVector tempLens = LAA(age,_);
    NumericVector tempNs = Ns(age,_);
    NumericVector xout = joinVec(LenBins, tempLens);
    std::sort(xout.begin(), xout.end());
    probLA(age,_) = calcprob(tempLens, tempNs, xout, LenBins);
  }

  // Calculate mean selectivity-at-age
  NumericVector Select_at_age(maxage);
  for (int age=0; age<maxage; age++) {
    Select_at_age(age) += sum(probLA(age,_) * Select_at_length);
  }

  // Calculate mean length-at-age
  NumericVector Len_at_age(maxage);
  for (int age=0; age<maxage; age++) {
    Len_at_age(age) = sum(LAA(age,_)*Ns(age,_))/(sum(Ns(age,_)));
  }

  // Calculate fished length composition
  NumericVector NAA(maxage);
  for (int age=0; age<maxage; age++) {
    NAA(age) = sum(Ns(age,_));
  }

  NumericVector LenComp(Nbins);
  for (int l=0; l<Nbins; l++) {
    LenComp(l) = sum(probLA(_,l) * Select_at_length(l) * NAA/sum(NAA));
  }
  LenComp = LenComp/sum(LenComp);

  List out(8);
  out(0) = probLA;
  out(1) = LenComp;
  out(2) = Select_at_age;
  out(3) = Select_at_length;
  out(4) = LAA;
  out(5) = LenMids;
  out(6) = LenBins;
  out(7) = NAA;
  return(out);
}



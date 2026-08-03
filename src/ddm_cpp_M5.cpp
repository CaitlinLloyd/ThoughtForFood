// Time-varying DDM simulation with health delay relative to umami and taste (3-factor model)
// Modified from original code adapted from Hsiang-Yu and Todd Hare (https://osf.io/k8g39/)

#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <time.h>

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator

// [[Rcpp::depends(RcppParallel)]]
struct ddm3hw : public Worker {

  // input vectors/matrices to read from
  const double d_v;        // taste drift coefficient
  const double d_h;        // health drift coefficient  
  const double d_p;        // umami drift coefficient
  const double thres;      // decision boundary
  const double nDT;        // non-decision time
  const double bias;       // starting bias
  const double vd;         // taste value
  const double hd;         // health value
  const double pd;         // umami value
  const double sd_n;       // noise standard deviation
  const double tIn_h;      // health timing relative to umami and taste (from Bayesian model)
  GEN gen;

  // output vector to write to
  RVector<double> vecOut;

  // initialize from Rcpp input and output matrices/vectors
  ddm3hw(const double d_v, const double d_h, const double d_p, const double thres, const double nDT,
        const double bias, const double vd, const double hd, const double pd, const double sd_n,
        const double tIn_h, NumericVector vecOut, GEN gen)
    : d_v(d_v), d_h(d_h), d_p(d_p), thres(thres), nDT(nDT), bias(bias), vd(vd), hd(hd), pd(pd),
      sd_n(sd_n), tIn_h(tIn_h), gen(gen), vecOut(vecOut) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    double T = 5.2, dt = 0.001, lt;
    lt = (int)(T/dt);

    // Initialize timing vectors - umami and taste are always active from t=0
    std::vector<double> vec_tTaste(lt, 1);    // taste is always active from t=0
    std::vector<double> vec_tumami(lt, 1);  // umami is always active from t=0
    std::vector<double> vec_tHealth(lt, 0);   // health starts inactive

    // Determine timing based on health onset parameter
    double abs_tIn_h = std::abs(tIn_h);
    
    // Set activation times
    int start_h = (int)(abs_tIn_h / dt);  // when health becomes active
    
    // Handle different timing scenarios
    if (tIn_h >= 0) {
      // Health comes after umami and taste
      for (int t = start_h; t < lt; t++) {
        vec_tHealth[t] = 1;
      }
    } else {
      // Health comes before umami and taste - umami and taste are delayed
      for (int t = 0; t < lt; t++) {
        vec_tHealth[t] = 1;  // health active from start
      }
      for (int t = 0; t < start_h && t < lt; t++) {
        vec_tTaste[t] = 0;   // taste delayed
        vec_tumami[t] = 0; // umami delayed
      }
    }

    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = T;
      double X = bias * thres;
      int flag = 0;
      double cont = 0;
      double noise = 0;

      while (flag == 0 && cont < lt) {
        
        int time_idx = (int)cont;
        
        // Apply time-varying weights based on activation vectors
        double taste_weight = vec_tTaste[time_idx] * d_v;
        double health_weight = vec_tHealth[time_idx] * d_h;
        double umami_weight = vec_tumami[time_idx] * d_p;

        noise = gen() * sqrt(dt);
        X = X + (taste_weight*vd + health_weight*hd + umami_weight*pd) * dt + noise;

        if (X > thres) {
          flag = 1;
          vecOut[i] = nDT + cont * dt;
        }
        else if (X < 0) {
          flag = 1;
          vecOut[i] = -nDT - cont * dt;
        }
        cont++;
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector ddm3h_parallel(double d_v, double d_h, double d_p, double thres, double nDT, double bias, 
                           double vd, double hd, double pd, double sd_n, double tIn_h, unsigned int N) {

  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0, sd_n);
  GEN  gen(eng, dist);

  // output vector
  NumericVector vecOut(N);

  // create the worker
  ddm3hw ddm3hw(d_v, d_h, d_p, thres, nDT, bias, vd, hd, pd, sd_n, tIn_h, vecOut, gen);

  // call the worker
  parallelFor(0, N, ddm3hw);

  return vecOut;
}
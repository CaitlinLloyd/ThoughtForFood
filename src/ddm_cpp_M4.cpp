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
struct ddm3w : public Worker {

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
  const double tIn_p;      // umami timing relative to health and taste (from Bayesian model)
  GEN gen;

  // output vector to write to
  RVector<double> vecOut;

  // initialize from Rcpp input and output matrices/vectors
  ddm3w(const double d_v, const double d_h, const double d_p, const double thres, const double nDT,
        const double bias, const double vd, const double hd, const double pd, const double sd_n,
        const double tIn_p, NumericVector vecOut, GEN gen)
    : d_v(d_v), d_h(d_h), d_p(d_p), thres(thres), nDT(nDT), bias(bias), vd(vd), hd(hd), pd(pd),
      sd_n(sd_n), tIn_p(tIn_p), gen(gen), vecOut(vecOut) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    double T = 5.2, dt = 0.001, lt;
    lt = (int)(T/dt);

    // Initialize timing vectors - health and taste are always active from t=0
    std::vector<double> vec_tTaste(lt, 1);    // taste is always active from t=0
    std::vector<double> vec_tHealth(lt, 1);   // health is always active from t=0
    std::vector<double> vec_tUmami(lt, 0);  // umami starts inactive

    // Determine timing based on umami onset parameter
    double abs_tIn_p = std::abs(tIn_p);
    
    // Set activation times
    int start_p = (int)(abs_tIn_p / dt);  // when umami becomes active
    
    // Handle different timing scenarios
    if (tIn_p >= 0) {
      // Umami comes after health and taste
      for (int t = start_p; t < lt; t++) {
        vec_tUmami[t] = 1;
      }
    } else {
      // Umami comes before health and taste - health and taste are delayed
      for (int t = 0; t < lt; t++) {
        vec_tUmami[t] = 1;  // umami active from start
      }
      for (int t = 0; t < start_p && t < lt; t++) {
        vec_tTaste[t] = 0;   // taste delayed
        vec_tHealth[t] = 0;  // health delayed
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
        double umami_weight = vec_tUmami[time_idx] * d_p;

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
NumericVector ddm3_parallel(double d_v, double d_h, double d_p, double thres, double nDT, double bias, 
                           double vd, double hd, double pd, double sd_n, double tIn_p, unsigned int N) {

  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0, sd_n);
  GEN  gen(eng, dist);

  // output vector
  NumericVector vecOut(N);

  // create the worker
  ddm3w ddm3w(d_v, d_h, d_p, thres, nDT, bias, vd, hd, pd, sd_n, tIn_p, vecOut, gen);

  // call the worker
  parallelFor(0, N, ddm3w);

  return vecOut;
}
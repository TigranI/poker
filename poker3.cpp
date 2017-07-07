#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "poker.hpp"

#include <nlopt.hpp>

//double optFn(const std::vector<double> &x, std::vector<double> &grad, void *data)
//{
//  objfunc_data *d = reinterpret_cast<objfunc_data*>(data);
//  double sum_singular = 0;
//  uint64_t k = d->k;

//  for(uint64_t i = 0; i < k; ++i) {
//      sum_singular += x[i];
//  }
//  return sum_singular;
//}


//double optCstr(const std::vector<double> &x, std::vector<double> &grad, void *data)
//{
//    constraint_data *d = reinterpret_cast<constraint_data*>(data);
//    double r = d->r; // entry r_ij
//    uint64_t row = d->row;
//    uint64_t col = d->col;
//    uint64_t k = d->k;
//    uint64_t m = d->m;

//    double constraint = 0;
//    for(uint64_t t = 0; t < k; ++t) {
//        constraint += x[t] * x[k + row*k + t] * x[k + m*k + col*k + t];
//    }
//    constraint -= r;
//    return constraint;
//}

void Game::playZ() {

  uint64_t P = 10;
  nlopt::opt opt(nlopt::LN_COBYLA, P);

  std::vector<double> lb(P);
  std::vector<double> ub(P);

  for(uint64_t q = 0; q < P/2; ++q) {
      lb[q] = -1;
      ub[q] = 1;
  }
  for(uint64_t q = P/2; q < P; ++q) {
      lb[q] = -1;
      ub[q] = 1;
  }

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

 // opt.set_min_objective(optFn, &fdata);

//    std::vector<constraint_data> cdata(pointsT.size() + rank*(rank+1));
//    uint64_t q = 0;
//    for( auto it = pointsT.begin(); it != pointsT.end(); ++it) {
//        cdata[q].row = (it->first).first;
//        cdata[q].col = (it->first).second;
//        cdata[q].r = it->second;
//        cdata[q].k = rank;
//        cdata[q].m = M;
//        cdata[q].n = N;
//        opt.add_equality_constraint(optCstr, &cdata[q++], 1e-8);
//    }

  opt.set_xtol_rel(1e-4);

  // initial
 // std::vector<double> x0 = setInitial();
//    std::vector<double> x(P);
//    for(uint64_t q = 0; q < P; ++q) {
//        x0[q] = unifd(gen);
//    }

  double minf = 0;
 // nlopt::result result = opt.optimize(x0, minf);

}

// computes ev(Y) = -ev(X)
double computeEV(double x1, double y1, double y2, double initialpot, double calledpot) {

  // check,check + bet,fold + bet,call
  double evcheck = - (y2 - y1) * (1 - y2) * initialpot / 2; // x>y2, y wins
  evcheck += (y2 - y1) * y1 * initialpot / 2;   // x<y1 , x wins

  double evcall = 0;
  if(x1 <= y1) {
      evcall += (y1 - x1) * x1 * calledpot / 2; // x<x1 & x1<y<y1, x wins
      evcall += (1 - y2) * x1 * calledpot / 2;  // x<x1 & y>y2  , x wins
    } else if(x1 <= y2) {
      evcall -= (x1 - y1) * y1 * calledpot / 2; // y1<x<x1 & y<y1, y wins
      evcall += (1 - y2) * x1 * calledpot / 2;   // y>y2 & x<x1 , x wins
    } else { // x1 > y2
      evcall += (1 - x1) * x1 * calledpot / 2;  // x<x1 & y>(1-x1) , x wins
      evcall -= y1 * (x1 - y1) * calledpot / 2; // <y1<x<x1 & y<y1 , y wins
    }

  double evfold = - (y1 + (1 - y2)) * (1 - x1) * initialpot / 2; // (y<y1 or y > y2) & x > x1 , y wins

  double ev = evcheck + evcall + evfold;  // ev(Y)

  return ev;
}

void Game::playAA() {

  double bet = 1;
  double initialpot = 2;
  double calledpot = 2 + 2 * bet;
  double step = 0.01;

  double maxminEV = -std::numeric_limits<double>::infinity();
  double y1bb = 0;
  double y2bb = 0; // best param set for player 1,2
  double x1bb = 0;

  // EV from X perspective: max min ev(Y)
  for(double x1 = 0.0; x1 <= 1.0; x1 += step) {  // x1
      double minEV = std::numeric_limits<double>::infinity(); // max of min EV is best for player 1
      double y1b = 0;
      double y2b = 0; // best param set for player 1,2
   for(double y1 = 0.0; y1 <= 0.7; y1 += step) {  // y1
      for(double y2 = y1; y2 <= 1.0; y2 += step) {  // y2

          double ev = computeEV(x1, y1, y2, initialpot, calledpot); // ev(Y)

          if(ev < minEV) {
              y1b = y1;
              y2b = y2;
              minEV = ev;
            }
        }  // y2
     }  // y1
   if(minEV > maxminEV) {
       y1bb = y1b;
       y2bb = y2b;
       x1bb = x1;
       maxminEV = minEV;
     }
   }  // x1

  std::cout << "MAX-MIN:  " << x1bb << "  " << y1bb << "  " << y2bb << "    EV =  " << maxminEV << "\n";


  double minmaxEV = std::numeric_limits<double>::infinity();

  // EV from Y perspective: min max ev(X)
  for(double y1 = 0.0; y1 <= 0.7; y1 += step) {  // y1
     for(double y2 = y1; y2 <= 1.0; y2 += step) {  // y2
      double maxEV = -std::numeric_limits<double>::infinity(); // max of min EV is best for player 1
      double x1b = 0; // best param set for player 1,2
      for(double x1 = 0.0; x1 <= 1.0; x1 += step) {  // x1

          double ev = computeEV(x1, y1, y2, initialpot, calledpot); // ev(Y)

          if(ev > maxEV) {
              x1b = x1;
              maxEV = ev;
            }
     }  // x1
   if(maxEV < minmaxEV) {
       y1bb = y1;
       y2bb = y2;
       x1bb = x1b;
       minmaxEV = maxEV;
     }
   }    // y2
  }   // y1

  std::cout << "MIN-MAX:  " << x1bb << "  " << y1bb << "  " << y2bb << "    EV =  " << minmaxEV << "\n";

  double x1o = 0.6;
  double y1o = 0.3;
  double y2o = 1;
  double evOpt = computeEV(x1o, y1o, y2o, initialpot, calledpot); // ev(Y)

  std::cout << "OPT?:  " << x1o << "  " << y1o << "  " << y2o << "    EV =  " << evOpt << "\n";

  double dz = 0.01;
  double eps = 0.000001;
  for(double y1 = 0.0; y1 <= 0.7; y1 += step) {  // y1
     for(double y2 = y1; y2 <= 1.0; y2 += step) {  // y2
      for(double x1 = 0.0; x1 <= 1.0; x1 += step) {  // x1

          double ev = computeEV(x1, y1, y2, initialpot, calledpot); // ev(Y)
          bool flag = false;

              for(double e2 = -dz; e2 <= dz+0.001; e2 += dz) {
                  for(double e3 = -dz; e3 <= dz+0.001; e3 += dz) {
                      if(e2 == 0 && e3 == 0)
                        continue;
                      double ev1 = computeEV(x1, y1+e2, y2+e3, initialpot, calledpot); // ev(Y)
                      if(ev1 < ev - eps) {
                          flag = true;
                          break;
                        }
                    }
                }
            for(double e1 = -dz; e1 <= dz+0.001; e1 += dz) {
                if(e1 == 0)
                  continue;
                double ev1 = computeEV(x1+e1, y1, y2, initialpot, calledpot); // ev(Y)
                if(ev1 > ev + eps) {
                    flag = true;
                    break;
                  }
              }
            if(!flag) {
                std::cout << "OPT!:  " << x1 << "  " << y1 << "  " << y2 << "    EV =  " << ev << "\n";
              }
        }
      }
    }

}

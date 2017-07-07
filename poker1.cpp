#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "poker.hpp"

std::vector<double> Game::computePB(double x1, double x2, double y1, double y2, double y3) {

  double ptotal = 0;
  double estep = 0.01;
  std::vector<double> pb;

  for(double x = 0; x <= 1 + 0.001; x += estep) {   // equity

      double p = 0;
      if(x >= x2) {
          p = (y2-y3)/(1-x2) * x + y2 - (y2-y3)/(1-x2);
      }
      else if(x >= x1)
          p = 0;
      else
          p = (y3-y1)/x1 * x + y1;
      ptotal += p;
      pb.push_back(p);
     // std::cout << x << "  " << db << "  " << prob << "\n";
    }

  for(int i=0; i < pb.size(); ++i) {
      pb[i] /= ptotal;
     // std::cout << pb[i] << "\n";
    }

  //  double sum = 0;
  //  for(double p : pb) {
  //      sum += p;
  //    }
  //  std::cout << sum << "\n";

  return pb;
}

bool Game::betStrategyC(double x1, double x2, double y1, double y2, double y3, double h1, double u) {

  double p = 0; // prob of bet
  if(h1 >= x2) {
      p = (y2-y3)/(1-x2) * h1 + y2 - (y2-y3)/(1-x2);
  }
  else if(h1 >= x1)
      p = 0;
  else
      p = (y3-y1)/x1 * h1 + y1;

  if(u < p)
       return true;
      else
       return false;

}

bool Game::callStrategyC(std::vector<double> pb, double h2, double odds) {

   double estep = 0.01;  // must be the same as in computePB()
   double hcdf = 0;
   int i = 0;
   int len = pb.size();
   double x = 0;
   while(x <= h2 && i < len) {
       hcdf += pb[i++];
       x += estep;
     }

//std::cout << " in callStrategyB() : "  <<  hcdf << "    " << odds << "  " << h << "\n";

   if(hcdf >= odds)
     return true;
   else
     return false;

}

void Game::playC() {

  std::mt19937 gen;
  gen.seed(7);
  std::uniform_real_distribution<double> unifh(0, 1);
  std::uniform_real_distribution<double> unifp(0, 1);
  double xstep = 0.1;
  //int I = 20000;
  double bet = 1;
  double pot = 2 + 2 * bet;

//  double minEV = std::numeric_limits<double>::infinity();
    double maxEV = -std::numeric_limits<double>::infinity();
    double x1b = 0;
    double x2b = 0; // best param set
    double y1b = 0;
    double y2b = 0;
    double y3b = 0;
    double betfreqb = 0;

  // param list (x1, x2, y1)  : p = step function
  for(double x1 = 0.5; x1 <= 0.5; x1 += xstep) {  // x1
      double x2 = x1;
     // for(double x2 = x1 + xstep; x2 <= 0.7; x2 += xstep) {  // x2
          for(double y1 = 0.5; y1 <= 0.8; y1 += xstep) {  // y1
              for(double y2 = 1.0; y2 <= 2.0; y2 += xstep) {  // y2
               for(double y3 = 0; y3 <= 0.5; y3 += xstep) {  // y3 at x2
                  double ev = 0;
                  double betfreq = 0;
                  std::vector<double> pb = computePB(x1, x2, y1, y2, y3);

                  //for(int iter = 0; iter < I; ++iter) {
                  int iter = 0;
                   for(double h1 = 0; h1 <= 1.0; h1 += 0.01) {  // h1
                       for(double h2x = 0; h2x <= 1.0; h2x += 0.01) { // h2x
                      ++iter;
//                      double h1 = unifh(gen);
//                      double h2 = unifh(gen);
                           double s = unifh(gen);
                           double h2 = h2x + sgnf(s - 0.5) * 0.001; // perturb
                      double u = unifp(gen);
                      if(betStrategyC(x1, x2, y1, y2, y3, h1, u)) {
                         // std::cout << h1 << "  " << h2 << "\n";
                          betfreq += 1;
                          double odds = bet / pot;
                          if(callStrategyC(pb, h2, odds)) {
                              ev += sgnf(h1 - h2) * pot / 2.0;  // bet, call
                            } else {
                              ev += 1.0;   // bet, fold
                            }
                        } else {
                          ev += sgnf(h1 - h2);  // no bet
                        }
                  //  } // I
                   }  // h2
                 }   // h1

                  betfreq /= iter;
                  ev /= iter;
                  std::cout << x1 << "  " << x2 << "   " << y1 << "  " << y2 << "  " << y3 << "  " << betfreq << "   EV =  " << ev << "\n";
//                  if(ev < minEV) {
//                      std::cout << x1 << "  " << x2 << "   " << p1 << "  " << p2 << "    EV =  " << ev << "\n";
//                      minEV = ev;
//                  }

                  if(ev > maxEV) {
                     // std::cout << x1 << "  " << x2 << "   " << p1 << "  " << p2 << "    EV =  " << ev << "\n";
                      maxEV = ev;
                      x1b = x2;
                      x2b = x2;
                      y1b = y1;
                      y2b = y2;
                      y3b = y3;
                      betfreqb = betfreq;
                    }
                 } // y3
                } // y2
                } // y1
          //  } // x2
        }  // x1

      std::cout << "best param set: " << x1b << "  " << x2b << "   " << y1b << "  " << y2b << "  " << y3b << "  " << betfreqb << "  EV =  " << maxEV << "\n";

}


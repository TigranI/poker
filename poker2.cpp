#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "poker.hpp"

std::vector<double> Game::pdfBet(std::vector<double> coefs) {

  double estep = 0.01;
  double ptotal = 0;
  std::vector<double> pdfb;

  for(double x = 0; x <= 1 + 0.001; x += estep) {   // equity

      double p = 0;
      if(x >= 0.5) {
          p = coefs[2] * x + coefs[3];
      }
      else
        p = coefs[0] * x + coefs[1];

      ptotal += p;
      pdfb.push_back(p);
     // std::cout << x << "  " << db << "  " << prob << "\n";
    }

  for(int i=0; i < pdfb.size(); ++i) {
      pdfb[i] /= ptotal;
    }

  //  double sum = 0;
  //  for(double p : pb) {
  //      sum += p;
  //    }
  //  std::cout << sum << "\n";

  return pdfb;

}

std::vector<double> Game::pdfCheck(std::vector<double> pdfb) {

  int N = pdfb.size();
  std::vector<double> pdfch(N);
  double ptotal = 0;

  for(int i = 0; i < N; ++i) {
      pdfch[i] = 1 - pdfb[i];
      ptotal += pdfch[i];
    }

  for(int i=0; i < pdfch.size(); ++i) {
      pdfch[i] /= ptotal;
    }

  return pdfch;

}

// cdf from pdf
std::vector<double> Game::cdf(std::vector<double> pdf) {

  int N = pdf.size();
  std::vector<double> cdf(N);
  double ptotal = 0;
  for(int i = 0; i < N; ++i) {
      ptotal += pdf[i];
      cdf[i] = ptotal;
    }

  return cdf;
}

// coefs from conditions: (a1, b1, a2, b2). y = ax + b
std::vector<double>  computeCoefs(double y1, double y2) {

  std::vector<double> coefs(4);
  coefs[0] = -2 * y1;
  coefs[1] = y1;
      coefs[2] = 2 * y2;
  coefs[3] = -y2;

  return coefs;
}

double computeProb(double h, double x1, double y1, double x2, double y2) {

  double eps = 0.01;

  if(x2 - x1 < eps)
    return 0;

  double p = 0;
  double mid = (x1 + x2) / 2;

  if(h < mid) {
      double t = (mid - h) / (mid - x1);
      p = t * y1;
    } else {
      double t = (h - mid) / (x2 - mid);
      p = t * y2;
    }

  return p;
}

// computes weakest calling hand
double Game::computeCallThreshold(double odds, std::vector<double> cdfBettor) {

  double estep = 0.01;
  double h = 0;
  int i = 0;
  int len = cdfBettor.size();
  while(cdfBettor[i] < odds && i < len) {
      ++i;
      h += estep;
    }

  return h;
}

bool Game::betStrategyT(double h, double u, double x1, double y1, double x2, double y2) {

//  double p = 0;
//  if(h >= 0.5)
//      p = coefs[2] * h + coefs[3];
//  else
//    p = coefs[0] * h + coefs[1];

  double p = computeProb(h, x1, y1, x2, y2);

  if(u < p)
    return true;
  else
    return false;

}

bool Game::callStrategyT(double hCallTh, double h) {

//  int i = 0;
//  int len = cdfBettor.size();
//  while(cdfBettor[i] < h && i < len) {
//      ++i;
//    }

//  if(cdfBettor[i] >= odds)
//    return true;
//  else
//    return false;

  if(h >= hCallTh)
    return true;
  else
    return false;
}

bool Game::raiseStrategyT(double h, double u, double x1, double y1, double x2, double y2) {

//  return false;
  double p = computeProb(h, x1, y1, x2, y2);

  if(u < p)
    return true;
  else
    return false;

}

// MIN-MAX:  0.4  1.3  0.1  1.1   EV =  -0.04235
void Game::playT() {

  // initial hand distributions for p1 and p2
  int N = 100;
  std::vector<double> pdf1(N, 1.0/N);
  std::vector<double> cdf1 = cdf(pdf1);
  std::vector<double> pdf2(N, 1.0/N);
  std::vector<double> cdf2 = cdf(pdf2);

  std::mt19937 genh;
  genh.seed(53);
  std::mt19937 genp;
  genp.seed(15);
  std::uniform_real_distribution<double> unifh(0, 1);
  std::uniform_real_distribution<double> unifp(0, 1);

  int I = 50000;
  double bet = 1;
  double raisebet = 3;
  double pot = 2 + 2 * bet;
  double raisedpot = 2 + 2 * raisebet;

  double maxminEV = -std::numeric_limits<double>::infinity(); // max of min EV is best for player 1
  double y1bb = 0;
  double y2bb = 0; // best param set for player 1
  double z1bb = 0;
  double z2bb = 0;
  double w1bb = 0;
  double w2bb = 0;
  //double bet1freqb = 0;
  double step = 0.1;

  for(double y1 = 0.4; y1 <= 0.4; y1 += step) {  // y1
      for(double y2 = 1.6; y2 <= 1.6; y2 += step) {  // y2
          for(double w1 = 0.0; w1 <= 0.5; w1 += step) {  // w1
              for(double w2 = 0; w2 <= 1; w2 += step)  {  // w2
          double minEV = std::numeric_limits<double>::infinity(); // max EV of player 2
          double y1b = 0;
          double y2b = 0; // best param set for player 2
          double z1b = 0;
          double z2b = 0;
          double w1b = 0;
          double w2b = 0;

          std::vector<double> coefs1 = computeCoefs(y1, y2);

          std::vector<double> pdfBet1 = pdfBet(coefs1);
          std::vector<double> cdfBet1 = cdf(pdfBet1);
          std::vector<double> pdfCheck1 = pdfCheck(pdfBet1); // for raise
          std::vector<double> cdfCheck1 = cdf(pdfCheck1);

//          for(int i=0;i<cdfBet1.size();i++) {
//              std::cout << "pcdfs:  " << pdfBet1[i] << "   " << pdfCheck1[i]  << "   " << cdfBet1[i] << "   " << cdfCheck1[i] << "\n";
//            }


          double odds = bet / pot;
          double rodds = (raisebet - bet) / raisedpot;

          double callTh1 = computeCallThreshold(odds, cdf2);  // odds are the same here for p1 and p2
          double callTh2 = computeCallThreshold(odds, cdfBet1);
          double callTh3 = computeCallThreshold(rodds, cdfCheck1);

          for(double z1 = 0.1; z1 <= 0.1; z1 += step) {  // z1
              for(double z2 = 1.1; z2 <= 1.1; z2 += step) {  // z2

            double EV = 0;

            for(int iter = 0; iter < I; ++iter) {
              double h1 = unifh(genh);
              double h2 = unifh(genh);

              double u1 = unifp(genp);
              if(betStrategyT(h1, u1, 0, y1, 1, y2)) {  // bet by p1
                  if(callStrategyT(callTh2, h2)) { // call by p2
                       EV += sgnf(h1 - h2) * pot / 2.0;  // bet, call
                    } else {
                      EV += 1.0;  // bet, fold
                    }
                } else {
                 // std::vector<double> coefs2 = computeCoefs(z1, z2);
                  double u2 = unifp(genp);
                  if(betStrategyT(h2, u2, 0, z1, 1, z2)) {  // bet by p2
                      if(callStrategyT(callTh1, h1)) {  // call by p1
                          double u3 = unifp(genp);
                          if(raiseStrategyT(h1, u3, callTh3, w1, 1, w2)) {     // raise by p1
                              // std::cout <<  odds << "  " << rodds << "  " << callTh1 << "  " << callTh2 << "  " << callTh3 << "   RAISED! "  << "\n";
                              if(callStrategyT(callTh3, h2)) {   // call by p2
                                  EV += sgnf(h1 - h2) * raisedpot / 2.0; // check, bet, raise, call
                                } else {
                                  EV += 1.0 + bet; // check, bet, raise, fold
                                }
                            }
                          else {
                            EV += sgnf(h1 - h2) * pot / 2.0;  // check, bet, call
                          }
                        } else {
                          EV += -1.0;  // check, bet, fold
                        }
                    } else {
                        EV += sgnf(h1 - h2);  // check, check
                    }
                }

            }  // iter

             EV /= I;

             std::cout << y1 << "   " << y2 << "  " << w1 << "  " << w2 << "  " << z1 << "  " << z2 << "   EV =  " << EV << "\n";

             // min EV is best for player 2
              if(EV < minEV) {
                  minEV = EV;
                  y1b = y1;
                  y2b = y2;
                  z1b = z1;
                  z2b = z2;
                  w1b = w1;
                  w2b = w2;
                }

                }  // z2
            }     // z1

              // maximizing min EV is best for player 1
              if(minEV > maxminEV) {
                  maxminEV = minEV;
                  y1bb = y1b;
                  y2bb = y2b;
                  z1bb = z1b;
                  z2bb = z2b;
                  w1bb = w1b;
                  w2bb = w2b;
                }

          } // w2
        } // w1
      }  // y2
    }     // y1

  std::cout << "MIN-MAX:  " << y1bb << "  " << y2bb << "  " << w1bb << "  " << w2bb << "  " << z1bb << "  " << z2bb << "   EV =  " << maxminEV << "\n";

}

void Game::playD() {

  double s = 3;
  double pot = 2;

  double y1 = 0.0;
  double y2 = 0.8;

  std::vector<double> cdf;
  double ptotal = y1 + (1 - y2);
  double h = 0;
  double hstep = 0.01;
  double odds = s / (pot + 2 * s);

  while(h <= 1.001) {
      if(h <= y1)
        cdf.push_back(h / ptotal);
      else if(h < y2) {
          cdf.push_back(y1 / ptotal);
        } else {
          cdf.push_back((y1 + h - y2) / ptotal);
        }
      h += hstep;
    }

   int N = cdf.size();

   double q = 0;
   for(int k = 0; k < N; k++) {
       std::cout << q << "  " << cdf[k] << "\n";
       q += hstep;
     }

   double z = 0;
   int i = 0;
   while(cdf[i] < odds) {
       ++i;
       z += hstep;
     }
   int I = i;
   std::cout << " I = " << I << "\n";

   std::mt19937 genh;
   genh.seed(53);
   std::uniform_real_distribution<double> unifh(0, 1);


   for(double x = 0; x <= 1; x += hstep) {

     double h = 0;
     i = 0;
     while(h < x) {
         h += hstep;
         ++i;
       }
     double evcall = 0.5 * pot * x;  // ev of call of player 1
     while(h <= 1) {
         evcall += (1 - 2*cdf[i]) * (s + 1) * hstep;
         h += hstep;
         ++i;
       }
     double evcheck = y1 - (1-y2); // ev of check

     double ev = (y1 + (1-y2)) * evcall + (1 - (y1 + (1-y2))) * evcheck;

     double ev1 = 0;
     double T = 5000;
     for(int k = 0; k < T; k++) {

         double h1 = unifh(genh);
         double h2 = unifh(genh);
         if(h1 <= y1 || h1 >= y2) {
             if(h2 >= x) {
                 ev1 += sgnf(h1 - h2) * (pot + 2 * s) / 2.0; // bet, call
               } else {
                 ev1 += pot / 2; // bet, fold
               }
           } else {
             ev1 += sgnf(h1 - h2) * pot / 2.0; // check
           }

       }

     ev1 /= T;

     std::cout << x << "  " << ev << "  " << ev1 << "\n";


     } // x

}


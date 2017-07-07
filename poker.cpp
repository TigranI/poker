#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "poker.hpp"


Game::Game(uint64_t n, uint64_t k, uint64_t s, double alpha) : N(n), K(k), S(s), p1b(N, 0), p2c(N, 0), ps(S+1,0), alpha(alpha) { }

void Game::setup() {

  p1b[0] = 1.0;
  p2c[0] = 1.0;
  p2c[N-1] = 0.0;

  for(int i = 0; i < S+1; ++i) {
      ps[i] = static_cast<double>(i) / S;
  }

}

void Game::play() {

//  std::vector<std::vector<double> > ev(S, std::vector<double>(S,0)); // 4-vector

//  for(int i = 0; i < N; ++i) {
//    for(int j = 0; j < N; ++j) {
//          if(i == j)
//            continue;

//      for(int s = 0; s < S + 1; ++s) {
//          for(int t = 0; t < S + 1; ++t) {
//              ev[s][t] += sgnf(j - i) * (alpha + 1) * ps[s] * ps[t] + sgnf(j - i) + (1 - ps[s]);
//          }
//      }

//    }
//  }

  std::mt19937 gen;
  gen.seed(17);

  std::uniform_real_distribution<double> unifh(0, 1);
 // std::uniform_real_distribution<double> unifp(0, 1);

  double alpha = 1;
  uint64_t T = 10000;
  double evmaxmin = -std::numeric_limits<double>::infinity();
  double cmaxmin = 0;
  double bmaxmin = 0;
  double cmin = 0;

  for(double b = 0;b <= 1.01; b += 0.5) {
      double evmin = std::numeric_limits<double>::infinity();

      for(double c = 0; c <= 1.01; c += 0.1) {
           double ev = 0;
           for(int i = 0; i < T; ++i) {
             double h1 = unifh(gen);
             double h2 = unifh(gen);
             double pb = (1-b)*h1 + b;
             double pc = (1-c)*h2 + c;
             ev += sgnf(h1 - h2) * (alpha + 1) * pb * pc + sgnf(h1 - h2) + pb * (1 - pc); // bet,call + nobet + bet,fold
           }  // i
           ev /= T;
           std::cout  << b << "  " << ev << "  " << evmin << "  " << cmaxmin << "\n";
             if(ev < evmin) {
                evmin = ev;
                cmin = c;
             }

      } // c

      if(evmin > evmaxmin) {
          evmaxmin = evmin;
          bmaxmin = b;
          cmaxmin = cmin;
      }
   } // b

  std::cout  << evmaxmin << "    " << bmaxmin << "  " << cmaxmin << "\n";

}

class Path {
public:
  double ev;
  int opp;
  std::vector<int> p;
  std::vector<int> q;

  Path() = default;

  Path(double ev, int opp,  std::vector<int> p, std::vector<int> q) : ev(ev), opp(opp), p(p), q(q) {}
};

class Tree {
public:
  double ev;
  int opp;
  std::vector<Path> paths;

  Tree() = default;

  Tree(int opp) : opp(opp) {
    paths.resize(10);
  }

  void addPath(Path path, int pos) {
    paths[pos] = path;
  }

  void setEV(double ev1) {
    ev = ev1;
  }

};

// rock paper scissors learning
void Game::rps() {

  int N = 3; // basis size
  int K = 3;
  double dp = 0.05; // step
  int T = 30;
  int pathtotal = 0;

  std::vector<double> P(K);  // 3 probabilities
  std::vector<std::vector<double> > Q(N, std::vector<double>(K));

  std::vector<Tree> trees(N); // N basic opponents

  // ev = (p0q1 + p1q2 + p2q0) - (p0q2 + p1q0 + p2q1)

  P[0] = 0.7; P[1] = 0.2; P[2] = 1 - (P[0] + P[1]);
  Q[0][0] = 0.8; Q[0][1] = 0.1; Q[0][2] = 1 - (Q[0][0] + Q[0][1]);
  Q[1][0] = 0.1; Q[1][1] = 0.8; Q[1][2] = 1 - (Q[1][0] + Q[1][1]);
  Q[2][0] = 0.1; Q[2][1] = 0.1; Q[2][2] = 1 - (Q[2][0] + Q[2][1]);

  // initial setup

  for(int i = 0; i < N; ++i) {
      double evTotal = 0;
      int pos = 0;
      int sign = 1;
      Tree tree(i);
      for(int j = 0; j < K; ++j) {
          for(int k = 0; k < K; ++k) {
              if(j == k) continue;  // ev = 0 case
              std::vector<int> p1(1, j);
              std::vector<int> q1(1, k);
              double ev = P[j] * Q[i][k];
              if(j - k == 1 || k - j == 2) {
                 sign = -1.0;
                 ev *= sign;
              } else
                sign = 1.0;
              evTotal += ev;
             // std::cout << i << "  " << j << "  " << k << "  " << sign << "\n";
              Path path(ev, i, p1, q1);
              tree.addPath(path, pos++);
          }
      }
      tree.setEV(evTotal);
      trees[i] = tree;
      pathtotal = pos;
  }

   std::cout << trees[0].ev << "  "<< trees[1].ev << "  "<< trees[2].ev << "\n";

  // learning phase


  for(int iter = 0; iter < T; ++iter) {

      double evmin = std::numeric_limits<double>::infinity();
      int minpos = 0;

      for(int i = 0; i < N; ++i) {
          if(trees[i].ev < evmin) {
              evmin = trees[i].ev;
              minpos = i;
          }
      }

      if(evmin >= 0) {
        std::cout << "all opponents beaten" << "\n";
        break;
      }

      std::cout << P[0] << "  "<< P[1] << "  "<< P[2] << "  " << evmin << "\n";

      double evpathmin = std::numeric_limits<double>::infinity();
      int minpathpos = 0;
      for(int j=0; j < pathtotal; ++j) {
          if(trees[minpos].paths[j].ev < evpathmin) {
              evpathmin = trees[minpos].paths[j].ev;
              minpathpos = j;
          }
      }

      int index = trees[minpos].paths[minpathpos].p[0];
      P[index] -= dp;
      double totalP = P[0] + P[1] + P[2];  // normalize prob
      for(int s=0; s < K; ++s) {
          P[s] /= totalP;
      }

      // play again with new P

      for(int i = 0; i < N; ++i) {
          double evTotal = 0;
          int sign = 1;
          for(int j = 0; j < K; ++j) {
              for(int k = 0; k < K; ++k) {
                  if(j == k) continue;  // ev = 0 case

                  double ev = P[j] * Q[i][k];
                  if(j - k == 1 || k - j == 2) {
                     sign = -1.0;
                     ev *= sign;
                  } else
                    sign = 1.0;
                  evTotal += ev;
              }
          }
          trees[i].setEV(evTotal);
      }
  }

}

void Game::onestreet() {

  // define B

  // define C(B)

  // play

  int I = 3;
  double maxB = 2; // 2*pot, size > 2*pot should be occasional, pot = 1
  double minBper = 0.025;
  double eps = 0.001;
  std::mt19937 gen;
  gen.seed(17);
  std::uniform_real_distribution<double> unifh(0, 1);
  std::uniform_real_distribution<double> unifs(0, 1);

  // param list (c1, c2, a)
  for(double c1 = 0.5; c1 <= 0.5; c1 += 0.2) {  // mean bet = step(c1)
 // for(double c1 = 1; c1 <= 4; c1 += 1.0) {  // mean bet = c1*x + c2
      double c2 = 1 - c1;
     // for(double c2 = 1; c2 <= -1; c2 += 1) {
          for(double a = 0.5; a <= 2; a += 2.9) {   // cdf = x^a, pdf = a * x^(a-1)

              double ev = 0;
              double avebet = 0;
              for(int iter = 0; iter < I; ++iter) {

                  double h1 = unifh(gen);
                  double h2 = unifh(gen);

                  // how much to bet
                 // double b0 = linearFn(c1, c2, h1);
                  double b0 = stepFn(c1, h1);
                  if(b0 < minBper) {
                      ev += sgnf(h1 - h2);  // no bet
                     // std::cout << "no bet: " << h1 << "  " << h2 << "  " << ev << "\n";
                      continue;
                    }
                  b0 = std::min(b0, 1 - minBper);
                  double y = unifs(gen);  // uniform for power cdf
                  double scale = std::min(1-b0, b0);
                  double db = std::pow(y, 1.0/a) * scale;
                  double s = sgnf(unifs(gen) - 0.5);  // symmetry around b0: b0 +- db
                  double betsize = b0 + s * db;  // % of maxB
                  if(betsize > 1)
                     betsize = 1;
                  else if(betsize < 0)
                     betsize = 0;

                  avebet += betsize;

                  // call or not?
                  double potsize = 1 + 2 * maxB * betsize; // potsize if called, original pot + bet + call
                  double odds = maxB * betsize / potsize;
                  bool isCall = callOdds(h2, betsize, c1, c2, a, odds, minBper);

                  std::string res = "wrong";
                   if((h1 >= h2 && !isCall) || (h1 < h2 && isCall))
                     res = "correct";
                // std::cout <<  " betinfo =  " << b0 << "   " << db << "   " << betsize << "\n";
                //  std::cout << h1 << "  " << h2 << "  " << betsize << "    " << odds << "     " << isCall << "\n";
                //std::cout << h1 << "  " << h2 << "    " << isCall << "   " << res << "\n";

                  if(isCall) {
                      ev += sgnf(h1 - h2) * potsize;  // bet, call
                     // std::cout << "call: " << h1 << "  " << h2 << "  " << ev << "\n";
                    }
                  else {
                      ev += 1;   // bet, fold
                    //  std::cout << "fold: " << h1 << "  " << h2 << "  " << ev << "\n";
                    }



                } // iter

              avebet /= I;
              ev /= I;
              std::cout << c1 << "  " << c2 << "  " << a << "  " << avebet << "  EV =  " << ev << "\n";
            }  // a
      //  } // c2
    } // c1

}

// computes distr for a fixed bet. e --> b0, betsize --> db = (betsize-b0) --> pdf(db)
// to call the cdf of computed prob distribution at point h needs to be >= odds
// which means that most hands with given betsize are "weaker" than h
bool Game::callOdds(double h, double b, double c1, double c2, double a, double odds, double minBper) {

  double eps = 0.01;
  double ptotal = 0;
  double estep = 0.05;
  std::vector<double> pb;
  for(double x = estep; x <= 1 - estep + 0.001; x += estep) {   // equity

      // double b0 = linearFn(c1, c2, x);
      double b0 = stepFn(c1, x);
      double prob = 0;
      if(b0 < minBper) {
        pb.push_back(prob);
        continue;
        }
      b0 = std::min(b0, 1 - minBper);
      double rawdb = std::fabs(b - b0);
      double scale = std::min(1 - b0, b0);
      double db = rawdb / scale;
      if(db < 1)
        prob = a * std::pow(db + eps, a - 1) / 2;  // pdf(db) / 2 since bet = b +- db
      ptotal += prob;
      pb.push_back(prob);
     // std::cout << x << "  " << db << "  " << prob << "\n";
    }

  for(int i=0; i < pb.size(); ++i) {
      pb[i] /= ptotal;
      std::cout << pb[i] << "\n";
    }

  std::cout << "\n \n";

//  double sum = 0;
//  for(double p : pb) {
//      sum += p;
//    }
//  std::cout << sum << "\n";


   double hcdf = 0;
   int i = 0;
   double x = estep;
   while(x <= h) {
       hcdf += pb[i++];
       x += estep;
     }

// std::cout << " in callOdds() : "  <<  hcdf << "    " << odds << "  " << h << "\n";


   if(hcdf >= odds)
     return true;
   else
     return false;
}

double Game::linearFn(double c1, double c2, double x) {

  return c1 * x + c2;
}

double Game::stepFn(double c, double x) {

  if(x < c)
    return 0.0;
  else
    return 1.0;
}

void Game::akq() {

  double a = 1;

  double pb = 0.4;
  double pc = 1;

  // ev = ace,bet,call + ace,bet,fold - queen,bet,call + queen,bet,fold - queen,check,check

  for(double pc = 0; pc <= 1 + 0.01; pc += 0.1) {

      double ev = 0.5 * 1 * pc * (a + 1) + 0.5 * 1 * (1 - pc) * 1 - 0.5 * pb * pc * (a + 1) + 0.5 * pb * (1 - pc) * 1 - 0.5 * (1 - pb) * 1;

      std::cout << pb << "  " << pc << "   ev = " << ev << "\n";
    }

}

bool Game::betStrategyA(double x1, double x2, double h1, double u) {

  double p = 0; // prob of bet
  if(h1 > x2) {
    p = 1.0;
    return true;
  }
  else if(h1 > x1) {
      p = 1.0 / (x2 - x1) * h1 + x1 / (x1 - x2);
      if(u < p)
       return true;
      else
        return false;
    }

  return 0;

}

bool Game::callStrategyA(double x1, double x2, double h, double odds) {

  double ptotal = 0;
  double estep = 0.025;
  std::vector<double> pb;

  for(double x = 0; x <= 1 + 0.001; x += estep) {   // equity

      double p = 0; // prob of bet
      if(x > x2) {
        p = 1.0;
      }
      else if(x > x1) {
          p = 1.0 / (x2 - x1) * x + x1 / (x1 - x2);
        }
      ptotal += p;
      pb.push_back(p);
     // std::cout << x << "  " << db << "  " << prob << "\n";
    }

  for(int i=0; i < pb.size(); ++i) {
      pb[i] /= ptotal;
     // std::cout << pb[i] << "\n";
    }

 // std::cout << "\n \n";

//  double sum = 0;
//  for(double p : pb) {
//      sum += p;
//    }
//  std::cout << sum << "\n";


   double hcdf = 0;
   int i = 0;
   double x = 0;
   while(x <= h) {
       hcdf += pb[i++];
       x += estep;
     }

// std::cout << " in callStrategyA() : "  <<  hcdf << "    " << odds << "  " << h << "\n";


   if(hcdf >= odds)
     return true;
   else
     return false;

}

void Game::playA() {

  std::mt19937 gen;
  gen.seed(17);
  std::uniform_real_distribution<double> unifh(0, 1);
  std::uniform_real_distribution<double> unifp(0, 1);
  double xstep = 0.1;
  int I = 20000;
  double bet = 1;
  double pot = 1 + 2 * bet;

  // param list (x1, x2)
  for(double x1 = 0; x1 <= 1 - xstep; x1 += xstep) {  // lower bound
    for(double x2 = x1 + xstep; x2 <= 1; x2 += xstep) {  // upper bound
        double ev = 0;
        for(int iter = 0; iter < I; ++iter) {

              double h1 = unifh(gen);
              double h2 = unifh(gen);
              double u = unifp(gen);
              if(betStrategyA(x1, x2, h1, u)) {
                  double odds = bet / pot;
                  if(callStrategyA(x1, x2, h2, odds)) {
                      ev += sgnf(h1 - h2) * pot;  // bet, call
                    } else {
                      ev += 1.0;   // bet, fold
                    }
                } else {
                  ev += sgnf(h1 - h2);  // no bet
                }
          } // I
         ev /= I;
         std::cout << x1 << "  " << x2 << "    EV =  " << ev << "\n";
      }  // x2
    }  // x1


}

bool Game::betStrategyB(double x1, double x2, double p1, double p2, double h1, double u) {

  double p = 0; // prob of bet
  if(h1 > x2) {
    p = 1.0;
    return true;
  }
  else if(h1 > x1)
      p = p2;
  else
      p = p1;

  if(u < p)
       return true;
      else
        return false;


  return 0;

}

bool Game::callStrategyB(double x1, double x2, double p1, double p2, double h, double odds) {

  double ptotal = 0;
  double estep = 0.025;
  std::vector<double> pb;

  for(double x = 0; x <= 1 + 0.001; x += estep) {   // equity

      double p = p1; // prob of bet
      if(x > x2) {
        p = 1.0;
      }
      else if(x > x1) {
          p = p2;
        }
      ptotal += p;
      pb.push_back(p);
     // std::cout << x << "  " << db << "  " << prob << "\n";
    }

  for(int i=0; i < pb.size(); ++i) {
      pb[i] /= ptotal;
     // std::cout << pb[i] << "\n";
    }

 // std::cout << "\n \n";

//  double sum = 0;
//  for(double p : pb) {
//      sum += p;
//    }
//  std::cout << sum << "\n";


   double hcdf = 0;
   int i = 0;
   double x = 0;
   while(x <= h) {
       hcdf += pb[i++];
       x += estep;
     }

//std::cout << " in callStrategyB() : "  <<  hcdf << "    " << odds << "  " << h << "\n";


   if(hcdf >= odds)
     return true;
   else
     return false;

}

void Game::playB() {

  std::mt19937 gen;
  gen.seed(17);
  std::uniform_real_distribution<double> unifh(0, 1);
  std::uniform_real_distribution<double> unifp(0, 1);
  double xstep = 0.25;
  double pstep = 0.33;
  int I = 50000;
  double bet = 1;
  double pot = 2 + 2 * bet;

  double minEV = std::numeric_limits<double>::infinity();
  double maxEV = -std::numeric_limits<double>::infinity();

  // param list (x1, x2, p1, p2)  : p = step function
  for(double x1 = 0.25; x1 <= 0.25; x1 += xstep) {  // lower bound
      for(double x2 = 0.75; x2 <= 0.75; x2 += xstep) {  // upper bound
          for(double p1 = 0.33; p1 <= 1; p1 += pstep) {  // first step value
              for(double p2 = 0; p2 <= 0; p2 += pstep) {  // second step value
                  double ev = 0;
                  for(int iter = 0; iter < I; ++iter) {

                      double h1 = unifh(gen);
                      double h2 = unifh(gen);
                      double u = unifp(gen);
                      if(betStrategyB(x1, x2, p1, p2, h1, u)) {
                         // std::cout << h1 << "  " << h2 << "\n";
                          double odds = bet / pot;
                          if(callStrategyB(x1, x2, p1, p2, h2, odds)) {
                              ev += sgnf(h1 - h2) * pot / 2.0;  // bet, call
                            } else {
                              ev += 1.0;   // bet, fold
                            }
                        } else {
                          ev += sgnf(h1 - h2);  // no bet
                        }
                    } // I
                  ev /= I;
                  std::cout << x1 << "  " << x2 << "   " << p1 << "  " << p2 << "    EV =  " << ev << "\n";
//                  if(ev < minEV) {
//                      std::cout << x1 << "  " << x2 << "   " << p1 << "  " << p2 << "    EV =  " << ev << "\n";
//                      minEV = ev;
//                  }

//                  if(ev > maxEV) {
//                      std::cout << x1 << "  " << x2 << "   " << p1 << "  " << p2 << "    EV =  " << ev << "\n";
//                      maxEV = ev;
//                    }
                } // p1
            } // p2
        }  // x2
    }  // x1


}



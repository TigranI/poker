#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <tuple>
#include <vector>

#include "poker.hpp"

void Player::makeIntervals(std::vector<std::tuple<double, uint64_t, uint64_t> > x, std::vector<std::tuple<double, double, uint64_t> >& intervals) {

//  for(auto t : intervals) {
//      std::cout << "before: " << std::get<0>(t) << "   " << std::get<1>(t) << "   " << std::get<2>(t) << "\n";
//    }

  for(auto t : x) {
      int i = 0;
      while(i < intervals.size()) {
          if(std::get<0>(t) >= std::get<0>(intervals[i]) && std::get<0>(t) <= std::get<1>(intervals[i])) {
              break;
            }
              ++i;
           }
          if(i < intervals.size()) {
                auto t1 = std::make_tuple(std::get<0>(intervals[i]), std::get<0>(t), std::get<1>(t));
                auto t2 = std::make_tuple(std::get<0>(t), std::get<1>(intervals[i]), std::get<2>(t));
                intervals.erase(intervals.begin() + i);
                intervals.push_back(t1);
                intervals.push_back(t2);
            }
      }

//    std::cout << "\n";
//    for(auto t : intervals) {
//        std::cout << "after: " << std::get<0>(t) << "   " << std::get<1>(t) << "   " << std::get<2>(t) << "\n";
//      }
}


double OneHand::play() {

  std::vector<std::shared_ptr<State> > treePointers(N);
  for(int i = 0; i < N; ++i) {
      treePointers[i] = std::shared_ptr<State>(players[i].root);
    }

  s = 1;
  pot = 2; // 1 bet for each player
  bool done = false;
  std::tuple<uint64_t, uint64_t> lastmove = std::make_tuple(0, 0);

  //std::cout << "game on \n";

  // one iteration per player's move, one treePointer is updated per iteration
  while(!done) {
      uint64_t i = 0;
      while(treePointers[turnToMove-1]->children[i]->lastmove != lastmove)
        ++i;
      treePointers[turnToMove-1] = treePointers[turnToMove-1]->children[i];
      uint64_t j = 0;
      while(h[turnToMove-1] < std::get<0>(treePointers[turnToMove-1]->intervals[j]) || h[turnToMove-1] > std::get<1>(treePointers[turnToMove-1]->intervals[j]))
        ++j;
      lastmove = std::make_tuple(turnToMove, std::get<2>(treePointers[turnToMove-1]->intervals[j]));

      // handle bet sizing and rules for ending the hand
      if(std::get<1>(lastmove) == 2 || std::get<1>(lastmove) == 3) // bet or call
        pot += s;
      if(std::get<1>(lastmove) == 3 || std::get<1>(lastmove) == 5) // call or fold
        done = true;
      if(std::get<1>(lastmove) == 4) // for this game if first player to move checks, it's over = no more decisions
        done = true;

      if(done) {
          if(std::get<1>(lastmove) == 5 && turnToMove == 1) // if p1 folds
            result = - (pot - s) / 2.0;
          if(std::get<1>(lastmove) == 5 && turnToMove == 2) // if p2 folds
            result = (pot - s) / 2.0;
          if(std::get<1>(lastmove) != 5)      // showdown
              result = sgnf(h[1] - h[0]) * pot / 2.0;
        }

    //  std::cout << done << "  " << turnToMove << "  " << h[0] << "  " << h[1] << "  " << result << "\n";

       turnToMove = (turnToMove == 1) ? 2 : 1;
    }

  return result;
}

void FullGame::play() {

  const uint64_t I = 1000;
  uint64_t firstToMove = 1;
  std::mt19937 genh;
  genh.seed(19);
  std::uniform_real_distribution<double> unifh(0, 1);

  for(uint64_t iter = 0; iter < I; ++iter) {
      std::vector<double> h = {unifh(genh), unifh(genh)};
      OneHand oneHand(players, firstToMove, h);
      result += oneHand.play();
      firstToMove = (firstToMove == 1) ? 2 : 1;
    }

  std::cout << "result = " << result <<"\n";
}

void Learning::process() {

  std::vector<Player> players;
  // create players here
  std::vector<double> p1 = {0.3, 0.9, 0.6};
  Player player1(1, p1); // optimal player
  players.push_back(player1);

  double step = 0.1;
  for(double x1 = 0.1; x1 <= 0.4; x1 += step) {
   for(double x2 = x1 + step; x2 <= 0.7; x2 += step) {
       for(double y = 0.4; y <= 0.8; y += step) {
        std::vector<double> p2 = {x1, x2, y};
        Player player2(2, p2);
        players.push_back(player2);

      //  player1.test();
      //  players[1].test();

        FullGame fullgame(players);
        fullgame.play(); // many hands played
        players.erase(players.begin() + 1);
    }  // y
   } // x2
  }  // x1

}



//************ test code below ************

// Player player0(p);
// players.push_back(player0);

//  std::cout << "poker 1: " << (player0.root->children)[0]->x[0] << "\n";

//  std::cout << "poker 2: " << (players[0].root->children)[0]->x[0] << "\n";

// player0.~Player();

//  std::vector<std::tuple<double, uint64_t, uint64_t> > x;
//  std::vector<std::tuple<double, double, uint64_t> > intervals;
//  auto i1 = std::make_tuple(0, 0.5, 2);
//  auto i2 = std::make_tuple(0.6, 0.8, 3);
//  auto i3 = std::make_tuple(0.85, 1, 3);
//  intervals.push_back(i1);
//  intervals.push_back(i2);
//  intervals.push_back(i3);
//  auto x1 = std::make_tuple(0.24, 1, 2);
//  auto x2 = std::make_tuple(0.75, 4, 3);
//  auto x3 = std::make_tuple(0.85, 4, 3);
//  x.push_back(x1);
//  x.push_back(x2);
//  x.push_back(x3);

//   void destroy_tree(State* State)
//   {
//     if( State != nullptr )  // it's not
//     {
//         for(int i = 0;i < State->children.size(); ++i)
//           destroy_tree(State->children[i]);
//         delete State;
//     }
//   }

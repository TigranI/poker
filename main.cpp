#include <cmath>
#include <fstream>
#include <iostream>

#include "poker.hpp"

int main()
{
 // std::cout << "Hello, Poker!" << endl;

  // Game game(4, 2, 4, 1.0);
 // Game game;
 // game.setup();
 // game.play();
 // game.rps();
 // game.onestreet();
// game.akq();
 // game.playA();
 // game.playB();
 // game.playC();
 // game.playT();
 // game.playD();
 // game.playAA();

 // Learning tournament;
 // tournament.process();

  uint64_t a = 16 + 32 + std::pow(2, 12) + std::pow(2, 13); // 2-pair, 3-pair same answer - bad
  uint64_t b = a % 15;
  std::cout << b << std::endl;

  uint64_t M = std::pow(10, 1);
  uint64_t w = 0;
  std::ofstream myfile;
  myfile.open ("example.txt");

  for(uint64_t i = 0; i < M; ++i) {
      if(w % 7 == 0)
        w += 1;
      else
        w += 4;
      myfile << w << "\n";
    }
 myfile.close();

  return 0;
}


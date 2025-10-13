#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>

#include "rlbwt.hpp"

int main() {
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT: TCGAT$ATA
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

    std::cout << "Size: " << move_lf.size() << std::endl;
    std::cout << "Move Runs: " << move_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_lf.permutation_runs() << std::endl;
    using Position = typename MoveLF<>::Position;
    Position pos = move_lf.get_position();
    for (size_t i = 0; i <= move_lf.size(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Run Data: " << move_lf.get_character() << std::endl;
        move_lf.next();
        pos = move_lf.get_position();
    }
    return 0;
}
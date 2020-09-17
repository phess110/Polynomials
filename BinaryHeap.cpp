#include "BinaryHeap.h"
#include <iostream>


int main() {
    //auto t1 = std::chrono::high_resolution_clock::now();
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    //std::cout << duration << std::endl;
    //std::cout << std::this_thread::get_id() << std::endl;

    //std::vector<PtValPair> v = { {1,2}, {2,3}, {5,7} };
    //PolyValGenerator p = PolyValGenerator(v);

    //std::cout << p.Eval(0.5) << std::endl;
    //std::cout << p.Eval(3) << std::endl;
    //std::cout << p.Eval(4) << std::endl;
    //std::cout << p.Eval(-10) << std::endl;

    //Polynomial::PolyPair r = Polynomial::PolyDiv(p, q);

    BinaryHeap<uint32_t, std::string, max_heap_comp<uint32_t>> heap;

    return 0;
}
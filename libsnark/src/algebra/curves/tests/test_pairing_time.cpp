#include "algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "algebra/curves/bls48/bls48_pp.hpp"
#include "iostream"
#include <chrono>

using namespace libsnark;



template<typename ppT>
void pairing_timing_test()
{
    int repeat_count = 100;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    double Timing_aP = 0;
    double Timing_bQ = 0;
    double Timing_pairing = 0;
    printf("Running timing tests:\n");


    for(int i=0; i<repeat_count; i++){
      //GT<ppT> GT_one = GT<ppT>::one();

      start = std::chrono::high_resolution_clock::now();
      G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
      end = std::chrono::high_resolution_clock::now();
      Timing_aP = Timing_aP +  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

      start = std::chrono::high_resolution_clock::now();
      G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();
      end = std::chrono::high_resolution_clock::now();
      Timing_bQ = Timing_bQ +  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      //std::cout << "Q is in G2"<< Q.is_well_formed() << std::endl;
      //assert(Q.Y.squared() - Q.X.squared()*Q.X == (Q.Z.squared()*Q.Z).squared());

      //G2<ppT> Q = Fr<ppT>("3") * G2<ppT>::one();

      //printf("GT_one:\n");
      //GT_one.print();

  /*
      printf("P:\n");
      P.print();
      P.print_coordinates();
      printf("Q:\n");
      Q.print();
      Q.print_coordinates();
      printf("\n\n");
  */

      //printf("e(sP, Q):\n");
      //GT<ppT> pre_ans1 = ppT::pairing(sP, Q);
      //pre_ans1.print();
      start = std::chrono::high_resolution_clock::now();
      ppT::reduced_pairing(P, Q);
      end = std::chrono::high_resolution_clock::now();
      Timing_pairing = Timing_pairing + std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();

      //std::cout << "Fr<ppT>::field_char()" << Fr<ppT>::field_char() << std::endl;
    }

    std::cout << "Average time of a*P = " << Timing_aP/repeat_count << " [ns]" << std::endl;
    std::cout << "Average time of b*Q = " << Timing_bQ/repeat_count << "[ns]" << std::endl;
    std::cout << "Average time of pairing " << Timing_pairing/repeat_count << "[ns]" << std::endl;
}


int main(void)
{


    std::cout << "alt_bn128 test start" << std::endl;
    alt_bn128_pp::init_public_params();
    pairing_timing_test<alt_bn128_pp>();



    std::cout << "bls test start" << std::endl;
    bls48_pp::init_public_params();
    pairing_timing_test<bls48_pp>();



}

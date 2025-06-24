//
// Created by 2646jiaxing on 2023/12/7.
//

#ifndef RHOMBUS_PARAMETERS_H
#define RHOMBUS_PARAMETERS_H

#include <cstdint>

namespace antchain::global{
    // GLOBAL Parameters.

    // logN: 12, 13
    constexpr uint32_t kLogPolyDegree = 13;

    // N: 4096, 8192
    constexpr uint32_t kPolyDegree = 8192;

    // total number of modulus, the last one is the special modulus used in key-switching
    constexpr uint32_t kNumModulus = 3; // only support 2 ~ 4 now

    // kSSBitLen means the plaintext modulus is 2^k, only support 2-power modulus now
    constexpr uint32_t kSSBitLen = 37; // support 8 ~ 64

    // +----------------------------------------------------------------------------+
    // Choose the modulus parameters

    // Case: N = 4096, q = {54, 55}
    // If the matrix and the vector are both 10 bits, the product is smaller than 25 bits
    // Then you can choose this modulus.
//     constexpr uint64_t kModulus[kNumModulus]{
//          36028797018652673, 18014398508400641};

    // Case N = 4096, q = {46, 46, 17}. If rotation is not required in your case, you can choose this parameter.
    // In fact, we don't need the third modulus which is only for the compatibility with SEAL decryption.
    //    constexpr uint64_t kModulus[kNumModulus]{
    //            70368743587841, 70368743669761, 114689};

    // Case N = 4096, q = {36, 36, 37}
//    constexpr uint64_t kModulus[kNumModulus]{
//            68719230977,
//            68719403009,
//            137438822401};

    // Case: N = 8192, q = {60, 60}
    // constexpr uint64_t kModulus[kNumModulus]{
    //     1152921504606748673, 1152921504606830593
    // };

    // Case: N = 8192, q = {50, 50, 60}.
    constexpr uint64_t kModulus[kNumModulus]{
            (1ULL << 50) - (1ULL << 14) + 1ULL,
            (1ULL << 50) - (1ULL << 38) + 1ULL,
            (1ULL << 59) | (1ULL << 47) | 1ULL};

    // Case: N = 8192, q = {50, 50, 50, 60}
    // constexpr uint64_t kModulus[kNumModulus]{
    //         (1ULL << 50) - (1ULL << 14) + 1ULL,
    //         (1ULL << 50) - (1ULL << 38) + 1ULL,
    //         (1ULL << 49) | (1ULL << 30) + 1ULL,
    //         (1ULL << 59) | (1ULL << 47) | 1ULL};

    // +----------------------------------------------------------------------------+
    // SPP: split-point
    // You can modify the split-points according to the actual test
    // If the elements of the matrix are smaller than 16 bits, you'd better use int16_t to represent it.
    // Similar, matrix elements < 32 bits --> int32_t, matrix elements < 64 bits --> int64_t,

    // N = 4096, matrix<int16_t> (means that element type of the matrix is int16_t)
    // constexpr int kSPPMap[kLogPolyDegree + 1]{
    //     0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5};

    // N = 4096, matrix<int32_t>
//     constexpr int kSPPMap[kLogPolyDegree + 1]{
//         0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};

    // N = 4096, matrix<int64_t>
    // constexpr int kSPPMap[kLogPolyDegree + 1]{
    //     0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4};

    // N = 8192, matrix<int16_t>
    // constexpr int kSPPMap[kLogPolyDegree + 1]{
    //     0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};

    // N = 8192, matrix<int32_t>
    constexpr int kSPPMap[kLogPolyDegree + 1]{
            0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4};

    // N = 8192, matrix<int64_t>
    // constexpr int kSPPMap[kLogPolyDegree + 1]{
    //     0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4};
}

#endif //RHOMBUS_PARAMETERS_H

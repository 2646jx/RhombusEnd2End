/*
Authors: Deevashwer Rathee, Mayank Rathee
Modified by Zhichong Huang
Copyright:
Copyright (c) 2021 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "BuildingBlocks/aux-protocols.h"

#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/value-extension.h"

using namespace sci;

AuxProtocols::AuxProtocols(int party, sci::NetIO *io,
                           OTPack<sci::NetIO> *otpack) {
  this->party = party;
  this->io = io;
  this->otpack = otpack;
  this->mill = new MillionaireProtocol<sci::NetIO>(party, io, otpack);
  this->mill_and_eq =
      new MillionaireWithEquality<sci::NetIO>(party, io, otpack);
}

AuxProtocols::~AuxProtocols() {
  delete mill;
  delete mill_and_eq;
}

void AuxProtocols::wrap_computation(uint64_t *x, uint8_t *y, int32_t size,
                                    int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *tmp_x = new uint64_t[size];
  for (int i = 0; i < size; i++) {
    if (party == sci::ALICE)
      tmp_x[i] = x[i] & mask;
    else
      tmp_x[i] = (mask - x[i]) & mask;  // 2^{bw_x} - 1 - x[i]
  }
  std::cout << "start compare" << std::endl;
  mill->compare(y, tmp_x, size, bw_x, true);  // computing greater_than
  std::cout << "end compare" << std::endl;

  delete[] tmp_x;
}

void AuxProtocols::multiplexer(uint8_t *sel, uint64_t *x, uint64_t *y,
                               int32_t size, int32_t bw_x, int32_t bw_y) {
  assert(bw_x <= 64 && bw_y <= 64 && bw_y <= bw_x);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

  uint64_t *corr_data = new uint64_t[size];
  uint64_t *data_S = new uint64_t[size];
  uint64_t *data_R = new uint64_t[size];

  // y = (sel_0 \xor sel_1) * (x_0 + x_1)
  // y = (sel_0 + sel_1 - 2*sel_0*sel_1)*x_0 + (sel_0 + sel_1 -
  // 2*sel_0*sel_1)*x_1 y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)]
  //     + [sel_1*x_1 + sel_0*(x_1 - 2*sel_1*x_1)]
  for (int i = 0; i < size; i++) {
    corr_data[i] = (x[i] * (1 - 2 * uint64_t(sel[i]))) & mask_y;
  }
  if (party == sci::ALICE) {
    otpack->iknp_straight->send_cot(data_S, corr_data, size, bw_y);
    otpack->iknp_reversed->recv_cot(data_R, (bool *)sel, size, bw_y);
  } else {  // party == sci::BOB
    otpack->iknp_straight->recv_cot(data_R, (bool *)sel, size, bw_y);
    otpack->iknp_reversed->send_cot(data_S, corr_data, size, bw_y);
  }
  for (int i = 0; i < size; i++) {
    y[i] = ((x[i] * uint64_t(sel[i]) + data_R[i] - data_S[i]) & mask_y);
  }

  delete[] corr_data;
  delete[] data_S;
  delete[] data_R;
}

void AuxProtocols::B2A(uint8_t *x, uint64_t *y, int32_t size, int32_t bw_y) {
  assert(bw_y <= 64 && bw_y >= 1);
  if (bw_y == 1) {
    for (int i = 0; i < size; i++) {
      y[i] = uint64_t(x[i]) & 1;
    }
    return;
  }
  uint64_t mask = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

  if (party == sci::ALICE) {
    uint64_t *corr_data = new uint64_t[size];
    for (int i = 0; i < size; i++) {
      corr_data[i] = (-2 * uint64_t(x[i])) & mask;
    }
    otpack->iknp_straight->send_cot(y, corr_data, size, bw_y);

    for (int i = 0; i < size; i++) {
      y[i] = (uint64_t(x[i]) - y[i]) & mask;
    }
    delete[] corr_data;
  } else {  // party == sci::BOB
    otpack->iknp_straight->recv_cot(y, (bool *)x, size, bw_y);

    for (int i = 0; i < size; i++) {
      y[i] = (uint64_t(x[i]) + y[i]) & mask;
    }
  }
}

/*
 ALICE holds the table spec, BOB holds the input x
 table: size * 2^bw_x
*/
template <typename T>
void AuxProtocols::lookup_table(T **spec, T *x, T *y, int32_t size,
                                int32_t bw_x, int32_t bw_y) {
  if (party == sci::ALICE) {
    assert(x == nullptr);
    assert(y == nullptr);
  } else {  // party == sci::BOB
    assert(spec == nullptr);
  }
  assert(bw_x <= 8 && bw_x >= 2);
  int32_t T_size = sizeof(T) * 8;
  assert(bw_y <= T_size);

  T mask_x = (bw_x == T_size ? -1 : ((1ULL << bw_x) - 1));
  T mask_y = (bw_y == T_size ? -1 : ((1ULL << bw_y) - 1));
  uint64_t N = 1 << bw_x;

  if (party == sci::ALICE) {
    PRG128 prg;
    T **data = new T *[size];
    for (int i = 0; i < size; i++) {
      data[i] = new T[N];
      for (uint64_t j = 0; j < N; j++) {
        data[i][j] = spec[i][j];
      }
    }

    // (bw_x-1)-th OT is bw_x bits, that is, 1oo2^{bw_x} OT
    otpack->kkot[bw_x - 1]->send(data, size, bw_y);

    for (int i = 0; i < size; i++) delete[] data[i];
    delete[] data;
  } else {  // party == sci::BOB
    uint8_t *choice = new uint8_t[size];
    for (int i = 0; i < size; i++) {
      choice[i] = x[i] & mask_x;
    }
    otpack->kkot[bw_x - 1]->recv(y, choice, size, bw_y);

    delete[] choice;
  }
}

void AuxProtocols::MSB(uint64_t *x, uint8_t *msb_x, int32_t size,
                       int32_t bw_x) {
  assert(bw_x <= 64);
  int32_t shift = bw_x - 1;
  uint64_t shift_mask = (shift == 64 ? -1 : ((1ULL << shift) - 1));

  uint64_t *tmp_x = new uint64_t[size];
  uint8_t *msb_xb = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    tmp_x[i] = x[i] & shift_mask;
    msb_xb[i] = (x[i] >> shift) & 1;
    if (party == sci::BOB) tmp_x[i] = (shift_mask - tmp_x[i]) & shift_mask;
  }

  mill->compare(msb_x, tmp_x, size, bw_x - 1, true);  // computing greater_than

  for (int i = 0; i < size; i++) {
    msb_x[i] = msb_x[i] ^ msb_xb[i];
  }

  delete[] tmp_x;
  delete[] msb_xb;
}

void AuxProtocols::MSB_to_Wrap(uint64_t *x, uint8_t *msb_x, uint8_t *wrap_x,
                               int32_t size, int32_t bw_x) {
  assert(bw_x <= 64);
  if (party == sci::ALICE) {
    PRG128 prg;
    prg.random_bool((bool *)wrap_x, size);
    uint8_t **spec = new uint8_t *[size];
    for (int i = 0; i < size; i++) {
      spec[i] = new uint8_t[4];
      uint8_t msb_xb = (x[i] >> (bw_x - 1)) & 1;
      for (int j = 0; j < 4; j++) {
        uint8_t bits_j[2];  // j0 || j1 (LSB to MSB)
        uint8_to_bool(bits_j, j, 2);
        spec[i][j] = (((1 ^ msb_x[i] ^ bits_j[0]) * (msb_xb ^ bits_j[1])) ^
                      (msb_xb * bits_j[1]) ^ wrap_x[i]) &
                     1;
      }
    }
    lookup_table<uint8_t>(spec, nullptr, nullptr, size, 2, 1);

    for (int i = 0; i < size; i++) delete[] spec[i];
    delete[] spec;
  } else {  // party == sci::BOB
    uint8_t *lut_in = new uint8_t[size];
    for (int i = 0; i < size; i++) {
      lut_in[i] = (((x[i] >> (bw_x - 1)) & 1) << 1) | msb_x[i];
    }
    lookup_table<uint8_t>(nullptr, lut_in, wrap_x, size, 2, 1);

    delete[] lut_in;
  }
}

void AuxProtocols::msb0_to_wrap(uint64_t *x, uint8_t *wrap_x, int32_t size,
                                int32_t bw_x) {
  assert(bw_x <= 64);
  if (party == sci::ALICE) {
    PRG128 prg;
    prg.random_bool((bool *)wrap_x, size);
    uint8_t **spec = new uint8_t *[size];
    for (int i = 0; i < size; i++) {
      spec[i] = new uint8_t[2];
      uint8_t msb_xb = (x[i] >> (bw_x - 1)) & 1;
      spec[i][0] = wrap_x[i] ^ msb_xb;
      spec[i][1] = wrap_x[i] ^ 1;
    }
#if USE_CHEETAH
    otpack->silent_ot->send_ot_cm_cc<uint8_t>(spec, size, 1);
#else
    otpack->iknp_straight->send(spec, size, 1);
#endif

    for (int i = 0; i < size; i++) delete[] spec[i];
    delete[] spec;
  } else {  // party == sci::BOB
    uint8_t *msb_xb = new uint8_t[size];
    for (int i = 0; i < size; i++) {
      msb_xb[i] = (x[i] >> (bw_x - 1)) & 1;
    }
#if USE_CHEETAH
    otpack->silent_ot->recv_ot_cm_cc<uint8_t>(wrap_x, msb_xb, size, 1);
#else
    otpack->iknp_straight->recv(wrap_x, msb_xb, size, 1);
#endif

    delete[] msb_xb;
  }
}

void AuxProtocols::msb1_to_wrap(uint64_t *x, uint8_t *wrap_x, int32_t size,
                                int32_t bw_x) {
  assert(bw_x <= 64);
  if (party == sci::ALICE) {
    PRG128 prg;
    prg.random_bool((bool *)wrap_x, size);
    uint8_t **spec = new uint8_t *[size];
    for (int i = 0; i < size; i++) {
      spec[i] = new uint8_t[2];
      uint8_t msb_xb = (x[i] >> (bw_x - 1)) & 1;
      spec[i][0] = wrap_x[i];
      spec[i][1] = wrap_x[i] ^ msb_xb;
    }
#if USE_CHEETAH
    otpack->silent_ot->send_ot_cm_cc<uint8_t>(spec, size, 1);
#else

    otpack->iknp_straight->send(spec, size, 1);
#endif

    for (int i = 0; i < size; i++) delete[] spec[i];
    delete[] spec;
  } else {  // party == sci::BOB
    uint8_t *msb_xb = new uint8_t[size];
    for (int i = 0; i < size; i++) {
      msb_xb[i] = (x[i] >> (bw_x - 1)) & 1;
    }
#if USE_CHEETAH
    otpack->silent_ot->recv_ot_cm_cc<uint8_t>(wrap_x, msb_xb, size, 1);
#else
    otpack->iknp_straight->recv(wrap_x, msb_xb, size, 1);
#endif

    delete[] msb_xb;
  }
}

void AuxProtocols::AND(uint8_t *x, uint8_t *y, uint8_t *z, int32_t size) {
  int old_size = size;
  size = ceil(size / 8.0) * 8;
  uint8_t *tmp_x = new uint8_t[size];
  uint8_t *tmp_y = new uint8_t[size];
  uint8_t *tmp_z = new uint8_t[size];
  memcpy(tmp_x, x, old_size * sizeof(uint8_t));
  memcpy(tmp_y, y, old_size * sizeof(uint8_t));

  // assert((size % 8) == 0);
  Triple triples_std(size, true);
  this->mill->triple_gen->generate(party, &triples_std, _16KKOT_to_4OT);

  uint8_t *ei = new uint8_t[(size) / 8];
  uint8_t *fi = new uint8_t[(size) / 8];
  uint8_t *e = new uint8_t[(size) / 8];
  uint8_t *f = new uint8_t[(size) / 8];

  // this->mill->AND_step_1(ei, fi, x, y, triples_std.ai, triples_std.bi, size);
  this->mill->AND_step_1(ei, fi, tmp_x, tmp_y, triples_std.ai, triples_std.bi,
                         size);

  int size_used = size / 8;
  if (party == sci::ALICE) {
    // Send share of e and f
    io->send_data(ei, size_used);
    io->send_data(ei, size_used);
    io->send_data(fi, size_used);
    io->send_data(fi, size_used);
    // Receive share of e and f
    io->recv_data(e, size_used);
    io->recv_data(e, size_used);
    io->recv_data(f, size_used);
    io->recv_data(f, size_used);
  } else  // party = sci::sci::BOB
  {
    // Receive share of e and f
    io->recv_data(e, size_used);
    io->recv_data(e, size_used);
    io->recv_data(f, size_used);
    io->recv_data(f, size_used);
    // Send share of e and f
    io->send_data(ei, size_used);
    io->send_data(ei, size_used);
    io->send_data(fi, size_used);
    io->send_data(fi, size_used);
  }

  // Reconstruct e and f
  for (int i = 0; i < size_used; i++) {
    e[i] ^= ei[i];
    f[i] ^= fi[i];
  }

  // this->mill->AND_step_2(z, e, f, nullptr, nullptr, triples_std.ai,
  //         triples_std.bi, triples_std.ci, size);
  this->mill->AND_step_2(tmp_z, e, f, nullptr, nullptr, triples_std.ai,
                         triples_std.bi, triples_std.ci, size);
  memcpy(z, tmp_z, old_size * sizeof(uint8_t));

  // cleanup
  delete[] tmp_x;
  delete[] tmp_y;
  delete[] tmp_z;
  delete[] ei;
  delete[] fi;
  delete[] e;
  delete[] f;

  return;
}

void AuxProtocols::reduce(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                          int32_t bw_y) {
  assert(bw_y <= bw_x);
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

  for (int i = 0; i < dim; i++) {
    y[i] = x[i] & mask_y;
  }
}

void AuxProtocols::digit_decomposition(int32_t dim, uint64_t *x,
                                       uint64_t *x_digits, int32_t bw_x,
                                       int32_t digit_size) {
  assert(false && "Inefficient version of digit decomposition called");
  int num_digits = ceil(double(bw_x) / digit_size);
  int last_digit_size = bw_x - (num_digits - 1) * digit_size;
  uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
  uint64_t last_digit_mask =
      (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);

  Truncation trunc(this->party, this->io, this->otpack);
  for (int i = 0; i < num_digits; i++) {
    trunc.truncate_and_reduce(dim, x, x_digits + i * dim, i * digit_size, bw_x);
    uint64_t mask = (i == (num_digits - 1) ? last_digit_mask : digit_mask);
    for (int j = 0; j < dim; j++) {
      x_digits[i * dim + j] &= mask;
    }
  }
}

void AuxProtocols::digit_decomposition_sci(int32_t dim, uint64_t *x,
                                           uint64_t *x_digits, int32_t bw_x,
                                           int32_t digit_size,
                                           bool all_digit_size) {
  int num_digits = ceil(double(bw_x) / digit_size);
  int last_digit_size = bw_x - (num_digits - 1) * digit_size;
  uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
  uint64_t last_digit_mask =
      (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);
  for (int i = 0; i < num_digits; i++) {
    for (int j = 0; j < dim; j++) {
      // y_{b,0}, y_{b,1}, ..., y_{b,c-1}
      x_digits[i * dim + j] = (x[j] >> (i * digit_size));
      x_digits[i * dim + j] &=
          (i == (num_digits - 1)) ? last_digit_mask : digit_mask;
    }
  }
  uint8_t *wrap_ = new uint8_t[dim * (num_digits - 1)];
  uint8_t *ones_ = new uint8_t[dim * (num_digits - 1)];
  uint8_t *dp_wrap_entering = new uint8_t[dim * num_digits];
  uint8_t *dp_temp = new uint8_t[dim * num_digits];
  uint64_t *dp_wrap_arith = new uint64_t[dim * num_digits];
  // Fill wrap_ and ones_
  uint64_t *temp_x_digits = new uint64_t[dim * (num_digits - 1)];

  for (int i = 0; i < (num_digits - 1); i++) {
    for (int j = 0; j < dim; j++) {
      // y_{b,0}, y_{b,1}, ..., y_{b,c-1}
      if (party == sci::ALICE)
        temp_x_digits[i * dim + j] = x_digits[i * dim + j] & digit_mask; // The 'AND' seems repeated work
      // D=2^digit_size, D-1-y_{b,0}, D-1-y_{b,1}, ..., D-1-y_{b,c-1}
      else
        temp_x_digits[i * dim + j] =
            (digit_mask - x_digits[i * dim + j]) & digit_mask;
    }
  }
  
  // compute Wrap(y_{0,i}, y_{1,i}, D), Eq(y_{0,i} + y_{1,i} ?= D-1 mod D)
  this->mill_and_eq->compare_with_eq(wrap_, ones_, temp_x_digits,
                                     (dim * (num_digits - 1)), digit_size);

  // DP steps proceed
  for (int i = 0; i < num_digits; i++) {
    // Lemma 1 in paper SIRNN
    if (i > 0) {
      // step 2: (e AND c)
      this->AND(ones_ + (i - 1) * dim, dp_wrap_entering + (i - 1) * dim,
                dp_temp + (i - 1) * dim, dim);
    }
    for (int j = 0; j < dim; j++) {
      if (i == 0) {
        // step 1: c = 0
        dp_wrap_entering[i * dim + j] = 0;
      } else {
        // step 3: d xor (e AND c)
        dp_wrap_entering[i * dim + j] =
            wrap_[(i - 1) * dim + j] ^ dp_temp[(i - 1) * dim + j];
      }
    }
  }
  // carry
  this->B2A(dp_wrap_entering, dp_wrap_arith, num_digits * dim, digit_size);
  for (int i = 0; i < num_digits; i++) {
    for (int j = 0; j < dim; j++) {
      x_digits[i * dim + j] += dp_wrap_arith[i * dim + j];
      uint64_t temp_mask =
          (i == (num_digits - 1)) ? last_digit_mask : digit_mask;
      x_digits[i * dim + j] &= temp_mask;
    }
    if (all_digit_size) {
      if (i == (num_digits - 1)) {
        XTProtocol *xt = new XTProtocol(this->party, this->io, this->otpack);
        uint64_t *temp_last_digs = new uint64_t[dim];
        xt->z_extend(dim, x_digits + (num_digits - 1) * dim, temp_last_digs,
                     last_digit_size, digit_size);
        for (int j = 0; j < dim; j++) {
          x_digits[i * dim + j] = temp_last_digs[j];
          x_digits[i * dim + j] &= digit_mask;
        }
        delete xt;
        delete[] temp_last_digs;
      }
    }
  }

  delete[] wrap_;
  delete[] ones_;
  delete[] dp_wrap_entering;
  delete[] dp_temp;
  delete[] dp_wrap_arith;
  delete[] temp_x_digits;
}

uint64_t lookup_msnzb(uint64_t index) {
  uint64_t ret = 0ULL;
  ret = floor(log2(index));
  if (index == 0) {
    ret = 0ULL;
  }
  // In the above step only at max log(64) = 6 bits are filled.
  ret <<= 1;
  // Last bit stores 1 if index is 0, else 0.
  if (index == 0) {
    ret ^= 1ULL;
  }
  return ret;
}

void AuxProtocols::msnzb_sci(uint64_t *x, uint64_t *msnzb_index, int32_t bw_x,
                             int32_t size, int32_t digit_size) {
  // The protocol only works when digit_size divides bw_x.
  int32_t last_digit_size = bw_x % digit_size;
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t digit_mask = (digit_size == 64 ? -1 : ((1ULL << digit_size) - 1));
  uint64_t last_digit_mask =
      (last_digit_size == 64 ? -1 : ((1ULL << last_digit_size) - 1));
  if (last_digit_size == 0) {
    last_digit_mask = digit_mask;
    last_digit_size = digit_size;
  }
  int32_t num_digits = ceil((bw_x * 1.0) / digit_size);
  uint64_t *x_digits = new uint64_t[num_digits * size];

  XTProtocol *xt = new XTProtocol(this->party, this->io, this->otpack);

  // Extract digits
  this->digit_decomposition_sci(size, x, x_digits, bw_x, digit_size);

  // Use LUTs for MSNZB on digits
  int D = (1 << digit_size);
  int DLast = (1 << last_digit_size);
  uint8_t *z_ = new uint8_t[num_digits * size];
  uint64_t *msnzb_ = new uint64_t[num_digits * size];
  uint64_t *msnzb_extended = new uint64_t[num_digits * size];
  int lookup_output_bits = (ceil(log2(digit_size))) + 1;
  int mux_bits = ceil(log2(bw_x));
  uint64_t msnzb_mask = (1ULL << (lookup_output_bits - 1)) - 1;
  uint64_t mux_mask = (1ULL << mux_bits) - 1;
  if (party == sci::ALICE) {
    uint64_t **spec;
    spec = new uint64_t *[num_digits * size];
    PRG128 prg;
    prg.random_data(z_, size * sizeof(uint8_t));
    prg.random_data(msnzb_, size * sizeof(uint64_t));
    for (int i = 0; i < (num_digits - 1) * size; i++) {
      spec[i] = new uint64_t[D];
      z_[i] &= 1;
      msnzb_[i] &= msnzb_mask;
      for (int j = 0; j < D; j++) {
        int idx = (x_digits[i] + j) & digit_mask;
        uint64_t lookup_val = lookup_msnzb(idx);
        spec[i][j] = ((lookup_val >> 1) - msnzb_[i]) & msnzb_mask;
        spec[i][j] <<= 1;
        spec[i][j] |=
            ((uint64_t)(((uint8_t)(lookup_val & 1ULL)) ^ z_[i]) & 1ULL);
      }
    }
    for (int i = (num_digits - 1) * size; i < num_digits * size; i++) {
      spec[i] = new uint64_t[DLast];
      z_[i] &= 1;
      msnzb_[i] &= msnzb_mask;
      for (int j = 0; j < DLast; j++) {
        int idx = (x_digits[i] + j) & last_digit_mask;
        uint64_t lookup_val = lookup_msnzb(idx);
        spec[i][j] = ((lookup_val >> 1) - msnzb_[i]) & msnzb_mask;
        spec[i][j] <<= 1;
        spec[i][j] |=
            ((uint64_t)(((uint8_t)(lookup_val & 1ULL)) ^ z_[i]) & 1ULL);
      }
    }
    if (last_digit_size == digit_size) {
      this->lookup_table<uint64_t>(spec, nullptr, nullptr, num_digits * size,
                                   digit_size, lookup_output_bits);
    } else {
      this->lookup_table<uint64_t>(spec, nullptr, nullptr,
                                   (num_digits - 1) * size, digit_size,
                                   lookup_output_bits);
      this->lookup_table<uint64_t>(spec + (num_digits - 1) * size, nullptr,
                                   nullptr, size, last_digit_size,
                                   lookup_output_bits);
    }

    // Zero extend to mux_bits
    xt->z_extend(num_digits * size, msnzb_, msnzb_extended,
                 lookup_output_bits - 1, mux_bits);

    for (int i = 0; i < num_digits * size; i++) {
      delete[] spec[i];
    }
    delete[] spec;
  } else {  // BOB
    if (last_digit_size == digit_size) {
      this->lookup_table<uint64_t>(nullptr, x_digits, msnzb_, num_digits * size,
                                   digit_size, lookup_output_bits);
    } else {
      this->lookup_table<uint64_t>(nullptr, x_digits, msnzb_,
                                   (num_digits - 1) * size, digit_size,
                                   lookup_output_bits);
      this->lookup_table<uint64_t>(nullptr, x_digits + (num_digits - 1) * size,
                                   msnzb_ + (num_digits - 1) * size, size,
                                   last_digit_size, lookup_output_bits);
    }

    for (int i = 0; i < (num_digits * size); i++) {
      z_[i] = (uint8_t)(msnzb_[i] & 1ULL);
      msnzb_[i] >>= 1;
    }

    // Zero extend to mux_bits
    xt->z_extend(num_digits * size, msnzb_, msnzb_extended,
                 lookup_output_bits - 1, mux_bits);

    for (int i = 0; i < num_digits; i++) {
      for (int j = 0; j < size; j++) {
        msnzb_extended[i * size + j] += (i * digit_size);
        msnzb_extended[i * size + j] &= mux_mask;
      }
    }
  }

  // Combine MSNZB of digits
  uint8_t *dp_zeros_ = new uint8_t[(num_digits - 1) * size];
  uint8_t *one_xor_zeros_ = new uint8_t[(num_digits - 1) * size];
  uint8_t *dp_zeros_final = new uint8_t[num_digits * size];

  if (party == sci::ALICE) {
    for (int i = 0; i < size; i++) {
      dp_zeros_final[(num_digits - 1) * size + i] =
          z_[(num_digits - 1) * size + i];
    }
    for (int i = 0; i < (num_digits - 1); i++) {
      for (int j = 0; j < size; j++) {
        one_xor_zeros_[i * size + j] = z_[i * size + j];
      }
    }
  } else {
    for (int i = 0; i < size; i++) {
      dp_zeros_final[(num_digits - 1) * size + i] =
          (1 ^ z_[(num_digits - 1) * size + i]);
    }
    for (int i = 0; i < (num_digits - 1); i++) {
      for (int j = 0; j < size; j++) {
        one_xor_zeros_[i * size + j] = (1 ^ z_[i * size + j]);
      }
    }
  }
  for (int i = (num_digits - 2); i >= 0; i--) {
    if (i == (num_digits - 2)) {
      for (int j = 0; j < size; j++) {
        dp_zeros_[i * size + j] = z_[(i + 1) * size + j];
      }
    } else {
      this->AND(dp_zeros_ + (i + 1) * size, z_ + (i + 1) * size,
                dp_zeros_ + i * size, size);
    }
  }
  this->AND(dp_zeros_, one_xor_zeros_, dp_zeros_final, (num_digits - 1) * size);

  uint64_t *msnzb_muxed = new uint64_t[num_digits * size];
  this->multiplexer(dp_zeros_final, msnzb_extended, msnzb_muxed,
                    num_digits * size, mux_bits, mux_bits);

  for (int i = 0; i < size; i++) {
    msnzb_index[i] = 0ULL;
    for (int j = 0; j < num_digits; j++) {
      msnzb_index[i] += msnzb_muxed[j * size + i];
      msnzb_index[i] &= mux_mask;
    }
  }

  delete xt;
  delete[] x_digits;
  delete[] z_;
  delete[] msnzb_;
  delete[] msnzb_extended;
  delete[] dp_zeros_;
  delete[] one_xor_zeros_;
  delete[] dp_zeros_final;
  delete[] msnzb_muxed;
  return;
}

void AuxProtocols::msnzb_one_hot(uint64_t *x, uint8_t *one_hot_vector,
                                 int32_t bw_x, int32_t size,
                                 int32_t digit_size) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  int msnzb_index_bits = ceil(log2(bw_x));
  uint64_t msnzb_index_mask = (1ULL << msnzb_index_bits) - 1;

  uint64_t *msnzb_index = new uint64_t[size];

  this->msnzb_sci(x, msnzb_index, bw_x, size, digit_size);

  // use LUT to get the one-hot representation
  int D = 1 << msnzb_index_bits;
  uint64_t *xor_mask = new uint64_t[size];
  if (party == sci::ALICE) {
    uint64_t **spec;
    spec = new uint64_t *[size];
    PRG128 prg;
    prg.random_data(one_hot_vector, size * bw_x * sizeof(uint8_t));
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < bw_x; j++) {
        one_hot_vector[i * bw_x + j] &= 1;
      }
      xor_mask[i] = 0ULL;
      for (int j = (bw_x - 1); j >= 0; j--) {
        xor_mask[i] <<= 1;
        xor_mask[i] ^= (uint64_t)one_hot_vector[i * bw_x + j];
      }
    }
    for (int i = 0; i < size; i++) {
      spec[i] = new uint64_t[D];
      for (int j = 0; j < D; j++) {
        int idx = (msnzb_index[i] + j) & msnzb_index_mask;
        uint64_t lookup_val = (1ULL << idx);
        lookup_val ^= xor_mask[i];
        spec[i][j] = lookup_val;
      }
    }
    this->lookup_table<uint64_t>(spec, nullptr, nullptr, size, msnzb_index_bits,
                                 bw_x);

    for (int i = 0; i < size; i++) {
      delete[] spec[i];
    }
    delete[] spec;
  } else {  // BOB
    uint64_t *temp = new uint64_t[size];
    this->lookup_table<uint64_t>(nullptr, msnzb_index, temp, size,
                                 msnzb_index_bits, bw_x);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < bw_x; j++) {
        one_hot_vector[i * bw_x + j] = (uint8_t)(temp[i] & 1ULL);
        temp[i] >>= 1;
      }
    }
    delete[] temp;
  }
  delete[] msnzb_index;
}

template void AuxProtocols::lookup_table(uint64_t **spec, uint64_t *x,
                                         uint64_t *y, int32_t size,
                                         int32_t bw_x, int32_t bw_y);
template void AuxProtocols::lookup_table(uint8_t **spec, uint8_t *x, uint8_t *y,
                                         int32_t size, int32_t bw_x,
                                         int32_t bw_y);

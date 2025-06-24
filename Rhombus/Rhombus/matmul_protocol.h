
#ifndef RHOMBUS_MATMUL_PROTOCOL_H
#define RHOMBUS_MATMUL_PROTOCOL_H

#include "utils/emp-tool.h"
#include "matmul.h"
#include "matvec.h"

namespace antchain{
    using namespace matmul;
    using namespace matvec;

    // Require: nrows, ncols <= N
    class RhombusLinear {
    public:
        RhombusLinear(sci::NetIO *io, int party, uint32_t N, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits);

        [[nodiscard]] RhombusMatMul& matvec_impl() const {return *matmul_impl_;}

        void send_keys();

        void recv_keys();

        void send_ciphertext(const seal::Ciphertext &ct);

        void send_ciphertext(const std::string &ct);

        void send_encrypted_vector(const std::vector<std::string> &ct);

        void recv_ciphertext(seal::Ciphertext &ct);

        void send_encrypted_vector(const std::vector<seal::Ciphertext> &ct_vec);

        void recv_encrypted_vector(std::vector<seal::Ciphertext> &ct_vec);

        void MatVecMul(uint32_t nr, uint32_t nc, uint64_t *mat, uint64_t *vec, uint64_t *mv);

        // server holds X, client holds Y, the output is secretly shared
        void MatMul(uint32_t n, uint32_t m, uint32_t k, uint64_t *matX, uint64_t *matY, uint64_t *outMat);

        // Point-wise convolution
        void Conv2D_1m1Filter(uint32_t H, uint32_t W, uint32_t CI, uint32_t CO, uint32_t stride, uint64_t *inputArr,
                              uint64_t *filterArr, uint64_t *outArr);

        [[nodiscard]] int64_t to_signed(uint64_t a) const;

        [[nodiscard]] uint64_t to_unsigned(int64_t a) const;

        void to_signed(const uint64_t *inputArr, size_t len, int64_t *outArr) const;

        void to_unsigned(const int64_t *inputArr, size_t len, uint64_t *outArr) const;

        void SignedMatMulMod(const int64_t *matX, const int64_t *matY, int64_t *matXY, uint32_t n, uint32_t m, uint32_t k) const;

        sci::NetIO * IO() {return io_;}

        [[nodiscard]] uint32_t mod_bits() const {return mod_bits_;}

        void set_method(bool use_PackRLWEs){
            use_PackRLWEs_based_method_ = use_PackRLWEs;
            matmul_impl_->set_method(use_PackRLWEs);
        }

        void set_n_thread(uint32_t n_thrd){threads_ = n_thrd;}

        [[nodiscard]] bool get_method() const {return use_PackRLWEs_based_method_;}

        [[nodiscard]] uint32_t get_remain_mod_num() const {return matmul_impl_->get_remain_mod_num();}

        void set_remain_mod_num(uint32_t n){
            matmul_impl_->set_remain_mod_num(n);
            matvec_impl_->set_remain_mod_num(n);
        }

    private:
        int party_{-1};
        bool use_PackRLWEs_based_method_ = true;
        std::shared_ptr<RhombusMatMul> matmul_impl_ = nullptr;
        std::shared_ptr<RhombusMatVec> matvec_impl_ = nullptr;
        uint32_t mod_bits_;
        uint64_t base_mod_;
        uint64_t half_mod_;
        uint64_t mod_mask_;
        sci::NetIO *io_;
        uint32_t threads_ = 4;
    };
}

#endif // RHOMBUS_MATMUL_PROTOCOL_H

#ifndef RHOMBUS_MATMUL_H
#define RHOMBUS_MATMUL_H

#include "seal/seal.h"
#include "seal_api.h"
#include "matrix.h"

namespace antchain::matmul{

#define THREAD_NUM_MAX 64

//#define RHOMBUS_UNIT_TESTS 1

    using namespace antchain::global;

    class RhombusMatMul{
    public:
        RhombusMatMul();

        RhombusMatMul(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits);

        // Generate keys: secret key, galois key and public key if in public key encryption mode
        // this function is only for test
        void GenKey();

        size_t GenPublicKey(std::string &out) const;

        size_t GenPublicKey(uint8_t *buffer, size_t buffer_size) const;

        size_t GenGaloisKey(std::string &out) const;

        size_t GenGaloisKey(uint8_t *buffer, size_t buffer_size) const;

        size_t GenGaloisKey(std::string &out, const std::vector<uint32_t> &galois_elts) const;

        size_t GenGaloisKey(uint8_t *buffer, size_t buffer_size, const std::vector<uint32_t> &galois_elts) const;

        size_t SetPublicKey(const std::string &in);

        size_t SetPublicKey(const uint8_t *buffer, size_t buffer_size);

        size_t SetGaloisKey(const std::string &in);

        size_t SetGaloisKey(const uint8_t *buffer, size_t buffer_size);

        void SetPublicKey(const seal::PublicKey &pk){
            public_key_ = pk;
        }

        void SetGaloisKey(const seal::GaloisKeys &gk){
            galois_keys_ = gk;
        }


        [[nodiscard]] const seal::SecretKey & secret_key() const{
            return secret_key_;
        }

        seal::PublicKey &public_key() {return public_key_;}

        seal::GaloisKeys &galois_key() {return galois_keys_;}

        // reset the secret key, and update the related members
        void reset_secret_key(const seal::SecretKey &new_sk);

        void SetSPPMap(const std::vector<int> &new_spp_table){
            kSPPMap_ = new_spp_table;
        }

        [[nodiscard]] const seal::SEALContext & seal_context() const{
            return *seal_context_;
        }

        [[nodiscard]] const AuxParms & get_aux_parms() const{
            return aux_parms_;
        }

        void reset_mod_bits(uint32_t mod_bits);

        // default: Expand-based method
        void set_method(bool use_PackRLWEs_based) {
            use_PackRLWEs_based_MatMul_ = use_PackRLWEs_based;
        }

        void set_remain_mod_num(uint32_t re_mn){
            remain_mod_num_ = re_mn;
        }

        uint32_t get_remain_mod_num() const {return remain_mod_num_;}

        // REQUIRED!!!
        void SetMatDims(uint32_t n, uint32_t m, uint32_t k, uint32_t X_bits, uint32_t Y_bits);

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                EncodeMatX_P(matrix, encoded_mat, threads);
            else
                EncodeMatX_E(matrix, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX(const std::vector<std::vector<T>> &matrix,
                        std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                EncodeMatX_P(matrix, encoded_mat, threads);
            else
                EncodeMatX_E(matrix, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                EncodeMatY_P(matrix, encoded_mat, threads);
            else
                EncodeMatY_E(matrix, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY(const std::vector<std::vector<T>> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                        uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                EncodeMatY_P(matrix, encoded_mat, threads);
            else
                EncodeMatY_E(matrix, encoded_mat, threads);
        }


        void EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat, std::vector<seal::Ciphertext> &encrypted_mat,
                           uint32_t threads = 4, bool is_symmetric = true) const;

        void EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat, std::vector<std::string> &encrypted_mat,
                           uint32_t threads = 4, bool is_symmetric = true) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const T *matrix, std::vector<std::string> &encrypted_mat, uint32_t threads = 4,
                         bool is_symmetric = true) const
        {
            // check
            if (matrix == nullptr)
                throw std::invalid_argument("null pointer of matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, encrypted_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const T *matrix, std::vector<seal::Ciphertext> &encrypted_mat, uint32_t threads = 4,
                         bool is_symmetric = true) const
        {
            // check
            if (matrix == nullptr)
                throw std::invalid_argument("null pointer of matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, encrypted_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const std::vector<std::vector<T>> &matrix, std::vector<seal::Ciphertext> &enc_mat,
                         uint32_t threads = 4, bool is_symmetric = true) const{
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, enc_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const std::vector<std::vector<T>> &matrix, std::vector<std::string> &enc_mat,
                         uint32_t threads = 4, bool is_symmetric = true) const
        {
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, enc_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatMul(const T* matX, const std::vector<seal::Ciphertext> &enc_matY,
                    std::vector<seal::Ciphertext> &result, uint32_t threads = 4) const
        {
            if (matX == nullptr)
                throw std::invalid_argument("empty matrix, nullptr");

            std::vector<seal::Plaintext> encoded_matX;
            EncodeMatX(matX, encoded_matX, threads);
            Compute(encoded_matX, enc_matY, result, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatMul(const std::vector<std::vector<T>> &matX, const std::vector<seal::Ciphertext> &enc_matY,
                    std::vector<seal::Ciphertext> &result, uint32_t threads = 4) const
        {
            if (matX.empty())
                throw std::invalid_argument("empty matrix");
            std::vector<seal::Plaintext> encoded_matX;
            EncodeMatX(matX, encoded_matX, threads);
            Compute(encoded_matX, enc_matY, result, threads);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatMulToSS(const T *matX, const std::vector<seal::Ciphertext> &enc_matY,
                        std::vector<seal::Ciphertext> &enc_matmul_share0, U *matmul_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> matmul;
            MatMul(matX, enc_matY, matmul, threads);
            H2A(matmul, enc_matmul_share0, matmul_share1, threads);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatMulToSS(const std::vector<std::vector<T>> &matX, const std::vector<seal::Ciphertext> &enc_matY,
                        std::vector<seal::Ciphertext> &enc_matmul_share0,  std::vector<std::vector<U>> &matmul_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> matmul;
            MatMul(matX, enc_matY, matmul, threads);
            H2A(matmul, enc_matmul_share0, matmul_share1, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void H2A(const std::vector<seal::Ciphertext> &encrypted_matXY,
                 std::vector<seal::Ciphertext> &encrypted_matXY_share0, std::vector<std::vector<T>> &matXY_share1,
                 uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                ConvToSS_P(encrypted_matXY, encrypted_matXY_share0, matXY_share1, threads);
            else
                ConvToSS_E(encrypted_matXY, encrypted_matXY_share0, matXY_share1, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void H2A(const std::vector<seal::Ciphertext> &encrypted_matXY,
                 std::vector<seal::Ciphertext> &encrypted_matXY_share0, T *matXY_share1,
                 uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                ConvToSS_P(encrypted_matXY, encrypted_matXY_share0, matXY_share1, threads);
            else
                ConvToSS_E(encrypted_matXY, encrypted_matXY_share0, matXY_share1, threads);
        }

        void Compute(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                     std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptMatXY(const std::vector<seal::Ciphertext> &enc_mat, T *result, uint32_t threads = 4) const
        {
            if (result == nullptr)
                throw std::invalid_argument("null pointer");

            std::vector<seal::Plaintext> dec_pt;
            DecryptMatXY(enc_mat, dec_pt, threads);
            DecodeMatXY(dec_pt, result, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY,
                          std::vector<std::vector<T>> &result, uint32_t threads = 4) const
        {
            std::vector<seal::Plaintext> dec_pt;
            DecryptMatXY(encrypted_matXY, dec_pt, threads);
            DecodeMatXY(dec_pt, result, threads);
        }

        void DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY, std::vector<seal::Plaintext> &encoded_matXY,
                          uint32_t threads = 4) const;

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY(const std::vector<seal::Plaintext> &encoded_matXY, T *matXY, uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                DecodeMatXY_P(encoded_matXY, matXY, threads);
            else
                DecodeMatXY_E(encoded_matXY, matXY, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY(const std::vector<seal::Plaintext> &encoded_matXY,
                         std::vector<std::vector<T>> &matXY, uint32_t threads = 4) const
        {
            if (use_PackRLWEs_based_MatMul_)
                DecodeMatXY_P(encoded_matXY, matXY, threads);
            else
                DecodeMatXY_E(encoded_matXY, matXY, threads);
        }

        void drop_unused_coeffs_P(std::vector<seal::Ciphertext> &ct, uint32_t threads = 4) const;

        // (de)serialization
        size_t CiphertextsToBytes(const std::vector<std::vector<seal::Ciphertext>> &ct, std::vector<std::vector<std::string>> &out) const;
        size_t CiphertextsToBytes(const std::vector<seal::Ciphertext> &ct, std::vector<std::string> &out) const;
        size_t BytesToCiphertexts(const std::vector<std::string> &ct_str, std::vector<seal::Ciphertext> &ct_out) const;
        size_t BytesToCiphertexts(const std::vector<std::vector<std::string>> &ct_str, std::vector<std::vector<seal::Ciphertext>> &ct_out) const;

    private:
        // Partition matrix X horizontally into X_0, X_1, ..., then get the pointers of index-th block
        template <typename T>
        void GetSubMat(const std::vector<const T *> &mat, const MatInfoLW &X_info, std::vector<const T *> &sub_mat, size_t index) const
        {
            size_t pad_ncols = size_t(1) << X_info.log_pad_ncols_;

            // ceil(nrows / pad_ncols): the number of sub-matrices
            size_t submat_num = (X_info.nrows_ + pad_ncols - 1) / pad_ncols;

            // Note: we split X into several sub-matrices of dimension 'pad_cols * ncols'
            // For example, '200 * 120' --> 2 sub-matrices: '128 * 120', '72 * 120'
            //              '1000 * 100' --> 8 sub-matrices: the first 7 matrices are of size '128 * 100', the last '104 * 100'
            // To sum up, [1, submat_num-1]-th sub-matrix is of dimension 'pad_ncols * ncols', the last sub-matrix is the remaining part

            // out of range
            if (index >= submat_num)
                throw std::invalid_argument("out of range");

            size_t row_bgn = index * pad_ncols;
            size_t row_end = std::min<size_t>(row_bgn + pad_ncols, X_info.nrows_);

            sub_mat.resize(row_end - row_bgn);
            for (size_t i = row_bgn; i < row_end; ++i){
                sub_mat[i - row_bgn] = mat[i];
            }
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX_E(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("matrix nullptr");

            std::vector<const T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matrix + i * mat_X_.ncols_;

            EncodeMatXInternal_E(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX_E(const std::vector<std::vector<T>> &matrix,
                          std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<const T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matrix[i].data();

            EncodeMatXInternal_E(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX_P(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("matrix nullptr");

            std::vector<const T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matrix + i * mat_X_.ncols_;

            EncodeMatXInternal_P(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX_P(const std::vector<std::vector<T>> &matrix,
                          std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<const T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matrix[i].data();

            EncodeMatXInternal_P(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY_E(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            // check
            if (matrix == nullptr)
                throw std::invalid_argument("nullptr of matrix");

            std::vector<const T *> mat(mat_Y_.nrows_);
            for (size_t i = 0; i < mat_Y_.nrows_; ++i)
                mat[i] = matrix + i * mat_Y_.ncols_;

            EncodeMatYInternal_E(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY_E(const std::vector<std::vector<T>> &matrix,
                          std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<const T *> mat(mat_Y_.nrows_);
            for (size_t i = 0; i < mat_Y_.nrows_; ++i)
                mat[i] = matrix[i].data();

            EncodeMatYInternal_E(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY_P(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("nullptr of matrix");

            std::vector<const T *> mat(mat_Y_.nrows_);
            for (size_t i = 0; i < mat_Y_.nrows_; ++i)
                mat[i] = matrix + i * mat_Y_.ncols_;

            EncodeMatYInternal_P(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY_P(const std::vector<std::vector<T>> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                          uint32_t threads = 4) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<const T *> mat(mat_Y_.nrows_);
            for (size_t i = 0; i < mat_Y_.nrows_; ++i)
                mat[i] = matrix[i].data();

            EncodeMatYInternal_P(mat, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSS_E(const std::vector<seal::Ciphertext> &encrypted_matXY,
                        std::vector<seal::Ciphertext> &encrypted_matXY_share0, T *matXY_share1, uint32_t threads = 4) const
        {
            if (matXY_share1 == nullptr)
                throw std::invalid_argument("matXY_share1 nullptr");

            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matXY_share1 + i * mat_Y_.ncols_;

            ConvToSSInternal_E(encrypted_matXY, encrypted_matXY_share0, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSS_E(const std::vector<seal::Ciphertext> &encrypted_matXY,
                        std::vector<seal::Ciphertext> &encrypted_matXY_share0, std::vector<std::vector<T>> &matXY_share1,
                        uint32_t threads = 4) const
        {
            matXY_share1.resize(mat_X_.nrows_);
            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i){
                matXY_share1[i].resize(mat_Y_.ncols_);
                mat[i] = matXY_share1[i].data();
            }

            ConvToSSInternal_E(encrypted_matXY, encrypted_matXY_share0, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSS_P(const std::vector<seal::Ciphertext> &encrypted_matXY,
                        std::vector<seal::Ciphertext> &encrypted_matXY_share0, T *matXY_share1, uint32_t threads = 4) const
        {
            if (matXY_share1 == nullptr)
                throw std::invalid_argument("matXY_share1 nullptr");

            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matXY_share1 + i * mat_Y_.ncols_;

            ConvToSSInternal_P(encrypted_matXY, encrypted_matXY_share0, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSS_P(const std::vector<seal::Ciphertext> &encrypted_matXY,
                        std::vector<seal::Ciphertext> &encrypted_matXY_share0, std::vector<std::vector<T>> &matXY_share1,
                        uint32_t threads = 4) const
        {
            matXY_share1.resize(mat_X_.nrows_);
            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i){
                matXY_share1[i].resize(mat_Y_.ncols_);
                mat[i] = matXY_share1[i].data();
            }

            ConvToSSInternal_P(encrypted_matXY, encrypted_matXY_share0, mat, threads);
        }

        void Compute_E(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                       std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        void Compute_P(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                       std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY_E(const std::vector<seal::Plaintext> &encoded_matXY, T *matXY, uint32_t threads = 4) const
        {
            if (matXY == nullptr)
                throw std::invalid_argument("nullptr");

            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matXY + i * mat_Y_.ncols_;

            DecodeMatXYInternal_E(encoded_matXY, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY_E(const std::vector<seal::Plaintext> &encoded_matXY,
                           std::vector<std::vector<T>> &matXY, uint32_t threads = 4) const
        {
            matXY.resize(mat_X_.nrows_);
            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i){
                matXY[i].resize(mat_Y_.ncols_);
                mat[i] = matXY[i].data();
            }

            DecodeMatXYInternal_E(encoded_matXY, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY_P(const std::vector<seal::Plaintext> &encoded_matXY, T *matXY, uint32_t threads = 4) const
        {
            if (matXY == nullptr)
                throw std::invalid_argument("nullptr");

            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i)
                mat[i] = matXY + i * mat_Y_.ncols_;

            DecodeMatXYInternal_P(encoded_matXY, mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXY_P(const std::vector<seal::Plaintext> &encoded_matXY,
                           std::vector<std::vector<T>> &matXY, uint32_t threads = 4) const
        {
            matXY.resize(mat_X_.nrows_);
            std::vector<T *> mat(mat_X_.nrows_);
            for (size_t i = 0; i < mat_X_.nrows_; ++i){
                matXY[i].resize(mat_Y_.ncols_);
                mat[i] = matXY[i].data();
            }

            DecodeMatXYInternal_P(encoded_matXY, mat, threads);
        }

        void hom_inner_product_E(std::vector<seal::Ciphertext> &temp, const seal::Plaintext *encoded_matX,
                                 const std::vector<std::vector<seal::Ciphertext>> &cached_ct, uint32_t threads = 4) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeSubMatXInternal(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols,
                                   seal::Plaintext *encoded_mat, uint32_t mat_bits, uint32_t threads) const
        {
            std::vector<std::vector<T>> pad_mat;
            matrix_col_padding(matrix, pad_mat, nrows, ncols);

            // determine if the sub-matrix is the last block of X
            size_t g0_size, g1_size;
            size_t spp_index = ((pad_mat.size() < ncols) ? 1 : 0);
            g0_size = (uint32_t)1 << (spp_PackRLWEs_[spp_index].u_ - spp_PackRLWEs_[spp_index].h_);
            g1_size = (uint32_t)1 << (spp_PackRLWEs_[spp_index].ell_ + spp_PackRLWEs_[spp_index].h_ - spp_PackRLWEs_[spp_index].u_);

            // After matrix-preprocessing, the elements will become larger
            uint32_t max_bits = spp_PackRLWEs_[spp_index].u_ - spp_PackRLWEs_[spp_index].h_ + mat_bits;

            auto encode_program = [&](size_t bgn, size_t end){
                if (max_bits <= sizeof(T) * 8){
                    std::vector<std::vector<std::vector<T>>> pp_mat;
                    global::matrix_row_combine64(pad_mat, spp_PackRLWEs_[spp_index], pp_mat, threads,
                                                 poly_modulus_degree_, ncols);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < g0_size; ++j){
                            encode_to_coeff(encoded_mat[i * g0_size + j], pp_mat[i][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER,
                                            aux_parms_, seal_context_->first_parms_id(), *seal_context_, global::kSSBitLen /*no use*/);
                        }
                    }
                }
                else if (max_bits <= 64){
                    std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                    global::matrix_row_combine64(pad_mat, spp_PackRLWEs_[spp_index], pp_mat, threads,
                                                 poly_modulus_degree_, ncols);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < g0_size; ++j){
                            encode_to_coeff(encoded_mat[i * g0_size + j], pp_mat[i][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER,
                                            aux_parms_, seal_context_->first_parms_id(), *seal_context_, global::kSSBitLen /* no use*/);
                        }
                    }
                }
                else if (max_bits <= 128){
                    std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                    global::matrix_row_combine128(pad_mat, spp_PackRLWEs_[spp_index], pp_mat, threads,
                                                  poly_modulus_degree_, ncols);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < g0_size; ++j){
                            encode_to_coeff128(encoded_mat[i * g0_size + j], pp_mat[i][j].data(), poly_modulus_degree_, aux_parms_,
                                               seal_context_->first_parms_id(), *seal_context_);
                        }
                    }
                }
            };

            uint32_t thread_block = (g1_size + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, g1_size);
                thread_pool[i] = std::thread(encode_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        // compute tau(ct) for all tau \in Gal(K_u/K_h)
        void hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                  const global::SPPParams &spp, uint32_t threads, bool remove_pack_factor = true) const;

        // DO NOT support multi-thread now. In most cases, there's no need to add multi-thread on this method.
        // The input cts will be modified !!!
        void hom_aut_and_add_group(std::vector<seal::Ciphertext> &cts, seal::Ciphertext &result,
                                   const global::SPPParams_Expand &spp, uint32_t threads) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatXInternal_P(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                                  uint32_t threads = 4) const
        {
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");
            if (threads == 0 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("invalid thread number");

            // partition windows
            size_t mw = size_t(1) << mat_X_.log_pad_ncols_;

            // ceil(n / mw), ceil(k / kw)
            size_t submatX_num = (mat_X_.nrows_ + mw - 1) / mw;

            uint32_t stride = 1ULL << spp_PackRLWEs_[0].ell_;

            // an upper bound
            encoded_mat.resize(submatX_num * stride);

            if (submatX_num >= threads){
                auto ecd_submat = [&](size_t bgn, size_t end){
                    std::vector<const T *> submatX;
                    for (size_t i = bgn; i < end; ++i){
                        GetSubMat(matrix, mat_X_, submatX, i);
                        EncodeSubMatXInternal(submatX, submatX.size(), mat_X_.ncols_,
                                              encoded_mat.data() + i * stride, mat_X_.mat_bits_, 1);
                    }
                };

                uint32_t thread_block = (submatX_num + threads - 1) / threads;
                std::vector<std::thread> thread_pool(threads);
                for (size_t i = 0; i < threads; ++i){
                    size_t bgn = i * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, submatX_num);
                    thread_pool[i] = std::thread(ecd_submat, bgn, end);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }
            else {
                std::vector<const T *> submatX;
                for (size_t i = 0; i < submatX_num; ++i){
                    GetSubMat(matrix, mat_X_, submatX, i);
                    EncodeSubMatXInternal(submatX, submatX.size(), mat_X_.ncols_,
                                          encoded_mat.data() + i * stride, mat_X_.mat_bits_, threads);
                }
            }
        }

        void MatMulBlock(const seal::Plaintext *ecd_submatX, const std::vector<seal::Ciphertext> &cached_subY,
                         seal::Ciphertext &enc_submatXY, const SPPParams &spp, uint32_t threads = 4, bool mul_factor = false) const;

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXYInternal_E(const std::vector<seal::Plaintext> &encoded_matXY, std::vector<T *> &matXY, uint32_t threads = 4) const
        {
            if (encoded_matXY.empty())
                throw std::invalid_argument("empty plaintexts");
            if (threads == 0 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num. invalid");

            uint32_t pad_k = 1UL << mat_Y_.log_pad_ncols_;

            uint32_t n = mat_X_.nrows_;
            uint32_t n_partition = poly_modulus_degree_ / pad_k;
            uint32_t nblock = (n + n_partition - 1) / n_partition;
            uint32_t last_block_nrows = n - n_partition * (nblock - 1);

            if (nblock != encoded_matXY.size())
                throw std::invalid_argument("ct num. mismatch");

            auto dcd_program = [&](size_t bgn, size_t end){
                std::vector<T> vec(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    uint32_t cur_block_nrows = n_partition;
                    uint32_t cur_block_ncols = mat_Y_.ncols_;
                    if (i == nblock - 1){
                        cur_block_nrows = last_block_nrows;
                    }
                    decode_from_coeff(vec.data(), poly_modulus_degree_, encoded_matXY[i], Dcd_SCALED_IN_ORDER,
                                      aux_parms_, *seal_context_, mod_bits_);
                    uint32_t row_id_start = i * n_partition;

                    for (size_t j = 0; j < cur_block_ncols; ++j){
                        for (size_t k = 0; k < cur_block_nrows; ++k)
                            matXY[k + row_id_start][j] = vec[k + j * n_partition];
                    }
                }
            };

            uint32_t thread_block = (nblock + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, nblock);
                thread_pool[i] = std::thread(dcd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXYInternal_P(const std::vector<seal::Plaintext> &encoded_matXY, std::vector<T *> &matXY, uint32_t threads = 4) const
        {
            // check
            if (encoded_matXY.empty())
                throw std::invalid_argument("empty matrix");

            size_t n = mat_X_.nrows_;
            size_t k = mat_Y_.ncols_;
            size_t mw = (size_t)1 << BitLength(mat_Y_.nrows_ - 1);
            size_t kw = poly_modulus_degree_ / mw;
            size_t nw = mw;

            // ceil(n / mw), ceil(k / kw)
            size_t row_block_num = (mat_X_.nrows_ + mw - 1) / mw;
            size_t col_block_num = (mat_Y_.ncols_ + kw - 1) / kw;

            size_t last_rblock_nrow = n - (row_block_num - 1) * nw;
            size_t last_cblock_ncol = k - (col_block_num - 1) * kw;

            auto get_sub_mat = [&](std::vector<T *> &sub_mat, size_t rblk_index, size_t cblk_index)
            {
                size_t submat_nrows;
                if (rblk_index == row_block_num - 1)
                    submat_nrows = last_rblock_nrow;
                else
                    submat_nrows = nw;

                sub_mat.resize(submat_nrows);
                for (size_t i = 0; i < submat_nrows; ++i){
                    sub_mat[i] = matXY[rblk_index * nw + i] + cblk_index * kw;
                }
            };

            auto dcd_block = [&](const seal::Plaintext &pt, std::vector<T *> &result, uint32_t nrows,
                    uint32_t ncols){
                std::vector<T> dcd_vec(poly_modulus_degree_);
                decode_from_coeff(dcd_vec.data(), poly_modulus_degree_, pt, Dcd_SCALED_IN_ORDER,
                                  aux_parms_, *seal_context_, mod_bits_);

                // stride
                size_t col_stride = poly_modulus_degree_ / kw;
                size_t pow2_nrows = 1ULL << BitLength(nrows - 1);
                size_t row_stride = col_stride / pow2_nrows;

                for (size_t i = 0; i < nrows; ++i){
                    for (size_t j = 0; j < ncols; ++j){
                        result[i][j] = dcd_vec[i * row_stride + j * col_stride];
                    }
                }
            };

            // decode
            auto dcd_program = [&](size_t bgn, size_t end){
                std::vector<T *> sub_mat;
                for (size_t j = bgn; j < end; ++j){
                    for (size_t i = 0; i < row_block_num; ++i){
                        get_sub_mat(sub_mat, i, j);
                        size_t nrows, ncols;
                        if (i == (row_block_num - 1))
                            nrows = last_rblock_nrow;
                        else
                            nrows = nw;
                        if (j == (col_block_num - 1))
                            ncols = last_cblock_ncol;
                        else
                            ncols = kw;
                        dcd_block(encoded_matXY[i * col_block_num + j], sub_mat, nrows, ncols);
                    }
                }
            };

            uint32_t thread_block = (col_block_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t j = 0; j < threads; ++j){
                size_t bgn = j * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, col_block_num);
                thread_pool[j] = std::thread(dcd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatYInternal_E(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num is not valid");
            if (matrix.size() != mat_Y_.nrows_)
                throw std::invalid_argument("matrix dimension mismatches");

            // Y: m * k, pad k to 2-power
            uint32_t pad_k = (uint32_t)1 << mat_Y_.log_pad_ncols_;

            // every mw rows will be packed into a plaintext
            size_t mw = poly_modulus_degree_ / pad_k;

            // Y will be encoded to pt_num = ceil(m/mw) plaintexts
            size_t pt_num = (mat_Y_.nrows_ + mw - 1) / mw;
            encoded_mat.resize(pt_num);

            // encode pt: [bgn, end) i-th pt will encode i*mw ~ (i+1)*mw - 1 rows of Y
            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    std::fill_n(temp.data(), poly_modulus_degree_, 0);
                    size_t row_bgn = i * mw;
                    size_t row_end = std::min<size_t>(row_bgn + mw, mat_Y_.nrows_);

                    // We use Ecd_1 to encode each block of Y:
                    // for a block of dimension mw * pad_k, the columns of this block will be concatenated to a length-N vector
                    for (size_t col_id = 0; col_id < mat_Y_.ncols_; ++col_id){ // for each column
                        for (size_t row_id = row_bgn; row_id < row_end; ++row_id){
                            temp[col_id * mw + row_id - row_bgn] = matrix[row_id][col_id];
                        }
                    }
                    // encode
                    antchain::encode_to_coeff(encoded_mat[i], temp.data(), poly_modulus_degree_, Ecd_SCALED_IN_ORDER,
                                              aux_parms_, seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                }
            };

            uint32_t thread_block = (pt_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, pt_num);
                thread_pool[i] = std::thread(ecd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatYInternal_P(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                                  uint32_t threads = 4) const
        {
            if (threads < 1 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num is not valid");

            // Y: m * k. pad m to 2-power number
            uint32_t pad_m = (uint32_t)1 << mat_Y_.log_pad_nrows_;

            // every kw columns will be packed into a plaintext
            size_t kw = poly_modulus_degree_ / pad_m;

            // Y will be encoded to pt_num = ceil(k/kw) plaintexts
            size_t pt_num = (mat_Y_.ncols_ + kw - 1) / kw;
            encoded_mat.resize(pt_num);

            // encode pt: [bgn, end), i-th pt will encode i*kw ~ (i+1)*kw-1 columns
            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    std::fill_n(temp.data(), poly_modulus_degree_, 0);
                    size_t col_bgn = i * kw;
                    size_t col_end = std::min<size_t>(col_bgn + kw, mat_Y_.ncols_);

                    // We use Ecd_2 to encode each column of Y
                    // The 0-th column: (a0, a1, ..., a_{m-1}) --> a0 - a1*X^{N-1} - ...- a_{m-1}*X^{N-m-1}
                    temp[0] = matrix[0][col_bgn];
                    for (size_t rows_id = 1; rows_id < mat_Y_.nrows_; ++rows_id){
                        temp[poly_modulus_degree_ - rows_id] = -matrix[rows_id][col_bgn];
                    }

                    // Other columns:
                    // i-th column (1 <= i < kw): (c0, c1, ..., c_{m-1}) --> c_{m-1}*X^{i*pad_m-m+1} + c_{m-2}*X^{i*pad_m-m+2} + ... + c0*X^{i*pad_m}
                    for (size_t j = col_bgn + 1; j < col_end; ++j){
                        for (size_t rows_id = 0; rows_id < mat_Y_.nrows_; ++rows_id){
                            temp[(j - col_bgn) * pad_m - rows_id] = matrix[rows_id][j];
                        }
                    }

                    // encode the block (pad_m * kw), with scale q/2^mod
                    antchain::encode_to_coeff(encoded_mat[i], temp.data(), poly_modulus_degree_, Ecd_SCALED_IN_ORDER, aux_parms_,
                                              seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                }
            };

            // multi-thread
            uint32_t thread_block = (pt_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, pt_num);
                thread_pool[i] = std::thread(ecd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        void EncryptMatYInternal(const std::vector<seal::Plaintext> &encoded_mat, seal::Ciphertext *encrypted_mat,
                                   std::string *serialized_enc_mat, uint32_t threads = 4, bool is_symmetric = true) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeSubMatX(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols, uint32_t block_size,
                           seal::Plaintext *ecd_mat, uint32_t mat_bits, uint32_t threads) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            uint32_t logN = BitLength(poly_modulus_degree_ - 1);
            // maximum possible value
            uint32_t max_bits = spp_Expand_.u_ + spp_Expand_.ell_ - logN + mat_bits;

            size_t g0_size, g1_size;
            g0_size = 1ULL << (spp_Expand_.u_ + spp_Expand_.ell_ - logN); // c*2^u/N
            g1_size = poly_modulus_degree_ >> spp_Expand_.u_; // N/2^u

            // matrix will be divided vertically into 'nblock' sub-matrix
            uint32_t nblock = (ncols + block_size - 1) / block_size;

            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<const T *> block(nrows);
                for (size_t i = bgn; i < end; ++i){
                    // get block
                    for (size_t j = 0; j < nrows; ++j){
                        block[j] = matrix[j] + i * block_size;
                    }
                    uint32_t cur_block_ncols = (i == nblock - 1) ? (ncols - i * block_size) : block_size;
                    if (max_bits <= sizeof(T) * 8){
                        std::vector<std::vector<std::vector<T>>> pp_mat;
                        global::matrix_col_combine64(block, spp_Expand_, cur_block_ncols, block_size, pp_mat,
                                                     1, poly_modulus_degree_);

                        for (size_t i1 = 0; i1 < g1_size; ++i1){
                            for (size_t j = 0; j < g0_size; ++j){
                                encode_to_coeff(ecd_mat[i * g1_size * g0_size + i1 * g0_size + j], pp_mat[i1][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER,
                                                aux_parms_, seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                            }
                        }
                    }else if (max_bits <= 64){
                        std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                        global::matrix_col_combine64(block, spp_Expand_, cur_block_ncols, block_size, pp_mat,
                                                     1, poly_modulus_degree_);
                        for (size_t i1 = 0; i1 < g1_size; ++i1){
                            for (size_t j = 0; j < g0_size; ++j){
                                encode_to_coeff(ecd_mat[i * g1_size * g0_size + i1 * g0_size + j], pp_mat[i1][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER,
                                                aux_parms_, seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                            }
                        }
                    }else{
                        std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                        global::matrix_col_combine128(block, spp_Expand_, cur_block_ncols, block_size,
                                                      pp_mat, 1, poly_modulus_degree_);
                        for (size_t i1 = 0; i1 < g1_size; ++i1){
                            for (size_t j = 0; j < g0_size; ++j){
                                encode_to_coeff128(ecd_mat[i * g1_size * g0_size + i1 * g0_size + j], pp_mat[i1][j].data(), poly_modulus_degree_, aux_parms_,
                                                   seal_context_->first_parms_id(), *seal_context_);
                            }
                        }
                    }
                }
            };

            uint32_t thread_block = (nblock + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, nblock);
                thread_pool[i] = std::thread(ecd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatXInternal_E(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                                  uint32_t threads = 4) const
        {

            uint32_t block_size = poly_modulus_degree_ >> mat_Y_.log_pad_ncols_;
            uint32_t nblock = (mat_X_.ncols_ + block_size - 1) / block_size;
            uint32_t logN = BitLength(poly_modulus_degree_ - 1);
            uint32_t g0_size = 1UL << (spp_Expand_.u_ + spp_Expand_.ell_ - logN);
            uint32_t g1_size = 1UL << (logN - spp_Expand_.u_);
            uint32_t submatX_num = (mat_X_.nrows_ + block_size - 1) / block_size;

            encoded_mat.resize(submatX_num * nblock * g1_size * g0_size);
            uint32_t stride = nblock * g1_size * g0_size;

            auto encode_slice = [&](size_t bgn, size_t end){
                std::vector<const T *> submatX;
                for (size_t i = bgn; i < end; ++i){
                    uint32_t nsubmat_row = (i != submatX_num - 1) ? block_size : (mat_X_.nrows_ - i * block_size);
                    submatX.resize(nsubmat_row);

                    // get submat
                    for (size_t j = 0; j < nsubmat_row; ++j){
                        submatX[j] = matrix[i * block_size + j];
                    }

                    // encode submat
                    EncodeSubMatX(submatX, nsubmat_row, mat_X_.ncols_, block_size,
                                  encoded_mat.data() + i * stride, mat_X_.mat_bits_, 1);
                }
            };

            uint32_t thread_block = (submatX_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, submatX_num);
                thread_pool[i] = std::thread(encode_slice, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSSInternal_E(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                                std::vector<T *> &matXY_share1, uint32_t threads = 4) const
        {
            auto coeff_modulus = seal_context_->get_context_data(encrypted_matXY[0].parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::parms_id_type parms_id = encrypted_matXY[0].parms_id();
            std::mt19937_64 gen(std::random_device{}());

            auto gen_rand_pt = [&](seal::Plaintext &plain){
                plain.parms_id() = seal::parms_id_zero;
                plain.resize(poly_modulus_degree_ * coeff_modulus_size);
                for (size_t k = 0; k < coeff_modulus_size; ++k){
                    std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[k].value() - 1);
                    std::generate_n(plain.data() + k * poly_modulus_degree_,
                                    poly_modulus_degree_, [&](){return dist(gen);});
                }
                plain.parms_id() = parms_id;
                plain.scale() = 1.;
            };

            size_t ct_num = encrypted_matXY.size();
            encrypted_matXY_share0.resize(ct_num);
            std::vector<seal::Plaintext> matXY_pt_share1(ct_num);

            auto H2A_program = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    gen_rand_pt(matXY_pt_share1[i]);
                    sub_plain(encrypted_matXY[i], matXY_pt_share1[i], encrypted_matXY_share0[i], *seal_context_);
                }
            };

            uint32_t thread_block = (threads - 1 + ct_num) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, ct_num);
                thread_pool[i] = std::thread(H2A_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            DecodeMatXYInternal_E(matXY_pt_share1, matXY_share1, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSSInternal_P(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                                std::vector<T *> &matXY_share1, uint32_t threads = 4) const
        {
            auto coeff_modulus = seal_context_->get_context_data(encrypted_matXY[0].parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::parms_id_type parms_id = encrypted_matXY[0].parms_id();
            std::mt19937_64 gen(std::random_device{}());

            uint32_t ct_num = encrypted_matXY.size();

            auto gen_rand_pt = [&](seal::Plaintext &plain){
                plain.parms_id() = seal::parms_id_zero;
                plain.resize(poly_modulus_degree_ * coeff_modulus_size);
                for (size_t k = 0; k < coeff_modulus_size; ++k){
                    std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[k].value() - 1);
                    std::generate_n(plain.data() + k * poly_modulus_degree_, poly_modulus_degree_, [&](){return dist(gen);});
                }
                plain.parms_id() = parms_id;
                plain.scale() = 1.;
            };

            std::vector<seal::Plaintext> matXY_pt_share1(ct_num);
            encrypted_matXY_share0.resize(ct_num);

            auto H2A_program = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    gen_rand_pt(matXY_pt_share1[i]);
                    sub_plain(encrypted_matXY[i], matXY_pt_share1[i], encrypted_matXY_share0[i], *seal_context_);
                }
            };

            uint32_t thread_block = (threads - 1 + ct_num) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, ct_num);
                thread_pool[i] = std::thread(H2A_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

            DecodeMatXYInternal_P(matXY_pt_share1, matXY_share1, threads);
        }

        bool use_PackRLWEs_based_MatMul_ = false;

        uint32_t poly_modulus_degree_;
        uint32_t mod_bits_;

        uint32_t remain_mod_num_ = 1;

        std::shared_ptr<seal::SEALContext> seal_context_;
        std::unique_ptr<seal::Decryptor> decryptor_;
        std::unique_ptr<seal::KeyGenerator> keygen_;

        seal::SecretKey secret_key_;
        seal::PublicKey my_public_key_;

        // other party's pk, gk
        seal::PublicKey public_key_;
        seal::GaloisKeys galois_keys_;

        // X * Enc(Y) The dimension of X, Y < N * N.
        MatInfoLW mat_X_;
        MatInfoLW mat_Y_;

        // Note: X_{n, m}, will be partitioned into t = ceil(n / pad_m) blocks, pad_m = 2^{ceil(log m)}
        // 1 ~ (t-1)-th blocks are of dimension pad_m * m, which correspond to spp_X_[0]
        // the t-th (last) block is of dimension [n - (pad_m * (t-1))] * m, corresponding to spp_X_[1]
        SPPParams spp_PackRLWEs_[2];

        SPPParams_Expand spp_Expand_;

        AuxParms aux_parms_;

        std::vector<int> kSPPMap_;

    };
}

#endif //RHOMBUS_MATMUL_H

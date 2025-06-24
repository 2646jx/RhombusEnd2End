
#ifndef RHOMBUS_MATVEC_H
#define RHOMBUS_MATVEC_H

#include "status.h"
#include "statusor.h"
#include "matrix.h"
#include "seal_api.h"

namespace antchain::matvec
{
    using namespace antchain::util;

    /*!
     * Matrix-vector multiplication:
     * The matrix is in plaintext, vector is encrypted.
     */
    class RhombusMatVec
    {
    public:
        RhombusMatVec() = default;

        RhombusMatVec(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits);

        // create a RhombusMatVec object
        static StatusOr<std::unique_ptr<RhombusMatVec>> Create();

        // If you want to set the HE parameters independently to the global parameters, call this method instead.
        static StatusOr<std::unique_ptr<RhombusMatVec>> Create(uint32_t poly_degree, uint32_t mod_bits,
                                                               const std::vector<int> &coeff_mod_bits);

        // Generate keys: secret key, galois key and public key if in public key encryption mode
        // this function is only for test
        Status GenKey();

        // Generate serialized public key, return the size (bytes) of the key
        StatusOr<uint64_t> GenPublicKey(std::string &out) const;

        // Make sure that the buffer has enough space to save the pk
        StatusOr<uint64_t> GenPublicKey(uint8_t *buffer, size_t buffer_size) const;

        // Generate galois key
        StatusOr<uint64_t> GenGaloisKey(std::string &out) const;

        // Make sure that the buffer has enough space to save the gk
        StatusOr<uint64_t> GenGaloisKey(uint8_t *buffer, size_t buffer_size) const;

        StatusOr<uint64_t> GenGaloisKey(std::string &out, const std::vector<uint32_t> &galois_elts) const;

        StatusOr<uint64_t> GenGaloisKey(uint8_t *buffer, size_t buffer_size, const std::vector<uint32_t> &galois_elts) const;

        // Load public key from a string
        StatusOr<uint64_t> SetPublicKey(const std::string &in);

        StatusOr<uint64_t> SetPublicKey(const uint8_t *buffer, size_t buffer_size);

        // Load galois key from a string
        StatusOr<uint64_t> SetGaloisKey(const std::string &in);

        StatusOr<uint64_t> SetGaloisKey(const uint8_t *buffer, size_t buffer_size);

        [[nodiscard]] const seal::SecretKey & get_secret_key() const{
            return secret_key_;
        }

        // Note: you should re-generate the pk, gk once the secret key is reset.
        void reset_secret_key(const seal::SecretKey &new_sk);

        // (re)set the spp table
        void SetSPPMap(const std::vector<int> &new_spp_table){
            kSPPMap_ = new_spp_table;
        }

        const global::LargeMatInfo & MatInfo() const{
            return mat_info_;
        }

        const AuxParms & aux_parms() const {return aux_parms_;}

        uint32_t get_remain_mod_num() const {return remain_mod_num_;}

        void set_remain_mod_num(uint32_t n) {remain_mod_num_ = n;}

        /*!
         * Set matrix rows and columns. You must set rows and columns before performing MVM!!!
         * @param rows : rows of matrix
         * @param cols : columns of matrix
         * @param mat_bits : bit size of matrix element
         * @return
         */
        void SetMatrixRowsCols(uint32_t rows, uint32_t cols, uint32_t mat_bits = 20);

        /*!
         * Encode a vector in Z_{2^k}. For a number A in Z_{2^k}, we can represent it as binary complement
         * or straight binary formats. For example, k = 48, A = -5, then we can use -5 or 2^{48} - 5 to represent A.
         * This encoding method is to compute round(q * A / 2^k), and therefore, if A exceeds 2^k, it will be reduced
         * back to Z_{2^k}.
         * For vector v = (v_0, v_1, ..., v_{l-1}), we encode each entry to corresponding coefficient of the plain polynomial.
         * Concretely, Ecd(v) --> q/2^k * (v_0 - v_{l-1}x^{N-1} - v_{l-2}x^{N-2} - ... - v_1x^{N-l+1})
         * @tparam T : Type T can be int(int32_t), uint32_t, int64_t or uint64_t
         * @param plain : the destination
         * @param vec : the vector to be encoded
         * @param vec_len : the length of the vector
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status EncodeVec(seal::Plaintext &plain, const T *vec, size_t vec_len) const
        {

            antchain::encode_to_coeff(plain, vec, vec_len, Ecd_SCALED_INVERSE_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                      *seal_context_ptr, mod_bits_);
            return Status::OK;
        }

        /*!
         * EncodeVec: An extended version. The parameter enc_role decides the encoding strategies.
         * @tparam T : the elements type in vector
         * @param plain : output
         * @param vec : input vector
         * @param vec_len : length of the input vector
         * @param enc_role : encoding strategy
         * @param parms_id : parms_id, which decides the HE levels
         * @return status indicates the execution state.
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status EncodeVec(seal::Plaintext &plain, const T *vec, size_t vec_len, EcdRole enc_role,
                         seal::parms_id_type parms_id) const
        {
            antchain::encode_to_coeff(plain, vec, vec_len, enc_role, aux_parms_, parms_id, *seal_context_ptr, mod_bits_);
            return Status::OK;
        }

        /*!
         * Encode a vector of any length.
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status EncodeVec(std::vector<seal::Plaintext> &plain, const T *vec, size_t vec_len, EcdRole enc_role,
                         seal::parms_id_type parms_id) const
        {
            if (vec == nullptr)
                return {util::error::NULL_PTR, "nullptr"};
            size_t pt_num = (vec_len + poly_modulus_degree_ - 1) / poly_modulus_degree_;
            plain.resize(pt_num);
            const T *vec_ptr = vec;
            size_t remain_size = vec_len;

            for (size_t i = 0; i < pt_num; ++i){
                size_t cur_vec_len = std::min<size_t>(poly_modulus_degree_, remain_size);
                Status status = EncodeVec(plain[i], vec_ptr, cur_vec_len, enc_role, parms_id);
                if (!status.ok()) return status;
                vec_ptr += cur_vec_len;
                remain_size -= cur_vec_len;
            }
            return Status::OK;
        }

        /*!
         * Decode the plain polynomial
         * @tparam T : type T can be int, uint32_t, int64_t or uint64_t
         * @param vec : destination
         * @param vec_len : vector length
         * @param plain : plain polynomial, it should be in ntt form
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status DecodeVec(T *vec, size_t vec_len, const seal::Plaintext &plain) const
        {
            antchain::decode_from_coeff(vec, vec_len, plain, Dcd_SCALED_STRIDE, aux_parms_, *seal_context_ptr, mod_bits_);
            return Status::OK;
        }

        /*!
         * DecodeVec: Decode to a vector from the plaintext polynomial, the parameter dcd_role decides the decoding strategies
         * @tparam T : the elements type in the vector, which can be signed or unsigned integer of 1, 2, 4, 8 bytes
         * @param vec : output
         * @param vec_len : vector length
         * @param plain : input polynomial
         * @param dcd_role : decoding strategy
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status DecodeVec(T *vec, size_t vec_len, const seal::Plaintext &plain, DcdRole dcd_role) const{
            antchain::decode_from_coeff(vec, vec_len, plain, dcd_role, aux_parms_, *seal_context_ptr, mod_bits_);
            return Status::OK;
        }

        /*!
         * EncryptVec: encrypt a vector of length smaller than N
         * @tparam T : elements type of the vector
         * @param vec : input vector
         * @param vec_len : vector length
         * @param ctxt : output, save as a string
         * @param ecd_role : encoding strategy
         * @param use_sym_enc : use symmetric encryption
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status EncryptVec(const T *vec, size_t vec_len, std::string &ctxt, EcdRole ecd_role = Ecd_SCALED_INVERSE_ORDER, bool use_sym_enc = true) const
        {
            if (vec == nullptr || vec_len == 0 || vec_len > poly_modulus_degree_)
                return {util::error::INVALID_VECTOR_SIZE, "Invalid vector size"};

            seal::Plaintext plain;
            Status status = EncodeVec(plain, vec, vec_len, ecd_role, seal_context_ptr->first_parms_id());

            if (!status.ok())
                return {util::error::INVALID_ARGUMENT, "EncodeVec went wrong: " + status.ToString()};
            try
            {
                if (use_sym_enc)
                    encrypt(plain, secret_key_, plain.is_ntt_form(), ctxt, *seal_context_ptr);
                else
                    encrypt(plain, my_public_key_, plain.is_ntt_form(), ctxt, *seal_context_ptr);
            }
            CATCH_EXCEPTION()
            return Status::OK;
        }

        // encrypt a large vector
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        Status EncryptVec(const T *vec, size_t vec_len, std::vector<std::string> &ctxt,
                          EcdRole ecd_role = Ecd_SCALED_INVERSE_ORDER, bool use_sym_enc = true) const
        {
            if (vec == nullptr || vec_len == 0)
                return {util::error::INVALID_VECTOR_SIZE, "Invalid vector size"};

            const T *vec_ptr = vec;
            size_t remain_vec_size = vec_len;
            size_t ctxt_num = (vec_len + poly_modulus_degree_ - 1) / poly_modulus_degree_;
            ctxt.resize(ctxt_num);

            for (size_t i = 0; i < ctxt_num; ++i)
            {
                size_t cur_vec_size = std::min(remain_vec_size, size_t(poly_modulus_degree_));
                Status status = EncryptVec(vec_ptr, cur_vec_size, ctxt[i], ecd_role, use_sym_enc);
                if (!status.ok()) return status;
                vec_ptr += cur_vec_size;
                remain_vec_size -= cur_vec_size;
            }
            return Status::OK;
        }

        /*!
         * DecryptVec: decrypt a ciphertext, then decode it to a vector
         * @tparam T : vector type
         * @param enc_vec : input ciphertext
         * @param vec_result : output
         * @param vec_len : vector length
         * @param dcd_role : decoding strategy
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        Status DecryptVec(const seal::Ciphertext &enc_vec, T *vec_result, size_t vec_len, DcdRole dcd_role = Dcd_SCALED_STRIDE) const
        {
            if (vec_result == nullptr || vec_len == 0 || vec_len > poly_modulus_degree_)
                return {util::error::INVALID_VECTOR_SIZE, "Invalid vector size or nullptr"};

            seal::Plaintext plain;
            if (enc_vec.is_ntt_form())
                decryptor_->decrypt(enc_vec, plain);
            else
            {
                seal::Ciphertext ctxt(enc_vec);
                transform_to_ntt_inplace(ctxt, *seal_context_ptr);
                decryptor_->decrypt(ctxt, plain);
            }
            DecodeVec(vec_result, vec_len, plain, dcd_role);
            return Status::OK;
        }

        /// decrypt a large encrypted vector
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        Status DecryptVec(const std::vector<seal::Ciphertext> &enc_vec, T *vec_result, size_t vec_len) const
        {
            return DecryptVec(enc_vec, vec_result, vec_len, Dcd_SCALED_STRIDE);
        }

        /// decrypt a large encrypted vector
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        Status DecryptVec(const std::vector<seal::Ciphertext> &enc_vec, T *vec_result, size_t vec_len, DcdRole dcd_role) const
        {
            if (vec_result == nullptr || vec_len == 0)
                return {util::error::INVALID_VECTOR_SIZE, "Invalid vector size or nullptr"};

            size_t remain_vec_size = vec_len;
            size_t enc_vec_size = enc_vec.size();
            T *vec_ptr = vec_result;
            for (size_t i = 0; i < enc_vec_size; ++i)
            {
                size_t cur_vec_size = std::min(remain_vec_size, (size_t)poly_modulus_degree_);
                Status status = DecryptVec(enc_vec[i], vec_ptr, cur_vec_size, dcd_role);
                if (!status.ok()) return status;
                vec_ptr += cur_vec_size;
                remain_vec_size -= cur_vec_size;
            }
            return Status::OK;
        }

        /*!
         * The number of rows and columns < N
         * @tparam T : signed type to represent the matrix element, e.g. int8_t, int16_t, int32_t, int64_t
         * @param ciphertext : input encrypted vector
         * @param matrix : input plaintext matrix
         * @param result : output the encrypted result
         * @param multi_thread_count : the number of threads
         * @return status, indicate the correctness of the function
         */
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status MatVecMul(const seal::Ciphertext &ciphertext, const T *matrix,
                         seal::Ciphertext &result, uint64_t threads = 1) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            return MatVecMul(encoded_mat, ciphertext, result, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status MatVecMul(const seal::Ciphertext &ciphertext, const std::vector<std::vector<T>> &matrix,
                         seal::Ciphertext &result, uint32_t threads = 4) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            return MatVecMul(encoded_mat, ciphertext, result, threads);
        }


        /*!
         * Matrix-vector multiplication, convert the result to secret sharing form.
         * M * E(v) = E(mv_share0) + mv_share1
         * @tparam T : matrix elements type
         * @tparam U : output shared vector type
         * @param enc_vec : input encrypted vector
         * @param matrix : input plaintext matrix
         * @param mv_share0 : encrypted share0
         * @param mv_share1 : share1, in plaintext
         * @param threads : the number of threads
         * @return status, indicate the correctness of the function
         */
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        Status MatVecMulToSS(const seal::Ciphertext &enc_vec, const T *matrix,
                                  seal::Ciphertext &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            seal::Ciphertext enc_mv;
            Status state;
            state = MatVecMul(enc_vec, matrix, enc_mv, threads);
            if (!state.ok())
                return {util::error::INVALID_ARGUMENT, "MatVecMul went wrong: " + state.error_message()};
            return ConvToSS(enc_mv, mv_share0, mv_share1, mat_info_.corner_block_array_[0].nrows_, Dcd_SCALED_STRIDE);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        Status MatVecMulToSS(const seal::Ciphertext &enc_vec, const std::vector<std::vector<T>> &matrix,
                             seal::Ciphertext &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            seal::Ciphertext enc_mv;
            Status state = MatVecMul(enc_vec, matrix, enc_mv, threads);
            if (!state.ok())
                return {util::error::INVALID_ARGUMENT, state.error_message()};
            return ConvToSS(enc_mv, mv_share0, mv_share1, mat_info_.corner_block_array_[0].nrows_, Dcd_SCALED_STRIDE);
        }

        // For large matrix (the number of rows (or columns) > N)
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status LargeMatVecMul(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                   std::vector<seal::Ciphertext> &enc_mv, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                return {util::error::NULL_PTR, "Empty matrix is not allowed"};
            if (enc_vec.size() != mat_info_.ncol_block_num_)
                return {util::error::INCORRECT_CIPHERTEXT_NUMS, "matrix, ciphertexts number mismatch"};

            enc_mv.resize(mat_info_.nrow_block_num_);
            std::vector<seal::Ciphertext> temp_ctxt;
            std::vector<const T *> mat_block;

            size_t info_index;
            auto get_info = [&](size_t row_blk_idx, size_t col_blk_idx) -> size_t {
                if (row_blk_idx < mat_info_.nrow_block_num_ - 1 && col_blk_idx < mat_info_.ncol_block_num_ - 1)
                    return 0;
                else if ((row_blk_idx == mat_info_.nrow_block_num_ - 1) && (col_blk_idx < mat_info_.ncol_block_num_ - 1))
                    return 2;
                else if ((row_blk_idx < mat_info_.nrow_block_num_ - 1) && (col_blk_idx == mat_info_.ncol_block_num_ - 1))
                    return 1;
                else
                    return 3;
            };

            for (size_t i = 0; i < mat_info_.ncol_block_num_; ++i){
                for (size_t j = 0; j < mat_info_.nrow_block_num_; ++j){
                    GetBlockFromLargeMat(matrix, j, i, mat_block);
                    info_index = get_info(j, i);
                    if (i == 0)
                        MatVecMulInternal(enc_vec[i], mat_block, mat_info_.corner_block_array_[info_index],
                                               enc_mv[j], threads);
                    else
                    {
                        temp_ctxt.resize(mat_info_.nrow_block_num_);
                        MatVecMulInternal(enc_vec[i], mat_block, mat_info_.corner_block_array_[info_index],
                                               temp_ctxt[j], threads);
                        add_inplace(enc_mv[j], temp_ctxt[j], *seal_context_ptr);
                    }
                }
            }
            return Status::OK;
        }

        // Large matrix, convert the result to secret sharing
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        Status LargeMatVecMulToSS(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                       std::vector<seal::Ciphertext> &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> enc_mv;
            Status state;
            state = LargeMatVecMul(enc_vec, matrix, enc_mv, threads);
            if (!state.ok())
                return {util::error::INVALID_ARGUMENT, "LargeMatVecMul went wrong"};
            return ConvToSS(enc_mv, mv_share0, mv_share1, mat_info_.nrows_, Dcd_SCALED_STRIDE);
        }


        /*!
         * A special case: If the matrix with dimension r * c, and r >> N, c <= N, which is a tall matrix. Then this method
         * could be used to deal with this case. This method use the Expand function but the PackRLWEs.
         * @tparam T : the elements type in the matrix
         * @param ciphertext : the input encrypted vector of length c <= N
         * @param matrix : the input matrix
         * @param result : output, which is a encrypted vector of length r
         * @param multi_thread_count : the number of thread
         * @return
         */
         /// TODO: Add multi-thread for plaintext-ciphertext multiplication, add the split-point picking
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status LargeTallMatVecMul(const seal::Ciphertext &ciphertext, const T* &matrix,
                                   std::vector<seal::Ciphertext> &result, uint32_t multi_thread_count = 1) const
        {
            if (multi_thread_count < 1)
                return {util::error::INVALID_ARGUMENT, "threads count must >= 1"};
            std::vector<seal::Ciphertext> expanded_ct;

            // remove Expand factor
            uint32_t ell = BitLength(mat_info_.ncols_ - 1);
            seal::Ciphertext dup_ct = ciphertext;
            mul_inv_pow2_inplace(dup_ct, *seal_context_ptr, 1 << ell);

            // Expand function only accept ciphertext with INTT form
            if (dup_ct.is_ntt_form()){
                transform_from_ntt_inplace(dup_ct, *seal_context_ptr);
            }

            // Expand to 2^ell ciphertexts, then transform them to NTT form
            expand(dup_ct, expanded_ct, (int)ell, galois_keys_, *seal_context_ptr, (int)multi_thread_count);

            // transform the outputs (ciphertexts) of expand to ntt form
            transform_batched_to_ntt_inplace(expanded_ct, *seal_context_ptr, multi_thread_count);

            result.resize(mat_info_.nrow_block_num_);

            auto get_info = [&](size_t row_blk_idx, size_t col_blk_idx) -> size_t
            {
                if (row_blk_idx < mat_info_.nrow_block_num_ - 1 && col_blk_idx < mat_info_.ncol_block_num_ - 1)
                    return 0;
                else if ((row_blk_idx == mat_info_.nrow_block_num_ - 1) && (col_blk_idx < mat_info_.ncol_block_num_ - 1))
                    return 2;
                else if ((row_blk_idx < mat_info_.nrow_block_num_ - 1) && (col_blk_idx == mat_info_.ncol_block_num_ - 1))
                    return 1;
                else
                    return 3;
            };

            auto mul_pt = [&](size_t bgn, size_t end){
                std::vector<const T *> mat_block;
                std::vector<seal::Plaintext> encoded_mat;
                seal::Ciphertext temp_ct;
                for (size_t i = bgn; i < end; ++i){
                    GetBlockFromLargeMat(matrix, i, 0, mat_block);
                    size_t info_index = get_info(i, 0);
                    encode_matrix_column_major(mat_block, mat_info_.corner_block_array_[info_index], encoded_mat, 1);
                    set_zero_ct(result[i], *seal_context_ptr, ciphertext.parms_id());

                    for (size_t j = 0; j < mat_info_.ncols_; ++j){
                        multiply_plain_ntt(expanded_ct[j], encoded_mat[j], temp_ct, *seal_context_ptr);
                        add_inplace(result[i], temp_ct, *seal_context_ptr);
                    }
                }
            };

            uint32_t thread_block = (mat_info_.nrow_block_num_ + multi_thread_count - 1) / multi_thread_count;
            std::vector<std::thread> thread_pool(multi_thread_count);
            for (size_t i = 0; i < multi_thread_count; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, mat_info_.nrow_block_num_);
                thread_pool[i] = std::thread(mul_pt, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

            return Status::OK;
        }

        /// Convert the result to secret sharing form
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        Status LargeTallMatVecMulToSS(const seal::Ciphertext &enc_vec, const T *matrix,
                                  std::vector<seal::Ciphertext> &mv_share0, U *mv_share1, uint64_t multi_thread_count = 4) const
        {
            std::vector<seal::Ciphertext> enc_mv;
            Status state;
            state = LargeTallMatVecMul(enc_vec, matrix, enc_mv, multi_thread_count);
            if (!state.ok())
                return {util::error::INVALID_ARGUMENT, "LargeTallMatVecMul went wrong"};
            return ConvToSS(enc_mv, mv_share0, mv_share1, mat_info_.nrows_, Dcd_SCALED_IN_ORDER);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, int>::type * = nullptr>
        Status ConvToSS(const seal::Ciphertext &mv, seal::Ciphertext &mv_share0, T *mv_share1, size_t vec_len, DcdRole dcd_role) const
        {
            auto coeff_modulus = seal_context_ptr->get_context_data(mv.parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::Plaintext plain;
            plain.parms_id() = seal::parms_id_zero;
            plain.resize(poly_modulus_degree_ * coeff_modulus_size);

            std::mt19937_64 gen(std::random_device{}());
            for (size_t i = 0; i < coeff_modulus_size; ++i){
                std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[i].value() - 1);
                std::generate_n(plain.data() + i * poly_modulus_degree_, poly_modulus_degree_, [&](){
                    return dist(gen);});
            }
            plain.parms_id() = mv.parms_id();
            plain.scale() = 1.;

            DecodeVec(mv_share1, vec_len, plain, dcd_role);
            sub_plain(mv, plain, mv_share0, *seal_context_ptr);
            return Status::OK;
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, int>::type * = nullptr>
        Status ConvToSS(const std::vector<seal::Ciphertext> &mv, std::vector<seal::Ciphertext> &mv_share0,
                        T *mv_share1, size_t vec_len, DcdRole dcd_role) const
        {
            size_t remain_vec_len = vec_len;
            size_t ct_num = mv.size();
            mv_share0.resize(ct_num);

            T *vec_ptr = mv_share1;
            for (size_t i = 0; i < ct_num; ++i){
                size_t cur_vec_len = std::min<size_t>(remain_vec_len, poly_modulus_degree_);
                // To check
                ConvToSS(mv[i], mv_share0[i], vec_ptr, cur_vec_len, dcd_role);
                vec_ptr += cur_vec_len;
                remain_vec_len -= cur_vec_len;
            }
            return Status::OK;
        }

        // Serialization
        [[nodiscard]] StatusOr<uint64_t> GetCiphertextSize(const seal::Ciphertext &ciphertext) const;
        [[nodiscard]] StatusOr<uint64_t> GetCiphertextSize(const std::vector<seal::Ciphertext> &ciphertext_vec) const;
        StatusOr<uint64_t> CiphertextToBytes(const seal::Ciphertext &ciphertext, uint8_t *out, uint64_t size) const;
        StatusOr<uint64_t> CiphertextToBytes(const seal::Ciphertext &ciphertext, std::string &out) const;
        StatusOr<uint64_t> CiphertextToBytes(const std::vector<seal::Ciphertext> &ciphertext_vec, uint8_t *out, uint64_t size) const;
        StatusOr<uint64_t> CiphertextToBytes(const std::vector<seal::Ciphertext> &ciphertext_vec, std::vector<std::string> &out) const;
        [[nodiscard]] StatusOr<seal::Ciphertext> BytesToCiphertext(uint8_t *in, uint64_t size) const;
        [[nodiscard]] StatusOr<seal::Ciphertext> BytesToCiphertext(const std::string &in) const;
        [[nodiscard]] StatusOr<std::vector<seal::Ciphertext>> BytesToCiphertextVec(uint8_t *in, uint64_t size, uint64_t ctx_num) const;
        [[nodiscard]] StatusOr<std::vector<seal::Ciphertext>> BytesToCiphertextVec(const std::vector<std::string> &in) const;


        // Only used after performing MVP, the method will remove unused coefficients of c0
        // It could compress the size of the result ciphertext.
        Status drop_unused_coeffs(seal::Ciphertext &ctxt, size_t vec_size)
        {
            drop_unrelated_coeffs(ctxt, vec_size, *seal_context_ptr);
            return Status::OK;
        }

        [[nodiscard]] const seal::SEALContext & seal_context() const {
            return *seal_context_ptr;
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status EncodeMat(const T *mat, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (mat == nullptr){
                return {util::error::NULL_PTR, "null pointer to mat"};
            }
            global::MatInfo mat_info = mat_info_.corner_block_array_[0];
            std::vector<const T *> matrix(mat_info.nrows_);
            for (size_t i = 0; i < mat_info.nrows_; ++i){
                matrix[i] = mat + i * mat_info.ncols_;
            }

            return EncodeMatInternal(matrix, mat_info, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status EncodeMat(const std::vector<std::vector<T>> &mat,
                         std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (mat.empty())
                return {util::error::INVALID_MATRIX_ROWS_COLS, "empty matrix"};
            global::MatInfo mat_info = mat_info_.corner_block_array_[0];
            std::vector<const T *> matrix(mat_info.nrows_);
            for (size_t i = 0; i < mat_info.nrows_; ++i){
                matrix[i] = mat[i].data();
            }
            return EncodeMatInternal(matrix, mat_info, encoded_mat, threads);
        }


        Status MatVecMul(const std::vector<seal::Plaintext> &encoded_mat, const seal::Ciphertext &ct,
                       seal::Ciphertext &result, uint32_t threads = 4) const;


    private:

        // compute tau(ct) for all tau \in Gal(K_u/K_h)
        void hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                  const global::MatInfo &mat_info, uint32_t threads, bool remove_pack_factor = true) const
        {
            size_t g0_size = (size_t)1 << (mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_);
            cached_ct.resize(g0_size);

            cached_ct[0] = ct;
            if (remove_pack_factor)
                mul_inv_pow2_inplace(cached_ct[0], *seal_context_ptr, mat_info.spp_parms_.PackRLWEs_factor_);

            auto hom_aut = [&](uint32_t galois_elt, size_t src_index, size_t dest_index){
                cached_ct[dest_index] = cached_ct[src_index];
                apply_galois_inplace(cached_ct[dest_index], galois_elt, galois_keys_, *seal_context_ptr);
            };

            for (uint32_t i = 0, j = mat_info.log_pad_ncols_ - 1; i < (mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_); ++i,  --j){
                uint32_t thread_count = threads;
                uint32_t galois_elt = (poly_modulus_degree_ >> j) + 1;
                uint32_t total_step = (uint32_t)1 << i;
                for (uint32_t k = 0; k < total_step; k += thread_count){
                    size_t step_last = total_step - k;
                    thread_count = ((step_last < threads) ? step_last : thread_count);
                    std::vector<std::thread> thread_pool(thread_count);
                    for (uint32_t l = 0; l < thread_count; ++l)
                        thread_pool[l] = std::thread(hom_aut, galois_elt, k + l, total_step + k + l);
                    std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){ t.join(); });
                }
            }
        }


        // T should be int16_t, int32_t or int64_t
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        bool encode_to_coeff(seal::Plaintext &plain, const T *plain_vec, size_t vec_size) const
        {
            using namespace seal::util;
#ifdef RHOMBUS_DEBUG
            if (vec_size < 1 || vec_size > poly_modulus_degree_)
                throw std::invalid_argument("plain vector size is invalid");
            if (plain_vec == nullptr)
                throw std::invalid_argument("plain vector size is invalid");
#endif
            bool is_all_zero = std::all_of(plain_vec, plain_vec + vec_size,
                                           [](auto item)
                                           { return (item == 0); });

            auto parms_id = seal_context_ptr->first_parms_id();
            auto context_data = seal_context_ptr->get_context_data(parms_id);
            auto &parms = context_data->parms();
            auto &coeff_modulus = parms.coeff_modulus();
            auto coeff_modulus_size = parms.coeff_modulus().size();
            auto ntt_tables = context_data->small_ntt_tables();

            plain.parms_id() = seal::parms_id_zero;
            size_t buffer_size = poly_modulus_degree_ * coeff_modulus_size;
            plain.resize(buffer_size);
            std::fill_n(plain.data(), buffer_size, 0);

            if (is_all_zero){
                plain.parms_id() = parms_id;
                return true;
            }

            for (size_t i = 0; i < coeff_modulus_size; ++i)
            {
                auto offset = i * poly_modulus_degree_;
                std::transform(plain_vec, plain_vec + vec_size, plain.data() + offset, [&](auto elt)
                               {
                    if (elt >= 0)
                        return seal::util::barrett_reduce_64((uint64_t) elt, coeff_modulus[i]);
                    else {
                        uint64_t v = seal::util::barrett_reduce_64(static_cast<uint64_t>(elt), coeff_modulus[i]);
                        return seal::util::sub_uint_mod(v, aux_parms_.pow_2_64_mod_qi_[i], coeff_modulus[i]);
                    } });
                // lazy ntt
                ntt_negacyclic_harvey_lazy(plain.data(offset), ntt_tables[i]);
            }
            plain.parms_id() = parms_id;
            plain.scale() = 1.;
            return false;
        }

        // encode an uint128_t vector with length N
        bool encode_to_coeff128(seal::Plaintext &plain, const uint64_t *plain_vec, size_t vec_size) const;

        // get sub-matrix at (block_row_index, block_col_index)
        template <typename T, typename std::enable_if<std::is_arithmetic<T>::value, T>::type * = nullptr>
        void GetBlockFromLargeMat(const T *large_mat, uint32_t block_row_index,
                                  uint32_t block_col_index, std::vector<const T *> &destination) const
        {
            if (block_row_index >= mat_info_.nrow_block_num_ || block_col_index >= mat_info_.ncol_block_num_)
                throw std::invalid_argument("block index should be in [0, block_row_num), [0, block_col_num)");
            uint32_t block_rows = ((((block_row_index + 1) << global::kLogPolyDegree) <= mat_info_.nrows_)
                    ? poly_modulus_degree_ : (mat_info_.nrows_ - (block_row_index << global::kLogPolyDegree)));

            destination.resize(block_rows);
            for (size_t i = 0; i < block_rows; ++i)
            {
                size_t row_offset = (block_row_index * poly_modulus_degree_) + i;
                size_t col_offset = block_col_index * poly_modulus_degree_;
                destination[i] = large_mat + row_offset * mat_info_.ncols_ + col_offset;
            }
        }


        /*!
         * For expand-based method. Encode the matrix by column, each column will be encoded into a polynomial
         * @tparam T
         * @param matrix
         * @param mat_info
         * @param encoded_mat
         * @param threads
         * @return
         */
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        Status encode_matrix_column_major(const std::vector<const T *> &matrix, const global::MatInfo &mat_info,
                                      std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1)
                return {util::error::INVALID_ARGUMENT, "threads number must >= 1"};
            size_t nrows = mat_info.nrows_;
            size_t ncols = mat_info.ncols_;

            if (nrows < 1 || nrows > poly_modulus_degree_)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "the number of the rows must be in [1, N]"};
            if (ncols < 1 || ncols > poly_modulus_degree_)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "the number of the columns must be in [1, N]"};

            encoded_mat.resize(ncols);

            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(nrows);
                for (size_t j = bgn; j < end; ++j){
                    std::fill_n(temp.data(), nrows, 0);
                    for (size_t i = 0; i < nrows; ++i)
                        temp[i] = matrix[i][j];
                    encode_to_coeff(encoded_mat[j], temp.data(), nrows);
                }
            };

            uint32_t thread_block = (ncols + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, ncols);
                thread_pool[i] = std::thread(ecd_program, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

            return Status::OK;
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status EncodeMatInternal(const std::vector<const T *> &matrix, const global::MatInfo &mat_info,
                                 std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            // check matrix
            if (mat_info.nrows_ > poly_modulus_degree_ || mat_info.ncols_ > poly_modulus_degree_)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "matrix rows or columns out of bound"};
            if (mat_info.nrows_ == 0 || mat_info.ncols_ == 0)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "the number of rows, columns should not be zero"};

            // input packing
            std::vector<std::vector<T>> con_mat;
            global::matrix_concat(matrix, poly_modulus_degree_, con_mat, mat_info);

            size_t g0_size = (size_t)1 << (mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_);
            size_t g1_size = (size_t)1 << (mat_info.log_pad_nrows_ - mat_info.spp_parms_.u_);

            encoded_mat.resize(g0_size * g1_size);

            if (mat_info.mat_bits_ + mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_ <= sizeof(T) * 8){
                std::vector<std::vector<std::vector<T>>> pp_mat;
                global::matrix_row_combine64(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            antchain::encode_to_coeff(encoded_mat[i1 * g0_size + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                      Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                                      *seal_context_ptr, mod_bits_);
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
                std::for_each(thread_pool.begin(),
                              thread_pool.end(), [](std::thread &t){t.join();});
            }
            else if (mat_info.mat_bits_ + mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_ <= 64){
                std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                global::matrix_row_combine64(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            antchain::encode_to_coeff(encoded_mat[i1 * g0_size + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                      Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                                      *seal_context_ptr, mod_bits_);
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
                std::for_each(thread_pool.begin(),
                              thread_pool.end(), [](std::thread &t){t.join();});
            }
            else {
                std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                global::matrix_row_combine128(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                              poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            antchain::encode_to_coeff128(encoded_mat[i1 * g0_size + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                         aux_parms_, seal_context_ptr->first_parms_id(), *seal_context_ptr);
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
                std::for_each(thread_pool.begin(),
                              thread_pool.end(), [](std::thread &t){t.join();});
            }
            return Status::OK;
        }

        Status MatVecMulInternal(const seal::Ciphertext &ct, const std::vector<seal::Plaintext> &encoded_mat,
                                      const global::MatInfo &mat_info, seal::Ciphertext &result, uint32_t threads = 4) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        Status MatVecMulInternal(const seal::Ciphertext &ciphertext, const std::vector<const T *> &matrix, const global::MatInfo &mat_info,
                                      seal::Ciphertext &result, uint32_t threads = 4) const
        {
            // check matrix
            if (mat_info.nrows_ > poly_modulus_degree_ || mat_info.ncols_ > poly_modulus_degree_)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "matrix rows or columns out of bound"};
            if (mat_info.nrows_ == 0 || mat_info.ncols_ == 0)
                return {util::error::INVALID_MATRIX_ROWS_COLS, "the number of rows, columns should not be zero"};

            // input packing: concatenate some rows of the matrix into some length-N vectors
            std::vector<std::vector<T>> con_mat;
            global::matrix_concat(matrix, poly_modulus_degree_, con_mat, mat_info);

            // group size of G0, G1
            size_t g0_size = (size_t)1 << (mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_);
            size_t g1_size = (size_t)1 << (mat_info.log_pad_nrows_ - mat_info.spp_parms_.u_);

            // MergeAut: compute {tau(ciphertext)} for tau in G0
            std::vector<seal::Ciphertext> cached_ct(g0_size);
            hom_aut_galois_group(ciphertext, cached_ct, mat_info, threads, true);

            // save the inner products
            std::vector<seal::Ciphertext> matvec(g1_size);

            // Here, we divide into three cases:
            // 1. the size of the elements of the preprocessed matrix <= bit_count(T)
            // 2. the size of the elements of the preprocessed matrix <= 64
            // 3. the size of the elements of the preprocessed matrix <= 128
            if (mat_info.mat_bits_ + mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_ <= sizeof(T) * 8){

                // matrix-preprocessing, corresponds to the SPP optimization
                // pp_mat: g1_size * g0_size, each element is a length-N vector
                std::vector<std::vector<std::vector<T>>> pp_mat;
                global::matrix_row_combine64(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                // compute the inner products
                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            encode_to_coeff(pt, pp_mat[i1][j0].data(), poly_modulus_degree_);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > 2){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                uint32_t thread_block = (g1_size + threads - 1) / threads;
                std::vector<std::thread> thread_pool(threads);
                for (size_t i = 0; i < threads; ++i){
                    size_t bgn = i * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, g1_size);
                    thread_pool[i] = std::thread(hom_inner_prod, bgn, end);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }
            else if (mat_info.mat_bits_ + mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_ <= 64){
                std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                global::matrix_row_combine64(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            encode_to_coeff(pt, pp_mat[i1][j0].data(), poly_modulus_degree_);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > 2){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                uint32_t thread_block = (g1_size + threads - 1) / threads;
                std::vector<std::thread> thread_pool(threads);
                for (size_t i = 0; i < threads; ++i){
                    size_t bgn = i * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, g1_size);
                    thread_pool[i] = std::thread(hom_inner_prod, bgn, end);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }
            else{
                std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                global::matrix_row_combine128(con_mat, mat_info.spp_parms_, pp_mat, threads,
                                              poly_modulus_degree_, poly_modulus_degree_);

                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < g0_size; ++j0){
                            encode_to_coeff128(pt, pp_mat[i1][j0].data(), poly_modulus_degree_);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > 2){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                uint32_t thread_block = (g1_size + threads - 1) / threads;
                std::vector<std::thread> thread_pool(threads);
                for (size_t i = 0; i < threads; ++i){
                    size_t bgn = i * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, g1_size);
                    thread_pool[i] = std::thread(hom_inner_prod, bgn, end);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }

            // PackRLWEs.
            PackRLWEs(matvec, mat_info.spp_parms_.u_, galois_keys_, result, *seal_context_ptr, threads);
            transform_to_ntt_inplace(result, *seal_context_ptr);
            return Status::OK;
        }

        uint32_t poly_modulus_degree_;
        uint32_t mod_bits_;

        uint32_t remain_mod_num_ = 2;

//        bool use_PackRLWEs_based_method_ = true;

        std::shared_ptr<seal::SEALContext> seal_context_ptr;
        std::unique_ptr<seal::KeyGenerator> keygen_;
        std::unique_ptr<seal::Decryptor> decryptor_;

        seal::SecretKey secret_key_;
        seal::PublicKey my_public_key_;

        // save other party's pk, gk.
        seal::PublicKey public_key_;
        seal::GaloisKeys galois_keys_;

        global::LargeMatInfo mat_info_;
        AuxParms aux_parms_;

        std::vector<int> kSPPMap_;
    };
} // namespace

#endif // RHOMBUS_MATVEC_H

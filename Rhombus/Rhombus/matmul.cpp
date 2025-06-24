
#include "matmul.h"
#include <memory>

namespace antchain::matmul{
    using namespace std;

    RhombusMatMul::RhombusMatMul() {

        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(kPolyDegree);

        switch (kNumModulus) {
            case 2:
                parms.set_coeff_modulus({kModulus[0], kModulus[1]});
                break;
            case 3:
                parms.set_coeff_modulus({kModulus[0], kModulus[1], kModulus[2]});
                break;
            case 4:
                parms.set_coeff_modulus({kModulus[0], kModulus[1], kModulus[2], kModulus[3]});
                break;
            default:
                throw std::invalid_argument("Modulus number must be 2, 3 or 4");
        }

        seal_context_ = std::make_shared<seal::SEALContext>(parms);

        // copy modulus
        std::vector<uint64_t> Allmods(kNumModulus);
        for (size_t i = 0; i < kNumModulus; ++i)
            Allmods[i] = kModulus[i];

        gen_aux_params(aux_parms_, kSSBitLen, Allmods, kNumModulus);

        keygen_ = make_unique<seal::KeyGenerator>(*seal_context_);
        secret_key_ = keygen_->secret_key();
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_, secret_key_);

        keygen_->create_public_key(my_public_key_);
        poly_modulus_degree_ = kPolyDegree;
        mod_bits_ = global::kSSBitLen;

        // use the global spp table
        kSPPMap_.resize(global::kLogPolyDegree + 1);
        std::copy_n(global::kSPPMap, global::kLogPolyDegree + 1, kSPPMap_.begin());
    }

    RhombusMatMul::RhombusMatMul(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits) {
        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_degree);
        size_t mod_num = coeff_mod_bits.size();
        if (mod_num < 2 || mod_num > 4)
            throw std::invalid_argument("Modulus number must be 2, 3 or 4");
        parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_degree, coeff_mod_bits));
        seal_context_ = std::make_shared<seal::SEALContext>(parms);

        // copy modulus
        std::vector<uint64_t> Allmods(mod_num);
        for (size_t i = 0; i < mod_num; ++i)
            Allmods[i] = parms.coeff_modulus()[i].value();

        gen_aux_params(aux_parms_, mod_bits, Allmods, mod_num);

        keygen_ = std::make_unique<seal::KeyGenerator>(*seal_context_);
        secret_key_ = keygen_->secret_key();
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_, secret_key_);

        // usually unused
        keygen_->create_public_key(my_public_key_);
        poly_modulus_degree_ = poly_degree;
        mod_bits_ = mod_bits;

        // we fix the SPP table as follows. You can also reset it by calling the SetSPPMap method.
        kSPPMap_ = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5};
    }

    void RhombusMatMul::reset_mod_bits(uint32_t mod_bits) {
        auto &parms = seal_context_->key_context_data()->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_mod_size = coeff_modulus.size();

        // copy modulus
        std::vector<uint64_t> Allmods(coeff_mod_size);
        for (size_t i = 0; i < coeff_mod_size; ++i){
            Allmods[i] = coeff_modulus[i].value();
        }
        gen_aux_params(aux_parms_, mod_bits, Allmods, coeff_mod_size);
        mod_bits_ = mod_bits;
    }

    void RhombusMatMul::GenKey() {
        size_t galois_elt_count = seal::util::get_power_of_two(poly_modulus_degree_);
        vector<uint32_t> galois_elt(galois_elt_count);

        // default: 13 galois elements {3,5,9,17,33,65,129,257,513,1025,2049,4097,8193} if N = 8192
        for (size_t i = 0; i < galois_elt_count; i++)
            galois_elt[i] = (1UL << (i + 1)) + 1;

        keygen_->create_galois_keys(galois_elt, galois_keys_);
    }

    void RhombusMatMul::SetMatDims(uint32_t n, uint32_t m, uint32_t k, uint32_t X_bits, uint32_t Y_bits) {
        gen_lw_mat_info(mat_X_, n, m, X_bits);
        gen_lw_mat_info(mat_Y_, m, k, Y_bits);

        // ----------------------- PackRLWEs based SPP parameters ---------------------------

        // block num of X: ceil(n/pad_m) = (n + pad_m - 1) / pad_m
        uint32_t pad_m = uint32_t(1) << mat_X_.log_pad_ncols_;
        uint32_t bnX = (mat_X_.nrows_ + pad_m - 1) >> mat_X_.log_pad_ncols_;
        uint32_t logN = BitLength(poly_modulus_degree_ - 1);

        // for the last block
        uint32_t lbrow = mat_X_.nrows_ - ((bnX - 1)  << mat_X_.log_pad_ncols_);
        uint32_t logpad_lrrow = BitLength(lbrow - 1);

        gen_spp_parms(spp_PackRLWEs_[1], logN - (mat_X_.log_pad_ncols_ - logpad_lrrow),
                      mat_X_.log_pad_ncols_, poly_modulus_degree_, kSPPMap_.data());

        if (bnX == 1)
            spp_PackRLWEs_[0] = spp_PackRLWEs_[1];
        else{
            // generate first spp parameters. If block_num = 1, then this spp won't be used.
            gen_spp_parms(spp_PackRLWEs_[0], logN, mat_X_.log_pad_ncols_,
                          poly_modulus_degree_, kSPPMap_.data());
        }

        // ----------------------- Expand based SPP parameters ---------------------------

        gen_spp_expand(spp_Expand_, BitLength((poly_modulus_degree_ >> mat_Y_.log_pad_ncols_) - 1) /* N/pad_k */,
                       poly_modulus_degree_, kSPPMap_.data());

    }

    void RhombusMatMul::reset_secret_key(const seal::SecretKey &new_sk) {
        keygen_ = make_unique<seal::KeyGenerator>(*seal_context_, new_sk);
        secret_key_ = keygen_->secret_key();
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_, secret_key_);
        keygen_->create_public_key(my_public_key_);
    }

    size_t RhombusMatMul::GenPublicKey(std::string &out) const {
        seal::Serializable<seal::PublicKey> pk = keygen_->create_public_key();

        std::ostringstream ostr;
        size_t save_size = pk.save(ostr);
        ostr.str().swap(out);
        return save_size;
    }

    size_t RhombusMatMul::GenPublicKey(uint8_t *buffer, size_t buffer_size) const {
        seal::Serializable<seal::PublicKey> pk = keygen_->create_public_key();

        size_t save_size = pk.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer is not enough to save the pk");
        return pk.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t RhombusMatMul::GenGaloisKey(std::string &out) const {
        size_t galois_elt_count = BitLength(poly_modulus_degree_ - 1);
        vector<uint32_t> galois_elt(galois_elt_count);

        // default: 13 galois elements {3,5,9,17,33,65,129,257,513,1025,2049,4097,8193} if N = 8192
        for (size_t i = 0; i < galois_elt_count; i++)
            galois_elt[i] = (1UL << (i + 1)) + 1;
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elt);
        std::ostringstream ostr;
        size_t save_size = gk.save(ostr);
        ostr.str().swap(out);
        return save_size;
    }

    size_t RhombusMatMul::GenGaloisKey(uint8_t *buffer, size_t buffer_size) const {
        size_t galois_elt_count = BitLength(poly_modulus_degree_ - 1);
        vector<uint32_t> galois_elt(galois_elt_count);

        for (size_t i = 0; i < galois_elt_count; ++i)
            galois_elt[i] = (1UL << (i + 1)) + 1;

        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elt);
        size_t save_size = gk.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer is not enough to save the gk");
        return gk.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t RhombusMatMul::GenGaloisKey(std::string &out, const std::vector<uint32_t> &galois_elts) const {
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elts);
        std::ostringstream ostr;
        size_t save_size = gk.save(ostr);
        ostr.str().swap(out);
        return save_size;
    }

    size_t RhombusMatMul::GenGaloisKey(uint8_t *buffer, size_t buffer_size,
                                       const std::vector<uint32_t> &galois_elts) const {
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elts);
        size_t save_size = gk.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer is not enough to save gk");
        return gk.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t RhombusMatMul::SetPublicKey(const std::string &in) {
        if (in.empty())
            throw invalid_argument("empty string");
        size_t load_size = public_key_.load(*seal_context_, (seal::seal_byte *)in.data(), in.size());
        return load_size;
    }

    size_t RhombusMatMul::SetPublicKey(const uint8_t *buffer, size_t buffer_size) {
        if (buffer == nullptr)
            throw std::invalid_argument("invalid buffer, null pointer");
        size_t load_size = public_key_.load(*seal_context_, (seal::seal_byte *)buffer, buffer_size);
        return load_size;
    }

    size_t RhombusMatMul::SetGaloisKey(const std::string &in) {
        if (in.empty())
            throw invalid_argument("empty string");
        size_t load_size = galois_keys_.load(*seal_context_, (seal::seal_byte *)in.data(), in.size());
        return load_size;
    }

    size_t RhombusMatMul::SetGaloisKey(const uint8_t *buffer, size_t buffer_size) {
        if (buffer == nullptr)
            throw std::invalid_argument("invalid buffer, null pointer");
        size_t load_size = galois_keys_.load(*seal_context_, (seal::seal_byte *)buffer, buffer_size);
        return load_size;
    }

    void RhombusMatMul::MatMulBlock(const seal::Plaintext *ecd_submatX,
                                    const std::vector<seal::Ciphertext> &cached_subY, seal::Ciphertext &enc_submatXY,
                                    const antchain::global::SPPParams &spp, uint32_t threads, bool mul_factor) const
    {
        size_t g1_size = size_t(1) << (spp.ell_ + spp.h_ - spp.u_);
        size_t g0_size = size_t(1) << (spp.u_ - spp.h_);

        std::vector<seal::Ciphertext> inn_prod(g1_size);
        seal::Ciphertext ct;

        // compute the matrix multiplication
        auto hom_compute = [&](size_t bgn, size_t end){
            seal::Ciphertext ct;
            for (size_t i1 = bgn; i1 < end; ++i1){
                set_zero_ct(inn_prod[i1], *seal_context_, seal_context_->first_parms_id(), true);
                for (size_t j = 0; j < g0_size; ++j){
                    multiply_plain_ntt(cached_subY[j], ecd_submatX[i1 * g0_size + j], ct, *seal_context_);
                    add_inplace(inn_prod[i1], ct, *seal_context_);
                }
                transform_from_ntt_inplace(inn_prod[i1], *seal_context_);
                
                while (inn_prod[i1].coeff_modulus_size() > remain_mod_num_)
                    rescale_to_next_inplace(inn_prod[i1], *seal_context_);
            }
        };

        // multi-thread
        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(hom_compute, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

        // pack the sparse encoded ciphertexts
        PackRLWEs(inn_prod, spp.u_, galois_keys_, enc_submatXY, *seal_context_, threads);
        if (mul_factor)
            mul_uint_inplace(enc_submatXY, spp_PackRLWEs_[0].PackRLWEs_factor_ / spp_PackRLWEs_[1].PackRLWEs_factor_, *seal_context_);
        transform_to_ntt_inplace(enc_submatXY, *seal_context_);
    }

    size_t RhombusMatMul::CiphertextsToBytes(const std::vector<std::vector<seal::Ciphertext>> &ct,
                                             std::vector<std::vector<std::string>> &out) const {
        size_t save_size = 0;
        size_t dim0 = ct.size();
        if (dim0 == 0)
            throw std::invalid_argument("empty ciphertext");
        size_t dim1 = ct[0].size();
        if (dim1 == 0)
            throw std::invalid_argument("empty ciphertext");

        out.resize(dim0);
        for (size_t i = 0; i < dim0; ++i){
            out[i].resize(dim1);
            for (size_t j = 0; j < dim1; ++j){
                std::ostringstream ostr;
                size_t cur_size = ct[i][j].save(ostr);
                save_size += cur_size;
                ostr.str().swap(out[i][j]);
            }
        }
        return save_size;
    }

    size_t RhombusMatMul::CiphertextsToBytes(const std::vector<seal::Ciphertext> &ct,
                                             std::vector<std::string> &out) const {
        size_t save_size = 0;
        if (ct.empty())
            throw std::invalid_argument("empty ct vector");
        size_t ct_num = ct.size();
        out.resize(ct_num);
        for (size_t i = 0; i < ct_num; ++i){
            std::ostringstream ostr;
            size_t cur_size = ct[i].save(ostr);
            save_size += cur_size;
            ostr.str().swap(out[i]);
        }
        return save_size;
    }

    size_t RhombusMatMul::BytesToCiphertexts(const std::vector<std::string> &ct_str,
                                             std::vector<seal::Ciphertext> &ct_out) const {
        size_t load_size = 0;
        size_t dim = ct_str.size();
        if (dim == 0)
            throw std::invalid_argument("empty ct string");

        ct_out.resize(dim);
        for (size_t i = 0; i < dim; ++i){
            size_t cur_size = ct_out[i].unsafe_load(*seal_context_, (seal::seal_byte *)ct_str[i].data(), ct_str[i].size());
            load_size += cur_size;
        }
        return load_size;
    }

    size_t RhombusMatMul::BytesToCiphertexts(const std::vector<std::vector<std::string>> &ct_str,
                                             std::vector<std::vector<seal::Ciphertext>> &ct_out) const {
        size_t load_size = 0;
        size_t dim0 = ct_str.size();
        if (dim0 == 0)
            throw std::invalid_argument("empty ct string");
        size_t dim1 = ct_str[0].size();
        if (dim1 == 0)
            throw std::invalid_argument("empty ct string");

        ct_out.resize(dim0);
        for (size_t i = 0; i < dim0; ++i){
            ct_out[i].resize(dim1);
            for (size_t j = 0; j < dim1; ++j){
                size_t cur_size = ct_out[i][j].unsafe_load(*seal_context_, (seal::seal_byte *)ct_str[i][j].data(), ct_str[i][j].size());
                load_size += cur_size;
            }
        }
        return load_size;
    }

    void RhombusMatMul::hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                             const global::SPPParams &spp, uint32_t threads,
                                             bool remove_pack_factor) const
    {
        size_t g0_size = (size_t)1 << (spp.u_ - spp.h_);
        cached_ct.resize(g0_size);

        cached_ct[0] = ct;
        if (remove_pack_factor)
            mul_inv_pow2_inplace(cached_ct[0], *seal_context_, spp.PackRLWEs_factor_);

        auto hom_aut = [&](uint32_t galois_elt, size_t src_index, size_t dest_index){
            cached_ct[dest_index] = cached_ct[src_index];
            apply_galois_inplace(cached_ct[dest_index], galois_elt, galois_keys_, *seal_context_);
        };

        uint32_t logN = BitLength(poly_modulus_degree_ - 1);
        for (uint32_t i = 0, j = logN - spp.h_ - 1; i < (spp.u_ - spp.h_); ++i,  --j){
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

    void RhombusMatMul::hom_aut_and_add_group(std::vector<seal::Ciphertext> &cts, seal::Ciphertext &result,
                                              const global::SPPParams_Expand &spp, uint32_t threads) const {
        size_t g0_size = (1ULL << (spp.u_ + spp.ell_)) / poly_modulus_degree_;
        if (g0_size != cts.size())
            throw std::invalid_argument("ct size and g0_size are mismatched");
        if (g0_size == 1){
            result = cts[0];
            return;
        }

        uint32_t logN = BitLength(poly_modulus_degree_ - 1);
        uint32_t cur_exp = logN - spp.ell_ + 1; // current
        uint32_t cur_galois_elt = ((uint32_t)1 << cur_exp) + 1;
        size_t cur_ct_size = g0_size;
        seal::Ciphertext temp;

        while (cur_ct_size > 1){
            for (size_t i = 0; i < (cur_ct_size >> 1); ++i){
                temp = cts[2 * i + 1];
                apply_galois_inplace(temp, cur_galois_elt, galois_keys_, *seal_context_);
                add(cts[2 * i], temp, cts[i], *seal_context_);
            }
            cur_ct_size >>= 1;
            ++cur_exp;
            cur_galois_elt = ((uint32_t)1 << cur_exp) + 1;
        }

        result = cts[0];
    }

    void RhombusMatMul::EncryptMatYInternal(const std::vector<seal::Plaintext> &encoded_mat,
                                              seal::Ciphertext *encrypted_mat, std::string *serialized_enc_mat,
                                              uint32_t threads, bool is_symmetric) const {
        uint32_t ct_num = encoded_mat.size();
        if ((encrypted_mat == nullptr) && (serialized_enc_mat == nullptr))
            throw std::invalid_argument("both destinations is nullptr");

        if ((encrypted_mat != nullptr) && (serialized_enc_mat != nullptr))
            std::cout << "Encrypt twice, slow case" << std::endl;

        auto enc_slice = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                if (is_symmetric){
                    if (encrypted_mat != nullptr){
                        encrypt(encoded_mat[i], secret_key_, true, encrypted_mat[i], *seal_context_);
                    }
                    if (serialized_enc_mat != nullptr){
                        encrypt(encoded_mat[i], secret_key_, true, serialized_enc_mat[i], *seal_context_);
                    }
                }else{
                    if (encrypted_mat != nullptr){
                        encrypt(encoded_mat[i], my_public_key_, true, encrypted_mat[i], *seal_context_);
                    }
                    if (serialized_enc_mat != nullptr){
                        encrypt(encoded_mat[i], my_public_key_, true, serialized_enc_mat[i], *seal_context_);
                    }
                }
            }
        };

        uint32_t thread_block = (ct_num + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, ct_num);
            thread_pool[i] = std::thread(enc_slice, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusMatMul::EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat,
                                      std::vector<seal::Ciphertext> &encrypted_mat, uint32_t threads,
                                      bool is_symmetric) const {
        size_t ct_num = encoded_mat.size();
        encrypted_mat.resize(ct_num);
        EncryptMatYInternal(encoded_mat, encrypted_mat.data(), nullptr, threads, is_symmetric);
    }

    void RhombusMatMul::EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat,
                                      std::vector<std::string> &encrypted_mat, uint32_t threads,
                                      bool is_symmetric) const {
        size_t ct_num = encoded_mat.size();
        encrypted_mat.resize(ct_num);
        EncryptMatYInternal(encoded_mat, nullptr, encrypted_mat.data(), threads, is_symmetric);
    }

    void RhombusMatMul::hom_inner_product_E(std::vector<seal::Ciphertext> &temp, const seal::Plaintext *encoded_matX,
                                            const std::vector<std::vector<seal::Ciphertext>> &cached_ct,
                                            uint32_t threads) const {
        uint32_t logN = BitLength(poly_modulus_degree_ - 1);
        uint32_t g0_size = 1UL << (spp_Expand_.u_ + spp_Expand_.ell_ - logN);
        uint32_t g1_size = 1UL << (logN - spp_Expand_.u_);
        uint32_t block_size = poly_modulus_degree_ >> mat_Y_.log_pad_ncols_;

        // the number of column blocks of X
        uint32_t nblock = (mat_X_.ncols_ + block_size - 1) / block_size;

        auto inn_prg = [&](size_t bgn, size_t end){
            seal::Ciphertext temp_ct;
            for (size_t k1 = bgn; k1 < end; ++k1){
                set_zero_ct(temp[k1], *seal_context_, seal_context_->first_parms_id());
                for (size_t k0 = 0; k0 < nblock; ++k0){
                    for (size_t k2 = 0; k2 < g1_size; ++k2){
                        multiply_plain_ntt(cached_ct[k0][k2], encoded_matX[k0 * g1_size * g0_size + k2 * g0_size + k1], temp_ct, *seal_context_);
                        add_inplace(temp[k1], temp_ct, *seal_context_);
                    }
                }
                while (temp[k1].coeff_modulus_size() > remain_mod_num_){
                    rescale_to_next_inplace(temp[k1], *seal_context_);
                }
            }
        };

        uint32_t thread_block = (g0_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g0_size);
            thread_pool[i] = std::thread(inn_prg, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusMatMul::Compute_E(const std::vector<seal::Plaintext> &encoded_matX,
                                  const std::vector<seal::Ciphertext> &encrypted_matY,
                                  std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        // N / pad_k, which means:
        // 1. 'block_size' rows of Y will be packed into one polynomial
        // 2. X will be divided into several 'block_size * block_size' sub-matrices.
        uint32_t block_size = poly_modulus_degree_ >> mat_Y_.log_pad_ncols_;
        // the number of column blocks of X
        uint32_t nblock = (mat_X_.ncols_ + block_size - 1) / block_size;
        uint32_t logN = BitLength(poly_modulus_degree_ - 1);

        uint32_t g0_size = 1UL << (spp_Expand_.u_ + spp_Expand_.ell_ - logN);
        uint32_t g1_size = 1UL << (logN - spp_Expand_.u_);
        uint32_t submatX_num = (mat_X_.nrows_ + block_size - 1) / block_size;

        // Expand Enc(Y)
        std::vector<std::vector<seal::Ciphertext>> cached_ct(nblock);

        auto expand_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                Expand(encrypted_matY[i], cached_ct[i], logN - spp_Expand_.u_, spp_Expand_.Expand_factor_,
                       galois_keys_, *seal_context_, 1);
            }
        };

        uint32_t thread_block = (nblock + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, nblock);
            thread_pool[i] = std::thread(expand_program, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

        encrypted_matXY.resize(submatX_num);
        uint32_t stride = nblock * g1_size * g0_size;

        auto hom_mul_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                set_zero_ct(encrypted_matXY[i], *seal_context_, seal_context_->first_parms_id());
                std::vector<seal::Ciphertext> temp(g0_size);
                hom_inner_product_E(temp, encoded_matX.data() + i * stride, cached_ct, 1);
                hom_aut_and_add_group(temp, encrypted_matXY[i], spp_Expand_, 1);
            }
        };

        thread_block = (submatX_num + threads - 1) / threads;
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, submatX_num);
            thread_pool[i] = std::thread(hom_mul_program, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusMatMul::Compute_P(const std::vector<seal::Plaintext> &encoded_matX,
                                  const std::vector<seal::Ciphertext> &encrypted_matY,
                                  std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        if (threads == 0 || threads > THREAD_NUM_MAX)
            throw std::invalid_argument("invalid thread number");

        // partition windows
        size_t mw = size_t(1) << mat_X_.log_pad_ncols_;
        size_t kw = poly_modulus_degree_ / mw;

        // ceil(n / mw), ceil(k / kw)
        size_t submatX_num = (mat_X_.nrows_ + mw - 1) / mw;
        size_t submatY_num = (mat_Y_.ncols_ + kw - 1) / kw;

        if (encrypted_matY.size() != submatY_num)
            throw std::invalid_argument("matrix Y mismatches");

        // homomorphic automorphisms for MergeAut stage
        std::vector<std::vector<seal::Ciphertext>> cached_ct(submatY_num);

        if (submatY_num >= threads){
            auto hom_aut = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    hom_aut_galois_group(encrypted_matY[i], cached_ct[i], spp_PackRLWEs_[0], 1, true);
                }
            };

            uint32_t thread_block = (threads + submatY_num - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, submatY_num);
                thread_pool[i] = std::thread(hom_aut, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        } else {
            for (size_t i = 0; i < submatY_num; ++i)
                hom_aut_galois_group(encrypted_matY[i], cached_ct[i], spp_PackRLWEs_[0], threads, true);
        }

        encrypted_matXY.resize(submatX_num * submatY_num);
        uint32_t stride = 1UL << spp_PackRLWEs_[0].ell_;
        if (submatY_num >= threads){
            auto matmul_block = [&](size_t bgn, size_t end, size_t i){
                for (size_t j = bgn; j < end; ++j){
                    if (i != submatX_num - 1){
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_PackRLWEs_[0], 1);
                    } else{
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_PackRLWEs_[1], 1, true);
                    }
                }
            };

            uint32_t thread_block = (submatY_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < submatX_num; ++i){
                for (size_t j = 0; j < threads; ++j){
                    size_t bgn = j * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, submatY_num);
                    thread_pool[j] = std::thread(matmul_block, bgn, end, i);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }
        } else {
            for (size_t i = 0; i < submatX_num; ++i){
                for (size_t j = 0; j < submatY_num; ++j){
                    if (i != submatX_num - 1)
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_PackRLWEs_[0], threads);
                    else
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_PackRLWEs_[1], threads, true);
                }
            }
        }
    }

    void RhombusMatMul::Compute(const std::vector<seal::Plaintext> &encoded_matX,
                                const std::vector<seal::Ciphertext> &encrypted_matY,
                                std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        if (use_PackRLWEs_based_MatMul_)
            Compute_P(encoded_matX, encrypted_matY, encrypted_matXY, threads);
        else
            Compute_E(encoded_matX, encrypted_matY, encrypted_matXY, threads);
    }

    void RhombusMatMul::DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                     std::vector<seal::Plaintext> &encoded_matXY, uint32_t threads) const {
        size_t ct_num = encrypted_matXY.size();
        if (ct_num == 0)
            throw std::invalid_argument("empty ciphertexts");

        encoded_matXY.resize(ct_num);
        auto dec_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                decryptor_->decrypt(encrypted_matXY[i], encoded_matXY[i]);
            }
        };

        uint32_t thread_block = (ct_num + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, ct_num);
            thread_pool[i] = std::thread(dec_program, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusMatMul::drop_unused_coeffs_P(std::vector<seal::Ciphertext> &ct, uint32_t threads) const {
        if (!use_PackRLWEs_based_MatMul_){
            return;
        }
        uint32_t mw = 1UL << BitLength(mat_X_.ncols_ - 1);
        uint32_t ncol_pack = poly_modulus_degree_ / mw;
        uint32_t submatX_num = (mat_X_.nrows_ + mw - 1) / mw;
        uint32_t submatY_num = (mat_Y_.ncols_ + ncol_pack - 1) / ncol_pack;
        uint32_t bgn_ct = (submatX_num - 1) * submatY_num;

        // compress the sparse output
        uint32_t nrows_last = mat_X_.nrows_ - (submatX_num - 1) * mw;
        uint32_t nrows_last_pow2 = 1UL << BitLength(nrows_last - 1);
        uint32_t result_num_ct = poly_modulus_degree_ / mw;
        // the number of useful coefficients in each output ct
        uint32_t output_ct_useful_coeff_num = result_num_ct * nrows_last_pow2;
        if (output_ct_useful_coeff_num < poly_modulus_degree_){

            auto drop = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    drop_unrelated_coeffs(ct[i + bgn_ct], output_ct_useful_coeff_num, *seal_context_);
                }
            };

            uint32_t thread_block = (submatY_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < threads; ++i){
                size_t bgn = i * thread_block;
                size_t end = std::min<size_t>(bgn + thread_block, submatY_num);
                thread_pool[i] = std::thread(drop, bgn, end);
            }
            std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
        }
    }

}

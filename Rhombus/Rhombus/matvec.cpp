

#include <memory>
#include <stdexcept>
#include <string>

#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintarithmod.h"
#include "matvec.h"

namespace antchain::matvec
{
    using namespace std;
    using namespace global;

    RhombusMatVec::RhombusMatVec(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits)
    : poly_modulus_degree_(poly_degree), mod_bits_(mod_bits)
    {
        size_t coeff_mod_size = coeff_mod_bits.size();
        if (coeff_mod_size < 2 || coeff_mod_size > 4)
            throw std::invalid_argument("coeff modulus size must be 2, 3, 4");

        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_degree);
        parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_degree, coeff_mod_bits));

        seal_context_ptr = std::make_shared<seal::SEALContext>(parms);

        // copy modulus
        std::vector<uint64_t> Allmods(coeff_mod_size);
        for (size_t i = 0; i < coeff_mod_size; ++i)
            Allmods[i] = parms.coeff_modulus()[i].value();

        gen_aux_params(aux_parms_, mod_bits, Allmods, coeff_mod_size);
        keygen_ = std::make_unique<seal::KeyGenerator>(*seal_context_ptr);
        secret_key_ = keygen_->secret_key();
        keygen_->create_public_key(my_public_key_);
        decryptor_ = std::make_unique<seal::Decryptor>(*seal_context_ptr, secret_key_);

        kSPPMap_ = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5};
    }

    StatusOr<unique_ptr<RhombusMatVec>> RhombusMatVec::Create()
    {
        std::unique_ptr<RhombusMatVec> Rhbmv(new RhombusMatVec());
        if (nullptr == Rhbmv)
            return Status(util::error::NULL_PTR, "Failed to create object.");
        seal::EncryptionParameters params(seal::scheme_type::ckks);
        params.set_poly_modulus_degree(global::kPolyDegree);

        switch (global::kNumModulus)
        {
        case 2:
            params.set_coeff_modulus({global::kModulus[0], global::kModulus[1]});
            break;
        case 3:
            params.set_coeff_modulus({global::kModulus[0], global::kModulus[1], global::kModulus[2]});
            break;
        case 4:
            params.set_coeff_modulus({global::kModulus[0], global::kModulus[1], global::kModulus[2], global::kModulus[3]});
            break;
        default:
            throw std::invalid_argument("Modulus number must be 2, 3 or 4");
        }

        Rhbmv->seal_context_ptr = std::make_shared<seal::SEALContext>(params);

        // copy modulus
        std::vector<uint64_t> Allmods(kNumModulus);
        for (size_t i = 0; i < kNumModulus; ++i)
            Allmods[i] = kModulus[i];

        // generate some auxiliary parameters for encoding and decoding
        gen_aux_params(Rhbmv->aux_parms_, kSSBitLen, Allmods, kNumModulus);

        // generate secret key and public key
        Rhbmv->keygen_ = std::make_unique<seal::KeyGenerator>(*Rhbmv->seal_context_ptr);
        Rhbmv->secret_key_ = Rhbmv->keygen_->secret_key();
        Rhbmv->keygen_->create_public_key(Rhbmv->my_public_key_);
        Rhbmv->decryptor_ = std::make_unique<seal::Decryptor>(*Rhbmv->seal_context_ptr, Rhbmv->secret_key_);
        Rhbmv->poly_modulus_degree_ = global::kPolyDegree;
        Rhbmv->mod_bits_ = global::kSSBitLen;

        Rhbmv->kSPPMap_.resize(global::kLogPolyDegree + 1);
        std::copy_n(global::kSPPMap, global::kLogPolyDegree + 1, Rhbmv->kSPPMap_.begin());
        return std::move(Rhbmv);
    }

    StatusOr<unique_ptr<RhombusMatVec>> RhombusMatVec::Create(uint32_t poly_degree, uint32_t mod_bits,
                                                              const std::vector<int> &coeff_mod_bits) {
        std::unique_ptr<RhombusMatVec> Rhbmv(new RhombusMatVec());
        if (Rhbmv == nullptr)
            return Status{util::error::NULL_PTR, "Failed to create object"};
        size_t coeff_mod_size = coeff_mod_bits.size();
        if (coeff_mod_size < 2 || coeff_mod_size > 4)
            return Status{util::error::INVALID_ARGUMENT, "coeff modulus size must be 2, 3, or 4"};

        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_degree);
        parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_degree, coeff_mod_bits));

        Rhbmv->seal_context_ptr = std::make_shared<seal::SEALContext>(parms);

        // copy modulus
        std::vector<uint64_t> Allmods(coeff_mod_size);
        for (size_t i = 0; i < coeff_mod_size; ++i)
            Allmods[i] = parms.coeff_modulus()[i].value();

        gen_aux_params(Rhbmv->aux_parms_, mod_bits, Allmods, coeff_mod_size);
        Rhbmv->keygen_ = std::make_unique<seal::KeyGenerator>(*Rhbmv->seal_context_ptr);
        Rhbmv->secret_key_ = Rhbmv->keygen_->secret_key();
        Rhbmv->keygen_->create_public_key(Rhbmv->my_public_key_);
        Rhbmv->decryptor_ = std::make_unique<seal::Decryptor>(*Rhbmv->seal_context_ptr, Rhbmv->secret_key_);
        Rhbmv->poly_modulus_degree_ = poly_degree;
        Rhbmv->mod_bits_ = mod_bits;

        // We fix the spp table in this case, you could modify it manually, or call SetSPPMap method
        Rhbmv->kSPPMap_ = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5};
        return std::move(Rhbmv);
    }

    // generate pk, gk
    Status RhombusMatVec::GenKey()
    {
        try
        {
            keygen_->create_public_key(my_public_key_);

            size_t galois_elt_count = seal::util::get_power_of_two(poly_modulus_degree_);
            vector<uint32_t> galois_elt(galois_elt_count);

            // default: 13 galois elements {3,5,9,17,33,65,129,257,513,1025,2049,4097,8193} if N = 8192
            for (size_t i = 0; i < galois_elt_count; i++)
                galois_elt[i] = (1UL << (i + 1)) + 1;
            keygen_->create_galois_keys(galois_elt, galois_keys_);
        }
        CATCH_EXCEPTION()
        return Status::OK;
    }

    void RhombusMatVec::reset_secret_key(const seal::SecretKey &new_sk) {
        secret_key_ = new_sk;
        keygen_ = make_unique<seal::KeyGenerator>(*seal_context_ptr, new_sk);
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_ptr, secret_key_);
        keygen_->create_public_key(my_public_key_);
    }

    StatusOr<uint64_t> RhombusMatVec::GenPublicKey(std::string &out) const
    {
        seal::Serializable<seal::PublicKey> pk = keygen_->create_public_key();

        std::ostringstream ostr;
        uint64_t save_size;
        try
        {
            save_size = pk.save(ostr);
        }
        CATCH_EXCEPTION()
        ostr.str().swap(out);
        return save_size;
    }

    StatusOr<uint64_t> RhombusMatVec::GenPublicKey(uint8_t *buffer, size_t buffer_size) const {
        seal::Serializable<seal::PublicKey> pk = keygen_->create_public_key();

        uint64_t save_size = pk.save_size();
        if (save_size > buffer_size)
            return Status{util::error::INVALID_ARGUMENT, "buffer is not enough to save the pk"};
        try{
            save_size = pk.save((seal::seal_byte *)buffer, buffer_size);
        }CATCH_EXCEPTION()
        return Status::OK;
    }

    StatusOr<uint64_t> RhombusMatVec::GenGaloisKey(std::string &out) const
    {
        uint64_t galois_elt_count = seal::util::get_power_of_two(poly_modulus_degree_);
        vector<uint32_t> galois_elt(galois_elt_count);

        // default: 13 galois elements {3,5,9,17,33,65,129,257,513,1025,2049,4097,8193} if N = 8192
        // 12 galois elements if N = 4096
        for (size_t i = 0; i < galois_elt_count; i++)
            galois_elt[i] = (1UL << (i + 1)) + 1;
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elt);
        std::ostringstream ostr;
        size_t save_size;
        try
        {
            save_size = gk.save(ostr);
        }
        CATCH_EXCEPTION()
        ostr.str().swap(out);
        return save_size;
    }

    StatusOr<uint64_t> RhombusMatVec::GenGaloisKey(uint8_t *buffer, size_t buffer_size) const {
        uint64_t galois_elt_count = seal::util::get_power_of_two(poly_modulus_degree_);
        vector<uint32_t> galois_elt(galois_elt_count);

        for (size_t i = 0; i < galois_elt_count; ++i)
            galois_elt[i] = (1UL << (i + 1)) + 1;
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elt);
        size_t save_size = gk.save_size();
        if (save_size > buffer_size)
            return Status{util::error::INVALID_ARGUMENT, "buffer size is not enough"};
        try{
            save_size = gk.save((seal::seal_byte *)buffer, buffer_size);
        } CATCH_EXCEPTION()
        return save_size;
    }

    StatusOr<uint64_t> RhombusMatVec::GenGaloisKey(std::string &out, const std::vector<uint32_t> &galois_elts) const {
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elts);
        std::ostringstream ostr;
        size_t save_size;
        try {
            save_size = gk.save(ostr);
        }
        CATCH_EXCEPTION()
        ostr.str().swap(out);
        return save_size;
    }

    StatusOr<uint64_t> RhombusMatVec::GenGaloisKey(uint8_t *buffer, size_t buffer_size,
                                         const std::vector<uint32_t> &galois_elts) const {
        seal::Serializable<seal::GaloisKeys> gk = keygen_->create_galois_keys(galois_elts);
        size_t save_size = gk.save_size();
        if (save_size > buffer_size)
            return Status{util::error::INVALID_ARGUMENT, "buffer size is not enough"};
        try{
            save_size = gk.save((seal::seal_byte *)buffer, buffer_size);
        } CATCH_EXCEPTION()
        return save_size;
    }

    StatusOr<uint64_t> RhombusMatVec::SetPublicKey(const std::string &in)
    {
        if (in.empty())
            return Status{util::error::INVALID_ARGUMENT, "input string cannot be null"};
        uint64_t load_size = public_key_.load(*seal_context_ptr, (seal::seal_byte *)in.data(), in.size());
        return load_size;
    }

    StatusOr<uint64_t> RhombusMatVec::SetPublicKey(const uint8_t *buffer, size_t buffer_size) {
        if (buffer == nullptr)
            return Status{util::error::INVALID_ARGUMENT, "buffer is null"};
        uint64_t load_size;
        try{
            load_size = public_key_.load(*seal_context_ptr, (seal::seal_byte *)buffer, buffer_size);
        } CATCH_EXCEPTION()
        return load_size;
    }

    StatusOr<uint64_t> RhombusMatVec::SetGaloisKey(const std::string &in)
    {
        if (in.empty())
            return Status{util::error::INVALID_ARGUMENT, "input string cannot be null"};
        uint64_t load_size = galois_keys_.load(*seal_context_ptr, (seal::seal_byte *)in.data(), in.size());
        return load_size;
    }

    StatusOr<uint64_t> RhombusMatVec::SetGaloisKey(const uint8_t *buffer, size_t buffer_size) {
        if (buffer == nullptr)
            return Status{util::error::INVALID_ARGUMENT, "empty buffer"};
        uint64_t load_size;
        try{
            load_size = galois_keys_.load(*seal_context_ptr, (seal::seal_byte *)buffer, buffer_size);
        }CATCH_EXCEPTION()
        return load_size;
    }

    void RhombusMatVec::SetMatrixRowsCols(uint32_t rows, uint32_t cols, uint32_t mat_bits)
    {
        global::gen_large_mat_info(mat_info_, rows, cols, mat_bits);
    }

    StatusOr<uint64_t> RhombusMatVec::GetCiphertextSize(const seal::Ciphertext &crt_ciphertext) const
    {
        uint64_t ciphertext_size = 0;
        try
        {
            ciphertext_size = crt_ciphertext.save_size(seal::Serialization::compr_mode_default);
        }
        CATCH_EXCEPTION()

        return ciphertext_size;
    }

    StatusOr<uint64_t> RhombusMatVec::GetCiphertextSize(const vector<seal::Ciphertext> &ciphertext_vec) const
    {
        uint64_t ciphertext_size = 0;
        uint64_t vec_size = ciphertext_vec.size();
        try
        {
            for (size_t i = 0; i < vec_size; ++i)
                ciphertext_size += ciphertext_vec[i].save_size(seal::Serialization::compr_mode_default);
        }
        CATCH_EXCEPTION()
        return ciphertext_size;
    }

    StatusOr<uint64_t> RhombusMatVec::CiphertextToBytes(const seal::Ciphertext &ciphertext,
                                              uint8_t *out, uint64_t size) const
    {
        if (nullptr == out)
        {
            return Status(util::error::INVALID_ARGUMENT, "The argument cannot be NULL.");
        }
        uint64_t ciphertext_size;
        try
        {
            ciphertext_size = ciphertext.save_size(seal::Serialization::compr_mode_default);
        }
        CATCH_EXCEPTION()
        uint8_t *key_ptr = out;

        try
        {
            ciphertext.save((seal::seal_byte *)key_ptr, (size_t)ciphertext_size);
        }
        CATCH_EXCEPTION()
        return ciphertext_size;
    }

    StatusOr<uint64_t> RhombusMatVec::CiphertextToBytes(const seal::Ciphertext &ciphertext, std::string &out) const
    {
        std::ostringstream ostr;
        uint64_t ctxt_size;
        ctxt_size = ciphertext.save(ostr);
        ostr.str().swap(out);
        return ctxt_size;
    }

    StatusOr<uint64_t> RhombusMatVec::CiphertextToBytes(const vector<seal::Ciphertext> &ciphertext_vec, uint8_t *out,
                                              uint64_t size) const
    {
        if (nullptr == out)
            return Status(util::error::INVALID_ARGUMENT, "The argument cannot be NULL");

        uint64_t ciphertext_size;
        uint64_t ciphertext_vec_size;
        uint64_t ctx_num = ciphertext_vec.size();

        ciphertext_vec_size = GetCiphertextSize(ciphertext_vec).ValueOrDie();

        if (size < ciphertext_vec_size)
            return Status(util::error::INVALID_ARGUMENT, "The parameter buffer size is invalid");

        uint8_t *key_ptr = out;

        try
        {
            uint64_t remain_size = size;
            for (size_t i = 0; i < ctx_num; i++)
            {
                ciphertext_size = ciphertext_vec[i].save((seal::seal_byte *)key_ptr, remain_size);
                key_ptr += ciphertext_size;
                remain_size -= ciphertext_size;
            }
        }
        CATCH_EXCEPTION()
        return ciphertext_vec_size;
    }

    StatusOr<uint64_t> RhombusMatVec::CiphertextToBytes(const std::vector<seal::Ciphertext> &ciphertext_vec,
                                              std::vector<std::string> &out) const
    {
        uint64_t save_size = 0;
        size_t ctxt_num = ciphertext_vec.size();
        if (ctxt_num == 0)
            return Status{util::error::INVALID_ARGUMENT, "ctxt number is zero"};
        out.resize(ctxt_num);
        for (size_t i = 0; i < ctxt_num; i++)
        {
            std::ostringstream ostr;
            uint64_t cur_save_size = ciphertext_vec[i].save(ostr);
            save_size += cur_save_size;
            ostr.str().swap(out[i]);
        }
        return save_size;
    }

    StatusOr<seal::Ciphertext> RhombusMatVec::BytesToCiphertext(uint8_t *in, uint64_t size) const
    {
        if (nullptr == in)
            return Status(util::error::INVALID_ARGUMENT, "The argument cannot be NULL.");
        uint8_t *key_ptr = in;
        seal::Ciphertext ciphertext;
        try
        {
            ciphertext.unsafe_load(*seal_context_ptr, (seal::seal_byte *)key_ptr, (size_t)size);
        }
        CATCH_EXCEPTION()
        return std::move(ciphertext);
    }

    StatusOr<seal::Ciphertext> RhombusMatVec::BytesToCiphertext(const std::string &in) const
    {
        seal::Ciphertext ctxt;
        ctxt.unsafe_load(*seal_context_ptr, (seal::seal_byte *)in.data(), in.size());
        return std::move(ctxt);
    }

    StatusOr<vector<seal::Ciphertext>> RhombusMatVec::BytesToCiphertextVec(uint8_t *in, uint64_t size, uint64_t ctx_num) const
    {
        if (nullptr == in)
            return Status(util::error::INVALID_ARGUMENT, "The argument cannot be NULL.");
        uint8_t *key_ptr = in;
        vector<seal::Ciphertext> ciphertext_vec(ctx_num);
        try
        {
            size_t remain_size = size;
            size_t one_ctx_size;
            for (size_t i = 0; i < ctx_num; i++)
            {
                one_ctx_size = ciphertext_vec[i].unsafe_load(*seal_context_ptr, (seal::seal_byte *)key_ptr, remain_size);
                remain_size -= one_ctx_size;
                key_ptr += one_ctx_size;
            }
        }
        CATCH_EXCEPTION()
        return std::move(ciphertext_vec);
    }

    StatusOr<vector<seal::Ciphertext>> RhombusMatVec::BytesToCiphertextVec(const std::vector<std::string> &in) const
    {
        if (in.empty())
            return Status{util::error::INVALID_ARGUMENT, "Empty vector"};
        vector<seal::Ciphertext> ciphertext_vec(in.size());
        try
        {
            for (size_t i = 0; i < in.size(); i++)
                ciphertext_vec[i].unsafe_load(*seal_context_ptr, (seal::seal_byte *)in[i].data(), in[i].size());
        }
        CATCH_EXCEPTION()
        return std::move(ciphertext_vec);
    }

    bool RhombusMatVec::encode_to_coeff128(seal::Plaintext &plain, const uint64_t *plain_vec, size_t vec_size) const
    {
#ifdef RHOMBUS_DEBUG
        if (vec_size < 1 || vec_size > poly_modulus_degree_)
            throw std::invalid_argument("plain vector size is invalid");
#endif
        bool is_all_zero = std::all_of(plain_vec, plain_vec + (vec_size << 1),
                                       [](auto item)
                                       { return (item == 0); });
        if (is_all_zero)
            return true;
        auto parms_id = seal_context_ptr->first_parms_id();
        auto context_data = seal_context_ptr->get_context_data(parms_id);
        auto &parms = context_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = parms.coeff_modulus().size();
        auto ntt_table = context_data->small_ntt_tables();

        plain.parms_id() = seal::parms_id_zero;
        size_t buffer_size = poly_modulus_degree_ * coeff_modulus_size;
        plain.resize(buffer_size);
        std::fill_n(plain.data(), buffer_size, 0);

        for (size_t i = 0; i < coeff_modulus_size; ++i)
        {
            auto offset = i * poly_modulus_degree_;
            // 2^128 mod qi
            uint64_t pow_2_128_mod_qi = seal::util::multiply_uint_mod(aux_parms_.pow_2_64_mod_qi_[i],
                                                                      aux_parms_.pow_2_64_mod_qi_[i], coeff_modulus[i]);
            for (size_t j = 0; j < poly_modulus_degree_; ++j)
            {
                // positive integer in Z_{2^128}
                if (plain_vec[(j << 1) + 1] < (1ULL << 63))
                {
                    plain[offset + j] = seal::util::barrett_reduce_128(plain_vec + (j << 1), coeff_modulus[i]);
                }
                else
                { // A negative integer elt will be represented as 2^128 + elt, therefore we subtract the 2^128 in this case.
                    uint64_t v = seal::util::barrett_reduce_128(plain_vec + (j << 1), coeff_modulus[i]);
                    plain[offset + j] = seal::util::sub_uint_mod(v, pow_2_128_mod_qi, coeff_modulus[i]);
                }
            }
            seal::util::ntt_negacyclic_harvey_lazy(plain.data(offset), ntt_table[i]);
        }
        plain.parms_id() = parms_id;
        plain.scale() = 1.;
        return false;
    }

    Status RhombusMatVec::MatVecMul(const std::vector<seal::Plaintext> &encoded_mat, const seal::Ciphertext &ct,
                                    seal::Ciphertext &result, uint32_t threads) const {
        if (encoded_mat.empty())
            return {util::error::INVALID_ARGUMENT, "encoded_mat is empty"};
        return MatVecMulInternal(ct, encoded_mat, mat_info_.corner_block_array_[0], result, threads);
    }

    Status RhombusMatVec::MatVecMulInternal(const seal::Ciphertext &ct, const std::vector<seal::Plaintext> &encoded_mat,
                                            const global::MatInfo &mat_info, seal::Ciphertext &result,
                                            uint32_t threads) const
    {
        // check matrix
        if (mat_info.nrows_ > poly_modulus_degree_ || mat_info.ncols_ > poly_modulus_degree_)
            return {util::error::INVALID_MATRIX_ROWS_COLS, "matrix rows or columns out of bound"};
        if (mat_info.nrows_ == 0 || mat_info.ncols_ == 0)
            return {util::error::INVALID_MATRIX_ROWS_COLS, "the number of rows, columns should not be zero"};


        size_t g0_size = (size_t)1 << (mat_info.spp_parms_.u_ - mat_info.spp_parms_.h_);
        size_t g1_size = (size_t)1 << (mat_info.log_pad_nrows_ - mat_info.spp_parms_.u_);

        // MergeAut: G0 = Gal(K_u/K_h)
        std::vector<seal::Ciphertext> cached_ct(g0_size);
        hom_aut_galois_group(ct, cached_ct, mat_info, threads, true);

        // save the inner products
        std::vector<seal::Ciphertext> matvec(g1_size);

        auto hom_mul_pt = [&](size_t bgn, size_t end){
            seal::Ciphertext tmp;
            for (size_t i1 = bgn; i1 < end; ++i1){
                set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id());
                for (size_t j = 0; j < g0_size; ++j){
                    multiply_plain_ntt(cached_ct[j], encoded_mat[i1 * g0_size + j], tmp, *seal_context_ptr);
                    add_inplace(matvec[i1], tmp, *seal_context_ptr);
                }
                // for PackRLWEs
                transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                while (matvec[i1].coeff_modulus_size() > remain_mod_num_){
                    rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                }
            }
        };

        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(hom_mul_pt, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});

        PackRLWEs(matvec, mat_info.spp_parms_.u_, galois_keys_, result, *seal_context_ptr, threads);
        transform_to_ntt_inplace(result, *seal_context_ptr);

        return Status::OK;
    }

} // namespace

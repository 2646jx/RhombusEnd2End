#include "matmul_protocol.h"

namespace antchain {
    using namespace std;

    RhombusLinear::RhombusLinear(sci::NetIO *io, int party, uint32_t N, uint32_t mod_bits,
                                 const std::vector<int> &coeff_mod_bits)
            : io_(io), party_(party), mod_bits_(mod_bits) {

        matmul_impl_ = std::make_shared<RhombusMatMul>(N, mod_bits, coeff_mod_bits);
        matmul_impl_->set_method(true);
        half_mod_ = 1ULL << (mod_bits_ - 1);
        base_mod_ = 1ULL << mod_bits_;
        mod_mask_ = base_mod_ - 1;

        matvec_impl_ = std::make_shared<RhombusMatVec>(N, mod_bits, coeff_mod_bits);

        if (mod_bits <= 40){
            set_remain_mod_num(1);
        }else{
            set_remain_mod_num(2);
        }

        if (party_ == 1)
        {
            recv_keys();
            send_keys();
        }
        else
        {
            send_keys();
            recv_keys();
        }
    }

    void RhombusLinear::send_keys() {
        std::string pk, gk;
        matmul_impl_->GenPublicKey(pk);
        matmul_impl_->GenGaloisKey(gk);

        // keep sk the same
        matvec_impl_->reset_secret_key(matmul_impl_->secret_key());

        uint32_t pk_size = pk.size();
        uint32_t gk_size = gk.size();
        io_->send_data(&pk_size, sizeof(uint32_t));
        io_->send_data(pk.c_str(), pk_size);
        io_->send_data(&gk_size, sizeof(uint32_t));
        io_->send_data(gk.c_str(), gk_size);
    }

    void RhombusLinear::recv_keys() {
        std::stringstream is_pk, is_gk;
        uint32_t pk_size, gk_size;
        io_->recv_data(&pk_size, sizeof(uint32_t));
        char *c_pk = new char[pk_size];
        io_->recv_data(c_pk, pk_size);
        is_pk.write(c_pk, pk_size);
        string pk = is_pk.str();
        matvec_impl_->SetPublicKey(pk);
        matmul_impl_->public_key().unsafe_load(matmul_impl_->seal_context(), is_pk);

        delete[] c_pk;

        io_->recv_data(&gk_size, sizeof(uint32_t));
        char *c_gk = new char[gk_size];
        io_->recv_data(c_gk, gk_size);
        is_gk.write(c_gk, gk_size);
        string gk = is_gk.str();
        matvec_impl_->SetGaloisKey(gk);
        matmul_impl_->galois_key().unsafe_load(matmul_impl_->seal_context(), is_gk);
        delete[] c_gk;
    }

    void RhombusLinear::send_ciphertext(const seal::Ciphertext &ct) {
        std::stringstream os;
        uint32_t ct_size;
        ct.save(os);
        ct_size = os.tellp();
        string ct_ser = os.str();
        io_->send_data(&ct_size, sizeof(uint32_t));
        io_->send_data(ct_ser.c_str(), ct_ser.size());
    }

    void RhombusLinear::send_ciphertext(const std::string &ct) {
        uint32_t ct_size = ct.size();
        io_->send_data(&ct_size, sizeof(uint32_t));
        io_->send_data(ct.c_str(), ct_size);
    }

    void RhombusLinear::send_encrypted_vector(const std::vector<std::string> &ct) {
        uint32_t ct_size = ct.size();
        uint32_t cur_ct_size;
        io_->send_data(&ct_size, sizeof(uint32_t));
        for (auto &c: ct) {
            cur_ct_size = c.size();
            io_->send_data(&cur_ct_size, sizeof(uint32_t));
            io_->send_data(c.c_str(), c.size());
        }
    }

    void RhombusLinear::recv_ciphertext(seal::Ciphertext &ct) {
        std::stringstream is;
        uint32_t ct_size;
        io_->recv_data(&ct_size, sizeof(uint32_t));
        char *c_enc_result = new char[ct_size];
        io_->recv_data(c_enc_result, ct_size);
        is.write(c_enc_result, ct_size);
        ct.unsafe_load(matmul_impl_->seal_context(), is);
        delete[] c_enc_result;
    }

    void RhombusLinear::send_encrypted_vector(const std::vector<seal::Ciphertext> &ct_vec) {
        uint32_t ncts = ct_vec.size();
        io_->send_data(&ncts, sizeof(uint32_t));
        for (size_t i = 0; i < ncts; ++i) {
            send_ciphertext(ct_vec[i]);
        }
    }

    void RhombusLinear::recv_encrypted_vector(std::vector<seal::Ciphertext> &ct_vec) {
        uint32_t ncts{0};
        io_->recv_data(&ncts, sizeof(uint32_t));
        ct_vec.resize(ncts);
        for (size_t i = 0; i < ncts; ++i) {
            recv_ciphertext(ct_vec[i]);
        }
    }

    int64_t RhombusLinear::to_signed(uint64_t a) const
    {
        return (a >= half_mod_) ? int64_t(a - base_mod_) : int64_t(a);
    }

    uint64_t RhombusLinear::to_unsigned(int64_t a) const
    {
        return ((uint64_t)a) & mod_mask_;
    }

    void RhombusLinear::to_signed(const uint64_t *inputArr, size_t len, int64_t *outArr) const
    {
        auto transform_func = [&](size_t bgn, size_t end)
        {
            for (size_t i = bgn; i < end; ++i)
            {
                uint64_t a = inputArr[i];
                outArr[i] = (a >= half_mod_) ? int64_t(a - base_mod_) : int64_t(a);
            }
        };

        // multi-thread
        uint64_t thread_block = (len + threads_ - 1) / threads_; // ceil(len/num_threads)
        std::vector<std::thread> thread_pool(threads_);
        for (size_t i = 0; i < threads_; ++i)
        {
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, len);
            thread_pool[i] = std::thread(transform_func, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusLinear::to_unsigned(const int64_t *inputArr, size_t len, uint64_t *outArr) const
    {
        auto transform_func = [&](size_t bgn, size_t end)
        {
            for (size_t i = bgn; i < end; ++i)
            {
                // implicit conversion to uint64_t
                outArr[i] = inputArr[i] & mod_mask_;
            }
        };

        // multi-thread
        uint64_t thread_block = (len + threads_ - 1) / threads_; // ceil(len/num_threads)
        std::vector<thread> thread_pool(threads_);
        for (size_t i = 0; i < threads_; ++i)
        {
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, len);
            thread_pool[i] = std::thread(transform_func, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    void RhombusLinear::SignedMatMulMod(const int64_t *matX, const int64_t *matY, int64_t *matXY, uint32_t n, uint32_t m, uint32_t k) const
    {
        auto matmul_func = [&](size_t bgn, size_t end)
        {
            for (size_t i = bgn; i < end; ++i)
            {
                for (size_t j = 0; j < k; ++j)
                {
                    int64_t tmp = 0;
                    for (size_t l = 0; l < m; ++l)
                    {
                        tmp += matX[i * m + l] * matY[l * k + j];
                    }
                    matXY[i * k + j] = to_signed((uint64_t)tmp & mod_mask_);
                }
            }
        };

        // multi-thread
        uint32_t thread_block = (n + threads_ - 1) / threads_;
        std::vector<thread> thread_pool(threads_);
        for (size_t i = 0; i < threads_; ++i)
        {
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, n);
            thread_pool[i] = std::thread(matmul_func, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    // mat is held by server, vec is shared
    void RhombusLinear::MatVecMul(uint32_t nr, uint32_t nc, uint64_t *mat, uint64_t *vec, uint64_t *mv)
    {
        if (party_ == 1) { // SERVER
            matvec_impl_->SetMatrixRowsCols(nr, nc, mod_bits_);
            vector<seal::Plaintext> encoded_mat;

            // to signed mat
            auto *signed_mat = new int64_t[nr * nc];
            to_signed(mat, nr * nc, signed_mat);

            // encode matrix
            matvec_impl_->EncodeMat(signed_mat, encoded_mat, threads_);

            seal::Plaintext encoded_vec_s;
            if (vec != nullptr){
                encode_to_coeff(encoded_vec_s, vec, nc, Ecd_SCALED_INVERSE_ORDER, matvec_impl_->aux_parms(),
                                matvec_impl_->seal_context().first_parms_id(), matvec_impl_->seal_context(), mod_bits_);
            }else{
                // encode zero
                vector<uint64_t> zero(nc, 0);
                encode_to_coeff(encoded_vec_s, zero.data(), nc, Ecd_SCALED_INVERSE_ORDER, matvec_impl_->aux_parms(),
                                matvec_impl_->seal_context().first_parms_id(), matvec_impl_->seal_context(), mod_bits_);
            }

            seal::Ciphertext enc_vec_c;
            recv_ciphertext(enc_vec_c);

            // merge share
            seal::Ciphertext enc_vec;
            add_plain(enc_vec_c, encoded_vec_s, enc_vec, matvec_impl_->seal_context());

            // compute
            seal::Ciphertext enc_mv, enc_mv_c;
            matvec_impl_->MatVecMul(encoded_mat, enc_vec, enc_mv, threads_);
            matvec_impl_->ConvToSS(enc_mv, enc_mv_c, mv, nr, Dcd_SCALED_STRIDE);
            matvec_impl_->drop_unused_coeffs(enc_mv_c, nr);

            send_ciphertext(enc_mv_c);

            delete[] signed_mat;

        }
        else{ // CLIENT
            matvec_impl_->SetMatrixRowsCols(nr, nc, mod_bits_);

            string enc_vec_c;
            matvec_impl_->EncryptVec(vec, nc, enc_vec_c, Ecd_SCALED_INVERSE_ORDER);

            send_ciphertext(enc_vec_c);

            seal::Ciphertext enc_mv_c;
            recv_ciphertext(enc_mv_c);

            matvec_impl_->DecryptVec(enc_mv_c, mv, nr, Dcd_SCALED_STRIDE);
        }
    }

    // X * Y
    void RhombusLinear::MatMul(uint32_t n, uint32_t m, uint32_t k, uint64_t *matX, uint64_t *matY,
                               uint64_t *outMat)
    {
        // X: n * m, holds by server,
        // Y: m * k, holds by client,
        // client encrypts Y, then sends it to server, who computes X*E(Y) homomorphically.

        // set the dimensions, required
        matmul_impl_->SetMatDims(n, m, k, mod_bits_, mod_bits_);

        if (party_ == 1) { // SERVER

            // use int64_t to represent matrix
            auto *signed_matX = new int64_t[n * m];

            to_signed(matX, n * m, signed_matX);
        
            vector<seal::Plaintext> encoded_mat;
            matmul_impl_->EncodeMatX(signed_matX, encoded_mat, threads_);

            vector<seal::Ciphertext> encrypted_matY;
            recv_encrypted_vector(encrypted_matY);

            vector<seal::Ciphertext> encrypted_matXY;
            vector<seal::Ciphertext> encrypted_matXY_c;
            matmul_impl_->Compute(encoded_mat, encrypted_matY, encrypted_matXY, threads_);
            matmul_impl_->H2A(encrypted_matXY, encrypted_matXY_c, outMat, threads_);
            matmul_impl_->drop_unused_coeffs_P(encrypted_matXY_c, threads_);

            send_encrypted_vector(encrypted_matXY_c);
            delete[] signed_matX;


        }else{ // CLIENT

            // use int64_t to represent matrix
            auto *signed_matY = new int64_t[m * k];
            to_signed(matY, m * k, signed_matY);
        
            vector<std::string> encrypted_matY;
            matmul_impl_->EncryptMatY(signed_matY, encrypted_matY, threads_);

            send_encrypted_vector(encrypted_matY);

            vector<seal::Ciphertext> encrypted_matXY_c;
            recv_encrypted_vector(encrypted_matXY_c);

            matmul_impl_->DecryptMatXY(encrypted_matXY_c, outMat, threads_);

            delete[] signed_matY;
        }
    }

    void RhombusLinear::Conv2D_1m1Filter(uint32_t H, uint32_t W, uint32_t CI, uint32_t CO, uint32_t stride,
                                         uint64_t *inputArr, uint64_t *filterArr, uint64_t *outArr)
    {

        uint32_t newH = (H + stride - 1) / stride;
        uint32_t newW = (W + stride - 1) / stride;
        uint32_t newHW = newH * newW;
        bool client_encrypt = true;
        if (newHW > CO){
            client_encrypt = false;
        }

        if (party_ == 1){

            // extract filter matrix
            vector<int64_t> filter_mat(CI * CO);
            if (client_encrypt){ // transpose
                // CO * CI matrix
                for (uint32_t i = 0; i < CO; ++i){
                    for (uint32_t j = 0; j < CI; ++j){
                        filter_mat[i * CI + j] = to_signed(filterArr[j * CO + i]);
                    }
                }
            }else {
                // CI * CO
                for (uint32_t i = 0; i < CI; ++i){
                    for (uint32_t j = 0; j < CO; ++j){
                        filter_mat[i * CO + j] = to_signed(filterArr[i * CO + j]);
                    }
                }
            }

            // extract image share matrix
            vector<int64_t> image_share_mat(newHW * CI);
            if (client_encrypt){
                // CI * newHW
                for (uint32_t i = 0; i < CI; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            image_share_mat[i * newHW + j0 * newW + j1] =
                                    to_signed(inputArr[stride * j0 * CI * W + stride * j1 * CI + i]);
                        }
                    }
                }
            } else{
                // newHW * CI
                for (uint32_t i0 = 0; i0 < newH; ++i0){
                    for (uint32_t i1 = 0; i1 < newW; ++i1){
                        for (uint32_t j = 0; j < CI; ++j){
                            image_share_mat[(i0 * newW + i1) * CI + j] =
                                    to_signed(inputArr[stride * i0 * CI * W + stride * i1 * CI + j]);
                        }
                    }
                }
            }

            if (client_encrypt){
                matmul_impl_->SetMatDims(CO, CI, newHW, mod_bits_, mod_bits_);
                vector<int64_t> filter_m_image(newHW * CO);
                SignedMatMulMod(filter_mat.data(), image_share_mat.data(), filter_m_image.data(), CO, CI, newHW);

                vector<seal::Plaintext> encoded_filter;
                matmul_impl_->EncodeMatX(filter_mat.data(), encoded_filter, threads_);

                vector<seal::Ciphertext> encrypted_image;
                recv_encrypted_vector(encrypted_image);

                vector<int64_t> filter_image_s(newHW * CO);
                vector<seal::Ciphertext> encrypted_filter_image_c;
                vector<seal::Ciphertext> encrypted_filter_image;
                
                matmul_impl_->Compute(encoded_filter, encrypted_image, encrypted_filter_image, threads_);
                matmul_impl_->H2A(encrypted_filter_image, encrypted_filter_image_c, filter_image_s.data(), threads_);
                matmul_impl_->drop_unused_coeffs_P(encrypted_filter_image_c, threads_);

                send_encrypted_vector(encrypted_filter_image_c);

                AddVecMod(filter_m_image.data(), filter_image_s.data(), filter_image_s.data(),
                          newHW * CO, mod_bits_);

                // convert
                for (uint32_t i = 0; i < CO; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            outArr[j0 * newW * CO + j1 * CO + i] =
                                    ((uint64_t)filter_image_s[i * newHW + j0 * newW + j1]) & mod_mask_;
                        }
                    }
                }
            }else{
                matmul_impl_->SetMatDims(newHW, CI, CO, mod_bits_, mod_bits_);
                vector<int64_t> image_m_filter(newHW * CO);
                SignedMatMulMod(image_share_mat.data(), filter_mat.data(), image_m_filter.data(), newHW, CI, CO);

                vector<string> encrypted_filter;
                matmul_impl_->EncryptMatY(filter_mat.data(), encrypted_filter, threads_);

                send_encrypted_vector(encrypted_filter);

                vector<seal::Ciphertext> encrypted_conv;
                recv_encrypted_vector( encrypted_conv);

                vector<int64_t> image_filter_s(newHW * CO);
                matmul_impl_->DecryptMatXY(encrypted_conv, image_filter_s.data(), threads_);

                AddVecMod(image_m_filter.data(), image_filter_s.data(), image_filter_s.data(),
                          newHW * CO, mod_bits_);

                // convert
                for (uint32_t i = 0; i < CO; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            outArr[j0 * newW * CO + j1 * CO + i] =
                                    ((uint64_t)image_filter_s[(j0 * newW + j1) * CO + i]) & mod_mask_;
                        }
                    }
                }
            }

        } else{

            // extract image matrix
            vector<int64_t> image_share_mat(newHW * CI);
            if (client_encrypt){
                // CI * newHW
                for (uint32_t i = 0; i < CI; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            image_share_mat[i * newHW + j0 * newW + j1] =
                                    to_signed(inputArr[stride * j0 * CI * W + stride * j1 * CI + i]);
                        }
                    }
                }
            } else{
                // newHW * CI
                for (uint32_t i0 = 0; i0 < newH; ++i0){
                    for (uint32_t i1 = 0; i1 < newW; ++i1){
                        for (uint32_t j = 0; j < CI; ++j){
                            image_share_mat[(i0 * newW + i1) * CI + j] =
                                    to_signed(inputArr[stride * i0 * CI * W + stride * i1 * CI + j]);
                        }
                    }
                }
            }

            if (client_encrypt){
                matmul_impl_->SetMatDims(CO, CI, newHW, mod_bits_, mod_bits_);
                vector<string> encrypted_image;
                matmul_impl_->EncryptMatY(image_share_mat.data(), encrypted_image, threads_);

                send_encrypted_vector(encrypted_image);

                vector<seal::Ciphertext> encrypted_filter_image_c;
                recv_encrypted_vector(encrypted_filter_image_c);

                vector<int64_t> filter_image_c(newHW * CO);
                matmul_impl_->DecryptMatXY(encrypted_filter_image_c, filter_image_c.data(), threads_);

                // convert
                for (uint32_t i = 0; i < CO; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            outArr[j0 * newW * CO + j1 * CO + i] =
                                    ((uint64_t)filter_image_c[i * newHW + j0 * newW + j1]) & mod_mask_;
                        }
                    }
                }
            }else{
                matmul_impl_->SetMatDims(newHW, CI, CO, mod_bits_, mod_bits_);

                vector<seal::Plaintext> encoded_image;
                matmul_impl_->EncodeMatX(image_share_mat.data(), encoded_image, threads_);

                vector<seal::Ciphertext> encrypted_filter;
                recv_encrypted_vector(encrypted_filter);

                vector<seal::Ciphertext> encrypted_image_filter_s;
                vector<int64_t> image_filter_c(newHW * CO);
                vector<seal::Ciphertext> encrypted_image_filter;
                matmul_impl_->Compute(encoded_image, encrypted_filter, encrypted_image_filter, threads_);
                matmul_impl_->H2A(encrypted_image_filter, encrypted_image_filter_s, image_filter_c.data(), threads_);
                
                matmul_impl_->drop_unused_coeffs_P(encrypted_image_filter_s, threads_);

                send_encrypted_vector(encrypted_image_filter_s);

                // convert
                for (uint32_t i = 0; i < CO; ++i){
                    for (uint32_t j0 = 0; j0 < newH; ++j0){
                        for (uint32_t j1 = 0; j1 < newW; ++j1){
                            outArr[j0 * newW * CO + j1 * CO + i] =
                                    ((uint64_t)image_filter_c[(j0 * newW + j1) * CO + i]) & mod_mask_;
                        }
                    }
                }
            }
        }
    }

}
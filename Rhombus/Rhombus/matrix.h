#ifndef RHOMBUS_MATINFO_H
#define RHOMBUS_MATINFO_H

#include <vector>
#include <thread>
#include <cassert>
#include "parameters.h"
#include "ring_op.h"


namespace antchain::global
{

    struct SPPParams
    {
        // Tr_{ell + h, h} = Tr_{ell + h, u} * Tr_{u, h}
        uint32_t ell_;
        uint32_t h_;
        uint32_t u_;

        uint32_t poly_degree_;

        uint32_t PackRLWEs_factor_;

        // Galois elements in Gal(K_u/K_h)
        std::vector<uint32_t> Gal_u_h;

        // A table of g0_size * g0_size, where g0_size = #Gal(K_u/K_h), used in matrix-preprocessing
        // (i, j)-entry is tau_{i}^{-1}((N/2^u) * j) mod 2N, where tau_i \in Gal(K_u/K_h), j \in [g0_size]
        std::vector<std::vector<uint32_t>> inv_ge_mul_powx_;
    };

    struct SPPParams_Expand
    {
        // Tr_{logN, logN - ell} = Tr_{u, logN - ell} * Tr_{logN, u}
        // ell = log(c)
        uint32_t ell_;
        uint32_t u_;

        uint32_t poly_degree_;
        uint32_t Expand_factor_;

        // Galois elements in Gal(K_u/K_{logN - ell})
        std::vector<uint32_t> Gal_u_logN_ell_;
        std::vector<uint32_t> inv_Gal_u_logN_ell_;

        // g0_size = c*2^u/N, g1_size = N/2^u
        // A table of g0_size * g0_size
        // (i, j)-entry is tau_{i}(X^{-(N/2^u) * j}), tau_i \in Gal(K_u/K_{logN - ell}), j \in [g0_size]
        std::vector<std::vector<uint32_t>> ge_mul_powx_;
    };

    // Matrix Information
    struct MatInfo
    {
        uint32_t nrows_;
        uint32_t ncols_;

        uint32_t log_pad_nrows_;
        uint32_t log_pad_ncols_;

        // bit-length of the elements
        uint32_t mat_bits_;

        SPPParams spp_parms_;
    };

    // A LightWeight matrix representation
    struct MatInfoLW
    {
        uint32_t nrows_;
        uint32_t ncols_;

        uint32_t log_pad_nrows_;
        uint32_t log_pad_ncols_;

        uint32_t mat_bits_;
    };


    struct LargeMatInfo
    {
        uint32_t nrows_;
        uint32_t ncols_;

        uint32_t mat_bits_;

        // ceil(nrows / N)
        uint32_t nrow_block_num_;

        // ceil(ncols / N)
        uint32_t ncol_block_num_;

        // 0, 1, 2, 3 <--> block[0][0], block[0][last_col], block[last_row][0], block[last_row][last_col]
        std::array<MatInfo, 4> corner_block_array_;
    };

    void gen_mat_info(MatInfo &mat_info, uint32_t rows, uint32_t cols,
                      uint32_t mat_bits = kSSBitLen, uint32_t poly_degree = kPolyDegree);

    void gen_lw_mat_info(MatInfoLW &mat_info, uint32_t rows, uint32_t cols,
                         uint32_t mat_bits = kSSBitLen, uint32_t poly_degree = kPolyDegree);

    void gen_spp_parms(SPPParams &spp_parms, uint32_t log_pad_rows, uint32_t log_pad_cols,
                       uint32_t poly_degree = kPolyDegree, const int *SPP_table = kSPPMap);

    // nrow = N
    void gen_spp_expand(SPPParams_Expand &spp_expand, uint32_t log_pad_cols, uint32_t poly_degree = kPolyDegree,
                        const int *SPP_table = kSPPMap);

    void gen_large_mat_info(LargeMatInfo &mat_info, uint32_t rows, uint32_t cols,
                            uint32_t mat_bits = kSSBitLen, uint32_t poly_degree = kPolyDegree);

    // Input packing
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_concat(const T *matrix, uint32_t N, std::vector<std::vector<T>> &con_mat, const MatInfo &mat_info)
    {
        size_t nrows_pad = size_t(1) << size_t(mat_info.log_pad_nrows_);
        size_t ncols_pad = size_t(1) << size_t(mat_info.log_pad_ncols_);

        size_t vec_num = nrows_pad * ncols_pad / N;
        size_t concat_num = N / ncols_pad;
        con_mat.resize(vec_num);
        for (size_t i = 0; i < vec_num; ++i){
            con_mat[i].resize(N);
            std::fill(con_mat[i].begin(), con_mat[i].end(), 0);
            for (size_t j = 0; j < concat_num; ++j){
                size_t index = i + j * vec_num;
                if (index < mat_info.nrows_)
                    std::copy_n(matrix + (i + j * vec_num) * mat_info.ncols_,
                            mat_info.ncols_, con_mat[i].begin() + j * ncols_pad);
            }
        }
    }

    // Input packing
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_concat(const std::vector<const T *> &matrix, uint32_t N, std::vector<std::vector<T>> &con_mat, const MatInfo &mat_info)
    {
        size_t nrows_pad = size_t(1) << size_t(mat_info.log_pad_nrows_);
        size_t ncols_pad = size_t(1) << size_t(mat_info.log_pad_ncols_);

        size_t vec_num = nrows_pad * ncols_pad / N;
        size_t concat_num = N / ncols_pad;
        con_mat.resize(vec_num);

        for (size_t i = 0; i < vec_num; ++i){
            con_mat[i].resize(N);
            std::fill(con_mat[i].begin(), con_mat[i].end(), 0);
            for (size_t j = 0; j < concat_num; ++j){
                size_t index = i + j * vec_num;
                if (index < mat_info.nrows_)
                    std::copy_n(matrix[index], mat_info.ncols_, con_mat[i].begin() + j * ncols_pad);
            }
        }
    }

    // Input: nrows * ncols matrix;
    // Output: Pow2Pad(nrows) * ncols matrix, where Pow2Pad(a) is the minimum 2-power integer >= a.
    // zero padding
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_col_padding(const std::vector<const T *> &matrix, std::vector<std::vector<T>> &pad_mat,
                            uint32_t nrows, uint32_t ncols)
    {
        uint32_t pad_rows = (uint32_t)1 << BitLength(nrows - 1);
        pad_mat.resize(pad_rows);
        for (size_t i = 0; i < pad_rows; ++i){
            pad_mat[i].resize(ncols);
            std::fill_n(pad_mat[i].begin(), ncols, 0);
            if (i < nrows)
                std::copy_n(matrix[i], ncols, pad_mat[i].begin());
        }
    }


    // For SPP optimization
    // Note: The elements will become slightly larger. Ensure that the capacity of type U is sufficient.
    //       If the elements of the preprocessed matrix may exceed 64-bit, use 'matrix_row_combine128' instead.
    // return: result is 3-dim vector shaped like:
    // {[-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----],
    //  [-----], [-----], [-----], [-----]}
    // which could be also viewed as an array of g1_size * g0_size, each element is a vector
    template <typename T, typename U,
            typename std::enable_if<std::is_signed_v<T> && std::is_signed_v<U>, T>::type * = nullptr>
    void matrix_row_combine64(const std::vector<std::vector<T>> &con_mat, const SPPParams &spp,
                              std::vector<std::vector<std::vector<U>>> &result, uint32_t threads,
                              uint32_t poly_degree = kPolyDegree, uint32_t ncols = kPolyDegree)
    {
        size_t g1_size = (size_t)1 << (spp.ell_ + spp.h_ - spp.u_);
        size_t g0_size = (size_t)1 << (spp.u_ - spp.h_);
        result.resize(g1_size);

        auto m_mask{(poly_degree << 1) - 1};
        auto n_mask{poly_degree - 1};

        // compute \sum_{i0=0}^{2^{u-h}-1} X^{(N/2^u)*i0} * tau_j0(y_{i0 * g1_size + i1}) for i1\in [bgn, end), tau_j0 \in Gal(K_u/K_h)
        // y_i is the i-th row of the concatenated matrix
        auto mat_pp_slice = [&](size_t bgn, size_t end){
            std::vector<U> acc_vec(poly_degree);
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(g0_size);
                for (size_t j0 = 0; j0 < g0_size; ++j0){
                    result[i1][j0].resize(poly_degree);
                    std::fill_n(acc_vec.data(), poly_degree, 0);
                    for (uint32_t i0 = 0; i0 < g0_size; ++i0){
                        // tau_{j0}^{-1}(X^{(N/2^u)*i0})
                        auto inv_powx = spp.inv_ge_mul_powx_[j0][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign)
                            inv_powx &= n_mask;
                        auto rem_inv_pow = poly_degree - inv_powx;


                        // y_{i0*g1_size+i1} * tau_{j0}^{-1}(X^{(N/2^u)*i0})
                        // [--------------------========] --> [        --------------------]
//                        std::transform(con_mat[i0 * g1_size + i1].cbegin(), con_mat[i0 * g1_size + i1].cbegin() + rem_inv_pow,
//                                       acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
//                                    return (sign ? -elt0 : elt0) + elt1;});
                        // [--------------------========] --> [========                    ]
//                        std::transform(con_mat[i0 * g1_size + i1].cbegin() + rem_inv_pow, con_mat[i0 * g1_size + i1].cend(),
//                                       acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
//                                    return (sign ? elt0 : -elt0) + elt1;});

                        if (ncols <= rem_inv_pow){
                            std::transform(con_mat[i0 * g1_size + i1].cbegin(), con_mat[i0 * g1_size + i1].cbegin() + ncols,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return (sign ? -elt0 : elt0) + elt1;});
                        }else{
                            std::transform(con_mat[i0 * g1_size + i1].cbegin(), con_mat[i0 * g1_size + i1].cbegin() + rem_inv_pow,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return (sign ? -elt0 : elt0) + elt1;});

                            std::transform(con_mat[i0 * g1_size + i1].cbegin() + rem_inv_pow, con_mat[i0 * g1_size + i1].cbegin() + ncols,
                                           acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
                                        return (sign ? elt0 : -elt0) + elt1;});
                        }
                    }

                    // final tau_j0:
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i)
                    {
                        uint32_t index = index_raw & m_mask;
                        if (index >= poly_degree)
                            result[i1][j0][index & n_mask] = -acc_vec[i];
                        else
                            result[i1][j0][index] = acc_vec[i];
                        index_raw += spp.Gal_u_h[j0];
                    }
                }
            }
        };

        // multi-thread implementation
        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(mat_pp_slice, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

        // transpose the matrix, then pad the number of rows (of the transposed matrix) to 2-power
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_transpose_and_padding(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols, uint32_t pad_cols,
                                      std::vector<std::vector<T>> &new_mat, uint32_t poly_degree = kPolyDegree)
    {
        assert(pad_cols >= ncols);
        assert(nrows <= poly_degree);
        assert(ncols <= poly_degree);

        if (matrix.empty())
            throw std::invalid_argument("empty matrix");
        new_mat.resize(pad_cols);
        for (size_t i = 0; i < pad_cols; ++i){
            // zero-padding
            new_mat[i].resize(nrows);
            std::fill_n(new_mat[i].begin(), nrows, 0);
            if (i < ncols){
                // transpose
                for (size_t j = 0; j < nrows; ++j){
                    new_mat[i][j] = matrix[j][i];
                }
            }
        }
    }

    /*!
     * Denote each column of the matrix (zero-pad to length-N) as w_i, i \in [pad_cols], if ncols < pad_cols,
     * zero-pad to pad_cols columns. This function outputs a table of dimension g1_size * g0_size, each element is an length-N vector.
     * Specifically, y_{i1, tau} = sum_{i0=0}^{c*2^u/N-1} tau^{-1}(w_{i0(N/2^u)+i1} * X^{-i0*(N/2^u)}), i1\in [g1_size],
     * tau\in[g0_size] (corresponds to Gal(Ku/K_{logN-ell})). Refer to Rhombus.
     * @param matrix : input
     * @param spp : required SPP parameters
     * @param ncols : the number of columns
     * @param pad_cols : padded columns
     * @param result : output
     * @param threads : threads
     * @param poly_degree : HE parameters, ring dimension.
     */
    template <typename T, typename U,
            typename std::enable_if<std::is_signed_v<T> && std::is_signed_v<U>, T>::type * = nullptr>
    void matrix_col_combine64(const std::vector<const T *> &matrix, const SPPParams_Expand &spp, uint32_t ncols, uint32_t pad_cols,
                              std::vector<std::vector<std::vector<U>>> &result, uint32_t threads, uint32_t poly_degree = kPolyDegree)
    {
        assert(pad_cols >= ncols);
        if (matrix.empty())
            throw std::invalid_argument("empty matrix");

        // G0 = Gal(Ku/K_{logN-ell}), G1 = Gal(K/K_{u})
        size_t g1_size = poly_degree >> spp.u_;
        size_t g0_size = (1UL << (spp.u_ + spp.ell_)) / poly_degree;

        // transpose and pad
        std::vector<std::vector<T>> t_mat;
        uint32_t nrows = matrix.size();
        matrix_transpose_and_padding(matrix, nrows, ncols, pad_cols, t_mat, poly_degree);

        result.resize(g1_size);
        auto m_mask{(poly_degree << 1) - 1};
        auto n_mask{poly_degree - 1};

        auto mat_pp_slice = [&](size_t bgn, size_t end){
            std::vector<U> acc_vec(poly_degree);
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(g0_size);
                for (size_t j = 0; j < g0_size; ++j){ // traverse G0
                    result[i1][j].resize(poly_degree);
                    std::fill_n(acc_vec.data(), poly_degree, 0);

                    for (uint32_t i0 = 0; i0 < g0_size; ++i0){
                        auto powx = spp.ge_mul_powx_[j][i0]; // X^powx = tau_j(X^{-i0 * (N/2^u)})
                        bool sign = (powx >= poly_degree);
                        if (sign) powx &= n_mask;
                        auto rem_pow = poly_degree - powx;

//                        std::transform(t_mat[i0 * g1_size + i1].cbegin(), t_mat[i0 * g1_size + i1].cbegin() + rem_pow,
//                                       acc_vec.cbegin() + powx, acc_vec.begin() + powx, [&](U elt0, U elt1)->U{
//                                    return (sign ? -elt0 : elt0) + elt1;});
//                        std::transform(t_mat[i0 * g1_size + i1].cbegin() + rem_pow, t_mat[i0 * g1_size + i1].cend(),
//                                       acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
//                                    return (sign ? elt0 : -elt0) + elt1;});

                        if (nrows <= rem_pow){
                            std::transform(t_mat[i0 * g1_size + i1].cbegin(), t_mat[i0 * g1_size + i1].cbegin() + nrows,
                                           acc_vec.cbegin() + powx, acc_vec.begin() + powx, [&](U elt0, U elt1)->U{
                                        return (sign ? -elt0 : elt0) + elt1;});
                        }else{
                            std::transform(t_mat[i0 * g1_size + i1].cbegin(), t_mat[i0 * g1_size + i1].cbegin() + rem_pow,
                                           acc_vec.cbegin() + powx, acc_vec.begin() + powx, [&](U elt0, U elt1)->U{
                                        return (sign ? -elt0 : elt0) + elt1;});

                            std::transform(t_mat[i0 * g1_size + i1].cbegin() + rem_pow, t_mat[i0 * g1_size + i1].cbegin() + nrows,
                                           acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
                                        return (sign ? elt0 : -elt0) + elt1;});
                        }
                    }

                    // final tau_j^{-1}
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & m_mask;
                        if (index >= poly_degree)
                            result[i1][j][index & n_mask] = -acc_vec[i];
                        else
                            result[i1][j][index] = acc_vec[i];
                        index_raw += spp.inv_Gal_u_logN_ell_[j];
                    }
                }
            }
        };

        // multi-thread
        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(mat_pp_slice, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

    // Use two uint64_t numbers to represent the 128-bit integers
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_row_combine128(const std::vector<std::vector<T>> &con_mat, const SPPParams &spp,
                               std::vector<std::vector<std::vector<uint64_t>>> &result, uint32_t threads,
                               uint32_t poly_degree = kPolyDegree, uint32_t ncols = kPolyDegree)
    {
        size_t g1_size = (size_t)1 << (spp.ell_ + spp.h_ - spp.u_);
        size_t g0_size = (size_t)1 << (spp.u_ - spp.h_);
        result.resize(g1_size);

        auto m_mask{(poly_degree << 1) - 1};
        auto n_mask{poly_degree - 1};

        auto set_to_uint128 = [](int64_t elt, uint64_t *res_128){
            res_128[0] = elt;
            res_128[1] = (elt >= 0) ? 0 : 0xFFFFFFFFFFFFFFFFULL;
        };

        auto add_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            result[0] += elt[0];
            uint8_t carry = (result[0] < elt[0]);
            result[1] += (elt[1] + carry);
        };

        auto sub_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            uint8_t borrow = (result[0] < elt[0]);
            result[0] -= elt[0];
            result[1] -= (elt[1] + borrow);
        };

        auto negate_uint128 = [](const uint64_t *elt, uint64_t *result){
            result[0] = ~elt[0] + 1;
            uint8_t carry = (result[0] == 0);
            result[1] = ~elt[1] + carry;
        };

        auto mat_pp_slice = [&](size_t bgn, size_t end){
            uint64_t temp[2];
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(g0_size);
                for (size_t j0 = 0; j0 < g0_size; ++j0){
                    result[i1][j0].resize(poly_degree << 1);
                    std::vector<std::array<uint64_t, 2>> acc_vec(poly_degree, {0, 0});
                    for (uint32_t i0 = 0; i0 < g0_size; ++i0){
                        auto inv_powx = spp.inv_ge_mul_powx_[j0][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign)
                            inv_powx &= n_mask;
                        auto rem_inv_pow = poly_degree - inv_powx;

                        if (ncols <= rem_inv_pow){
                            for (size_t j = 0; j < ncols; ++j){
                                set_to_uint128(con_mat[i0 * g1_size + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                            }
                        }else{
                            for (size_t j = 0; j < rem_inv_pow; ++j){
                                set_to_uint128(con_mat[i0 * g1_size + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                        }
                            for (size_t j = rem_inv_pow; j < ncols; ++j){
                                set_to_uint128(con_mat[i0 * g1_size + i1][j], temp);
                                if (sign)
                                    add_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                                else
                                    sub_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                            }
                        }


//                        for (size_t j = 0; j < rem_inv_pow; ++j){
//                            set_to_uint128(con_mat[i0 * g1_size + i1][j], temp);
//                            if (sign)
//                                sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
//                            else
//                                add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
//                        }
//                        for (size_t j = rem_inv_pow; j < poly_degree; ++j){
//                            set_to_uint128(con_mat[i0 * g1_size + i1][j], temp);
//                            if (sign)
//                                add_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
//                            else
//                                sub_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
//                        }
                    }

                    // final tau_j0
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & m_mask;
                        if (index >= poly_degree){
                            negate_uint128(acc_vec[i].data(), temp);
                            result[i1][j0][(index & n_mask) << 1] = temp[0];
                            result[i1][j0][((index & n_mask) << 1) + 1] = temp[1];
                        }else{
                            result[i1][j0][index << 1] = acc_vec[i][0];
                            result[i1][j0][(index << 1) + 1] = acc_vec[i][1];
                        }
                        index_raw += spp.Gal_u_h[j0];
                    }
                }
            }
        };

        // multi-thread implementation
        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(mat_pp_slice, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_col_combine128(const std::vector<const T *> &matrix, const SPPParams_Expand &spp, uint32_t ncols, uint32_t pad_cols,
                               std::vector<std::vector<std::vector<uint64_t>>> &result, uint32_t threads, uint32_t poly_degree = kPolyDegree)
    {
        assert(pad_cols >= ncols);
        if (matrix.empty())
            throw std::invalid_argument("empty matrix");

        size_t g0_size = (1ULL << (spp.u_ + spp.ell_)) / poly_degree;
        size_t g1_size = poly_degree >> spp.u_;

        // transpose and pad
        std::vector<std::vector<T>> t_mat;
        uint32_t nrows = matrix.size();
        matrix_transpose_and_padding(matrix, nrows, ncols, pad_cols, t_mat, poly_degree);

        result.resize(g1_size);

        auto m_mask{(poly_degree << 1) - 1};
        auto n_mask{poly_degree - 1};

        auto set_to_uint128 = [](int64_t elt, uint64_t *res_128){
            res_128[0] = elt;
            res_128[1] = (elt >= 0) ? 0 : 0xFFFFFFFFFFFFFFFFULL;
        };

        auto add_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            result[0] += elt[0];
            uint8_t carry = (result[0] < elt[0]);
            result[1] += (elt[1] + carry);
        };

        auto sub_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            uint8_t borrow = (result[0] < elt[0]);
            result[0] -= elt[0];
            result[1] -= (elt[1] + borrow);
        };

        auto negate_uint128 = [](const uint64_t *elt, uint64_t *result){
            result[0] = ~elt[0] + 1;
            uint8_t carry = (result[0] == 0);
            result[1] = ~elt[1] + carry;
        };

        auto mat_pp_slice = [&](size_t bgn, size_t end){
            uint64_t temp[2];
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(g0_size);
                for (size_t j = 0; j < g0_size; ++j){
                    result[i1][j].resize(poly_degree << 1);
                    std::vector<std::array<uint64_t, 2>> acc_vec(poly_degree, {0, 0});
                    for (uint32_t i0 = 0; i0 < g0_size; ++i0){
                        auto powx = spp.ge_mul_powx_[j][i0];
                        bool sign = (powx >= poly_degree);
                        if (sign) powx &= n_mask;
                        auto rem_pow = poly_degree - powx;

                        if (nrows <= rem_pow){
                            for (size_t k = 0; k < nrows; ++k){
                                set_to_uint128(t_mat[i0 * g1_size + i1][k], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[powx + k].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[powx + k].data(), temp);
                            }
                        }else{
                            for (size_t k = 0; k < rem_pow; ++k){
                                set_to_uint128(t_mat[i0 * g1_size + i1][k], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[powx + k].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[powx + k].data(), temp);
                            }

                            for (size_t k = rem_pow; k < nrows; ++k){
                                set_to_uint128(t_mat[i0 * g1_size + i1][k], temp);
                                if (sign)
                                    add_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
                                else
                                    sub_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
                            }
                        }


//                        for (size_t k = 0; k < rem_pow; ++k){
//                            set_to_uint128(t_mat[i0 * g1_size + i1][k], temp);
//                            if (sign)
//                                sub_uint128_inplace(acc_vec[powx + k].data(), temp);
//                            else
//                                add_uint128_inplace(acc_vec[powx + k].data(), temp);
//                        }
//                        for (size_t k = rem_pow; k < poly_degree; ++k){
//                            set_to_uint128(t_mat[i0 * g1_size + i1][k], temp);
//                            if (sign)
//                                add_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
//                            else
//                                sub_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
//                        }
                    }

                    // final tau_j^{-1}
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & m_mask;
                        if (index >= poly_degree){
                            negate_uint128(acc_vec[i].data(), temp);
                            result[i1][j][(index & n_mask) << 1] = temp[0];
                            result[i1][j][((index & n_mask) << 1) + 1] = temp[1];
                        }else{
                            result[i1][j][index << 1] = acc_vec[i][0];
                            result[i1][j][(index << 1) + 1] = acc_vec[i][1];
                        }
                        index_raw += spp.inv_Gal_u_logN_ell_[j];
                    }
                }
            }
        };

        // multi-thread
        uint32_t thread_block = (g1_size + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, g1_size);
            thread_pool[i] = std::thread(mat_pp_slice, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }
}

#endif // RHOMBUS_MATINFO_H
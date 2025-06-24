
#include "matrix.h"

namespace antchain::global{    

    void gen_mat_info(MatInfo &mat_info, uint32_t rows, uint32_t cols, uint32_t mat_bits, uint32_t poly_degree)
    {
        if (rows == 0 || cols == 0 || rows > kPolyDegree || cols > kPolyDegree)
            throw std::invalid_argument("Invalid rows or columns");
        if (mat_bits > 64)
            throw std::invalid_argument("Matrix bits too big");

        mat_info.nrows_ = rows;
        mat_info.ncols_ = cols;

        mat_info.mat_bits_ = mat_bits;

        mat_info.log_pad_nrows_ = BitLength(rows - 1);
        mat_info.log_pad_ncols_ = BitLength(cols - 1);

        // zero padding to ensure r * c >= N
        if (mat_info.log_pad_ncols_ + mat_info.log_pad_nrows_ < kLogPolyDegree)
            mat_info.log_pad_ncols_ = kLogPolyDegree - mat_info.log_pad_nrows_;

        gen_spp_parms(mat_info.spp_parms_, mat_info.log_pad_nrows_,
                      mat_info.log_pad_ncols_, poly_degree, kSPPMap);
    }

    void gen_lw_mat_info(MatInfoLW &mat_info, uint32_t rows, uint32_t cols,
                         uint32_t mat_bits, uint32_t poly_degree)
    {
        if (rows == 0 || cols == 0 || rows > kPolyDegree || cols > kPolyDegree)
            throw std::invalid_argument("Invalid rows or columns");
        if (mat_bits > 64)
            throw std::invalid_argument("Matrix bits too big");

        mat_info.nrows_ = rows;
        mat_info.ncols_ = cols;

        mat_info.mat_bits_ = mat_bits;

        // Note: we don't require nrows * ncols >= N
        mat_info.log_pad_nrows_ = BitLength(rows - 1);
        mat_info.log_pad_ncols_ = BitLength(cols - 1);
    }

    void gen_large_mat_info(LargeMatInfo &mat_info, uint32_t rows, uint32_t cols, uint32_t mat_bits, uint32_t poly_degree)
    {
        if (rows == 0 || cols == 0)
            throw std::invalid_argument("Invalid matrix rows or columns");

        mat_info.nrows_ = rows;
        mat_info.ncols_ = cols;
        mat_info.mat_bits_ = mat_bits;

        mat_info.nrow_block_num_ = (rows + kPolyDegree - 1) >> kLogPolyDegree;
        mat_info.ncol_block_num_ = (cols + kPolyDegree - 1) >> kLogPolyDegree;

        // block[0][0]
        uint32_t block_nrows = std::min(kPolyDegree, rows);
        uint32_t block_ncols = std::min(kPolyDegree, cols);
        gen_mat_info(mat_info.corner_block_array_[0], block_nrows, block_ncols, mat_bits, poly_degree);

        // block[0][last_col]
        // only one block in column
        if (mat_info.ncol_block_num_ == 1)
            mat_info.corner_block_array_[1] = mat_info.corner_block_array_[0];
        else{
            // block_rows equals to block[0][0]
            block_ncols = cols - ((mat_info.ncol_block_num_ - 1) << kLogPolyDegree);
            gen_mat_info(mat_info.corner_block_array_[1], block_nrows, block_ncols, mat_bits, poly_degree);
        }

        // block[last_row][0]
        if (mat_info.nrow_block_num_ == 1)
            mat_info.corner_block_array_[2] = mat_info.corner_block_array_[0];
        else{
            block_nrows = rows - ((mat_info.nrow_block_num_ - 1) << kLogPolyDegree);
            block_ncols = mat_info.corner_block_array_[0].ncols_;
            gen_mat_info(mat_info.corner_block_array_[2], block_nrows, block_ncols, mat_bits, poly_degree);
        }

        // block[last_row][last_col]
        if (mat_info.nrow_block_num_ == 1){
            mat_info.corner_block_array_[3] = mat_info.corner_block_array_[1];
            return;
        }
        if (mat_info.ncol_block_num_ == 1){
            mat_info.corner_block_array_[3] = mat_info.corner_block_array_[2];
            return;
        }
        block_nrows = rows - ((mat_info.nrow_block_num_ - 1) << kLogPolyDegree);
        block_ncols = cols - ((mat_info.ncol_block_num_ - 1) << kLogPolyDegree);
        gen_mat_info(mat_info.corner_block_array_[3], block_nrows, block_ncols, mat_bits, poly_degree);
    }

    void gen_spp_parms(SPPParams &spp_parms, uint32_t log_pad_rows, uint32_t log_pad_cols,
                       uint32_t poly_degree, const int *SPP_table)
    {
        uint32_t logN = BitLength(poly_degree - 1);
        uint32_t ell = log_pad_rows + log_pad_cols - logN;
        spp_parms.h_ = logN - log_pad_cols;
        spp_parms.ell_ = ell;
        spp_parms.u_ = spp_parms.h_ + SPP_table[spp_parms.ell_];

        spp_parms.PackRLWEs_factor_ = (uint32_t)1 << spp_parms.ell_;

        size_t g0_size = size_t(1) << (spp_parms.u_ - spp_parms.h_);
        spp_parms.Gal_u_h.resize(g0_size);
        std::vector<uint32_t> inv_Gal_u_h(g0_size);

        // galois elements in Gal(K_u/K_h): combination ~ 2^{h+1}+1, 2^{h+2}+1, ..., 2^{u}+1
        // init
        spp_parms.Gal_u_h[0] = 1;
        inv_Gal_u_h[0] = 1;

        uint32_t iter = spp_parms.h_ + 1;
        uint32_t M = poly_degree << 1; // 2N
        uint32_t M_mask = M - 1;

        for (uint32_t i = 0; i < (spp_parms.u_ - spp_parms.h_); ++i){
            for (uint32_t j = 0; j < ((uint32_t)1 << i); ++j){
                size_t index = (1 << i) + j;
                spp_parms.Gal_u_h[index] = (((1 << iter) + 1) * spp_parms.Gal_u_h[j]) & M_mask;
                inv_Gal_u_h[index] = ModInv(spp_parms.Gal_u_h[index], M);
            }
            ++iter;
        }

        // power of X: [0, N/2^u, 2*(N/2^u), 3*(N/2^u), ...]
        std::vector<uint32_t> powx(g0_size);
        uint32_t base_pow = poly_degree >> spp_parms.u_; // N/2^u
        for (uint32_t i = 0; i < g0_size; ++i)
            powx[i] = i * base_pow;

        // tau^{-1}(X^k) for tau in Gal(K_u/K_h), k in powx
        // g0_size * g0_size table
        spp_parms.inv_ge_mul_powx_.resize(g0_size);
        for (size_t i = 0; i < g0_size; ++i){
            spp_parms.inv_ge_mul_powx_[i].resize(g0_size);
            std::transform(powx.cbegin(), powx.cend(), spp_parms.inv_ge_mul_powx_[i].begin(),
                           [&](uint32_t elt){return (elt * inv_Gal_u_h[i]) & M_mask;});
        }
    }

    void gen_spp_expand(SPPParams_Expand &spp_expand, uint32_t log_pad_cols,
                        uint32_t poly_degree, const int *SPP_table)
    {
        // Tr_{logN, logN - ell} = Tr_{u, logN - ell} * Tr_{logN, u}

        uint32_t logN = BitLength(poly_degree - 1);
        uint32_t ell = log_pad_cols;
        spp_expand.ell_ = ell;
        spp_expand.u_ = logN - ell + SPP_table[ell];
        spp_expand.Expand_factor_ = (uint32_t)1 << ell;

        size_t g0_size = size_t(1) << (spp_expand.u_ - logN + ell);
        spp_expand.Gal_u_logN_ell_.resize(g0_size);
        spp_expand.inv_Gal_u_logN_ell_.resize(g0_size);

        // init
        // Gal(Ku/K_{logN-ell})
        spp_expand.Gal_u_logN_ell_[0] = 1;
        spp_expand.inv_Gal_u_logN_ell_[0] = 1;

        uint32_t iter = logN - ell + 1;
        uint32_t M = poly_degree << 1; // 2N
        uint32_t M_mask = M - 1;

        for (uint32_t i = 0; i < (spp_expand.u_ - logN + ell); ++i){
            for (uint32_t j = 0; j < ((uint32_t)1 << i); ++j){
                size_t index = (1 << i) + j;
                spp_expand.Gal_u_logN_ell_[index] = (((1 << iter) + 1) * spp_expand.Gal_u_logN_ell_[j]) & M_mask;
                spp_expand.inv_Gal_u_logN_ell_[index] = ModInv(spp_expand.Gal_u_logN_ell_[index], M);
            }
            ++iter;
        }

        // power of X: [0, -(N/2^u), -2(N/2^u), -3(N/2^u) ... ]
        std::vector<uint32_t> powx(g0_size);
        uint32_t base_pow = poly_degree >> spp_expand.u_;
        for (uint32_t i = 0; i < g0_size; ++i){
            powx[i] = (M - i * base_pow);
        }

        // tau_i(X^{-(N/2^u) * j})
        spp_expand.ge_mul_powx_.resize(g0_size);
        for (size_t i = 0; i < g0_size; ++i){
            spp_expand.ge_mul_powx_[i].resize(g0_size);
            std::transform(powx.cbegin(), powx.cend(), spp_expand.ge_mul_powx_[i].begin(),
                           [&](uint32_t elt){return (elt * spp_expand.Gal_u_logN_ell_[i]) & M_mask;});
        }
    }

}


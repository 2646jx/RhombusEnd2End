#ifndef RHOMBUS_RING_OP_H
#define RHOMBUS_RING_OP_H

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>

namespace antchain{

    // generate random vector, all elements are in [-2^(bits-1), 2^(bits-1))
    // the type T should be int8_t, int16_t, int32_t or int64_t.
    template <typename T, typename std::enable_if<std::is_signed<T>::value, T>::type * = nullptr>
    void GenRandIntVector(T *vec, size_t vec_len, uint32_t bits)
    {
        if (vec == nullptr || vec_len == 0)
            throw std::invalid_argument("empty vector is not allowed");
        if (bits > sizeof(T) * 8)
            throw std::invalid_argument("bits out of range");

        // randomly
        if (bits == sizeof(T) * 8)
        {
            auto gen_rand_vec = [&](size_t bgn, size_t end){
                std::mt19937 gen(std::random_device{}());
                if (sizeof(T) == 8){ // int64_t
                    std::uniform_int_distribution<T> dist(INT64_MIN, INT64_MAX);
                    std::generate_n(vec + bgn, end - bgn, [&](){return dist(gen);});
                }else if (sizeof(T) == 4) { // int32_t or int
                    std::uniform_int_distribution<T> dist(INT32_MIN, INT32_MAX);
                    std::generate_n(vec + bgn, end - bgn, [&](){return dist(gen);});
                }else if (sizeof(T) == 2) { // int16_t
                    std::uniform_int_distribution<T> dist(INT16_MIN, INT16_MAX);
                    std::generate_n(vec + bgn, end - bgn, [&](){return dist(gen);});
                }
                else{ // int8_t
                    std::uniform_int_distribution<T> dist(INT8_MIN, INT8_MAX);
                    std::generate_n(vec + bgn, end - bgn, [&](){return dist(gen);});
                }
            };
            gen_rand_vec(0, vec_len);
        }
        else
        {
            auto gen_rand_vec = [&](size_t bgn, size_t end){
                std::mt19937 gen(std::random_device{}());
                std::uniform_int_distribution<T> dist(-(T(1) << (T)(bits - 1)), (T(1) << (T)(bits - 1)) - (T)1);
                std::generate(vec + bgn, vec + end, [&](){return dist(gen);});
            };
            gen_rand_vec(0, vec_len);
        }
    }

    template <typename T, typename std::enable_if<std::is_floating_point_v<T>, T>::type * = nullptr>
    void GenRandRealVector(T *vec, size_t vec_len, T lower_bound, T upper_bound)
    {
        if (vec == nullptr || vec_len == 0)
            throw std::invalid_argument("empty vector is not allowed");
        if (lower_bound >= upper_bound)
            throw std::invalid_argument("invalid range");

        auto gen_rand_vec = [&](size_t bgn, size_t end){
            std::mt19937 gen(std::random_device{}());
            std::uniform_real_distribution<T> dist(lower_bound, upper_bound);
            std::generate_n(vec + bgn, end - bgn, [&](){return dist(gen);});
        };
        gen_rand_vec(0, vec_len);
    }

    template <typename T, typename std::enable_if<std::is_signed<T>::value, T>::type * = nullptr>
    void GenRandIntMat(std::vector<std::vector<T>> &mat, size_t nrows, size_t ncols, uint32_t bits)
    {
        if (nrows == 0 || ncols == 0)
            throw std::invalid_argument("invalid matrix dimension");

        // resize
        mat.resize(nrows);
        for (size_t i = 0; i < nrows; ++i)
        {
            mat[i].resize(ncols);
            GenRandIntVector(mat[i].data(), ncols, bits);
        }
    }


    // generate unsigned random vector, [0, 2^bits)
    template <typename T, typename std::enable_if<std::is_unsigned<T>::value, T>::type * = nullptr>
    void GenRandUintVector(T *vec, size_t vec_len, uint32_t bits){
        if (vec == nullptr || vec_len == 0)
            throw std::invalid_argument("empty vector is not allowed");
        if (bits > sizeof(T) * 8)
            throw std::invalid_argument("bits out of range");

        auto gen_rand_vec = [&](size_t bgn, size_t end){
            std::mt19937 gen(std::random_device{}());
            T upbd = (bits == sizeof(T) * 8) ? static_cast<T>(-1) : ((T)1 << (T)bits) - (T)1;
            std::uniform_int_distribution<T> dist(0, upbd);
            std::generate(vec + bgn, vec + end, [&](){return dist(gen);});
        };
        gen_rand_vec(0, vec_len);
    }

    template <typename T, typename std::enable_if<std::is_unsigned<T>::value, T>::type * = nullptr>
    void GenRandUintMat(std::vector<std::vector<T>> &mat, size_t nrows, size_t ncols, uint32_t bits){
        if (nrows == 0 || ncols == 1)
            throw std::invalid_argument("invalid matrix dimension");

        // resize
        mat.resize(nrows);
        for (size_t i = 0; i < nrows; ++i){
            mat[i].resize(ncols);
            GenRandUintVector(mat[i].data(), ncols, bits);
        }
    }


    template <typename T, typename U,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U>, int>::type = 0>
    void IntVecToUint(const T *signed_vec, U *unsigned_vec, size_t vec_len, uint32_t mod_bits)
    {
        if (sizeof(T) != sizeof(U))
            throw std::invalid_argument("sizeof(T) should equal to sizeof(U) to avoid the precision loss");
        if (mod_bits > sizeof(U) * 8)
            throw std::invalid_argument("The sharing(mod) bits are out of range");
        U mask = (mod_bits == sizeof(U) * 8) ? static_cast<U>(-1) : ((U)1 << (U)mod_bits) - (U)1;
        std::transform(signed_vec, signed_vec + vec_len, unsigned_vec, [&](auto elt){
            return (U)elt & mask;
        });
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U>, int>::type = 0>
    void IntMatToUint(const std::vector<std::vector<T>> &signed_mat, std::vector<std::vector<U>> &unsigned_mat,
                      size_t nrows, size_t ncols, uint32_t mod_bits)
    {
        if (signed_mat.size() != nrows)
            throw std::invalid_argument("mismatch");

        // resize
        unsigned_mat.resize(nrows);
        for (size_t i = 0; i < nrows; ++i){
            unsigned_mat[i].resize(ncols);
            IntVecToUint(signed_mat[i].data(), unsigned_mat[i].data(), ncols, mod_bits);
        }
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U>, int>::type = 0>
    void UintVecToInt(const T *unsigned_vec, U *signed_vec, size_t vec_len, uint32_t mod_bits){
        if (sizeof(T) != sizeof(U))
            throw std::invalid_argument("sizeof(T) should equal to sizeof(U) to avoid the precision loss");
        if (mod_bits > sizeof(U) * 8)
            throw std::invalid_argument("The sharing(mod) bits are out of range");
        if (mod_bits == sizeof(T) * 8){
            std::transform(unsigned_vec, unsigned_vec + vec_len, signed_vec, [&](auto elt){return (U)elt;});
        }else{
            T mask = (T(1) << (T)mod_bits) - (T)1;
            std::transform(unsigned_vec, unsigned_vec + vec_len, signed_vec, [&](auto elt)->U{
                T half_mod = T(1) << T(mod_bits - 1);
                if (elt >= half_mod) return -(U)((-elt) & mask);
                else return (U)elt;
            });
        }
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_unsigned_v<T> && std::is_signed_v<U>, int>::type = 0>
    void UintMatToInt(const std::vector<std::vector<T>> &unsigned_mat, std::vector<std::vector<U>> &signed_mat,
                      size_t nrows, size_t ncols, uint32_t mod_bits)
    {
        if (unsigned_mat.size() != nrows)
            throw std::invalid_argument("dimension mismatch");

        // resize
        signed_mat.resize(nrows);
        for (size_t i = 0; i < nrows; ++i){
            signed_mat[i].resize(ncols);
            UintVecToInt(unsigned_mat[i].data(), signed_mat[i].data(), ncols, mod_bits);
        }
    }


    template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
    void AddVecMod(const T *vec_0, const T *vec_1, T *vec_result, size_t vec_len, uint32_t mod_bits)
    {
        if (vec_0 == nullptr || vec_1 == nullptr || vec_result == nullptr || vec_len == 0)
            throw std::invalid_argument("invalid vector size or nullptr");
        if (mod_bits > sizeof(T) * 8)
            throw std::invalid_argument("mod bits are too large");

        // First, convert to unsigned type, then add (with mod), and convert back if necessary
        using unsigned_T = typename std::make_unsigned_t<T>;
        unsigned_T mask = (mod_bits == sizeof(T) * 8) ? static_cast<unsigned_T>(-1)
                                                      : ((unsigned_T)1 << (unsigned_T)mod_bits) -  (unsigned_T)1;
        // add directly if the vector is already unsigned.
        if (std::is_unsigned_v<T>){
            std::transform(vec_0, vec_0 + vec_len, vec_1, vec_result, [&](auto elt0, auto elt1){
                return (elt0 + elt1) & mask;});
        }else{
            std::vector<unsigned_T> temp(vec_len);
            std::transform(vec_0, vec_0 + vec_len, vec_1, temp.data(), [&](auto elt0, auto elt1){
                return (unsigned_T)(elt0 + elt1) & mask;});
            UintVecToInt(temp.data(), vec_result, vec_len, mod_bits);
        }
    }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
    void AddMatMod(const std::vector<std::vector<T>> &mat_0, const std::vector<std::vector<T>> mat_1,
                   std::vector<std::vector<T>> &mat_result, size_t nrows, size_t ncols, uint32_t mod_bits)
    {
        // check
        if (mat_0.size() != nrows || mat_1.size() != nrows)
            throw std::invalid_argument("dimension mismatch");

        // resize
        mat_result.resize(nrows);
        for (size_t i = 0; i < nrows; ++i){
            mat_result[i].resize(ncols);
            AddVecMod(mat_0[i].data(), mat_1[i].data(), mat_result[i].data(), ncols, mod_bits);
        }
    }

    template <typename T, typename U, typename V,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U> && std::is_integral_v<V>, int>::type = 0>
    void MatVecMulMod(const T *mat, const U *vec, V *result, size_t nrows, size_t ncols, uint32_t mod_bits){
        using unsigned_V = typename std::make_unsigned_t<V>;
        unsigned_V mask = (mod_bits == sizeof(V) * 8) ? static_cast<unsigned_V>(-1)
                                                      : ((unsigned_V)1 << (unsigned_V)mod_bits) - (unsigned_V)1;

        std::vector<unsigned_V> temp_vec(nrows);
        for (size_t i = 0; i < nrows; ++i){
            unsigned_V temp = 0;
            for (size_t j = 0; j < ncols; ++j){
                temp += (unsigned_V)mat[i * ncols + j] * (unsigned_V)vec[j];
            }
            temp_vec[i] = temp & mask;
        }
        if (std::is_unsigned_v<V>){
            std::copy_n(temp_vec.data(), nrows, result);
        }else{
            UintVecToInt(temp_vec.data(), result, nrows, mod_bits);
        }
    }

    template <typename T, typename U, typename V,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U> && std::is_integral_v<V>, int>::type = 0>
    void MatVecMulMod(const std::vector<std::vector<T>> &mat, const std::vector<U> &vec,
                      std::vector<V> &result, size_t nrows, size_t ncols, uint32_t mod_bits){
        using unsigned_V = typename std::make_unsigned_t<V>;
        unsigned_V mask = (mod_bits == sizeof(V) * 8) ? static_cast<unsigned_V>(-1)
                                                      : ((unsigned_V)1 << (unsigned_V)mod_bits) - (unsigned_V)1;

        std::vector<unsigned_V> temp_vec(nrows);
        for (size_t i = 0; i < nrows; ++i){
            unsigned_V temp = 0;
            for (size_t j = 0; j < ncols; ++j){
                temp += (unsigned_V)mat[i][j] * (unsigned_V)vec[j];
            }
            temp_vec[i] = temp & mask;
        }

        // resize
        result.resize(ncols);

        if (std::is_unsigned_v<V>){
            std::copy_n(temp_vec.data(), nrows, result.data());
        }else{
            UintVecToInt(temp_vec.data(), result.data(), nrows, mod_bits);
        }
    }

    // X: n * m, Y: m * k
    template <typename T, typename U, typename V,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U> && std::is_integral_v<V>, int>::type = 0>
    void MatMulMod(const T *matX, const U *matY, V *result, size_t n, size_t m, size_t k, uint32_t mod_bits){
        using unsigned_V = typename std::make_unsigned_t<V>;
        unsigned_V mask = (mod_bits == sizeof(V) * 8) ? static_cast<unsigned_V>(-1)
                                                      : ((unsigned_V)1 << (unsigned_V)mod_bits) - (unsigned_V)1;

        std::vector<unsigned_V> temp_mat(n * k);
        for (size_t i = 0 ; i < n; ++i){
            for (size_t j = 0; j < k; ++j){
                unsigned_V temp = 0;
                for (size_t l = 0; l < m; ++l){
                    temp += (unsigned_V)matX[i * m + l] * (unsigned_V)matY[l * k + j];
                }
                temp_mat[i * k + j] = temp & mask;
            }
        }
        if (std::is_unsigned_v<V>){
            std::copy_n(temp_mat.data(), n * k, result);
        }else{
            UintVecToInt(temp_mat.data(), result, n * k, mod_bits);
        }
    }

    template <typename T, typename U, typename V,
            typename std::enable_if<std::is_integral_v<T> && std::is_integral_v<U> && std::is_integral_v<V>, int>::type = 0>
    void MatMulMod(const std::vector<std::vector<T>> &matX, const std::vector<std::vector<U>> &matY,
                   std::vector<std::vector<V>> &result, size_t n, size_t m, size_t k, uint32_t mod_bits)
    {
        using unsigned_V = typename std::make_unsigned_t<V>;
        unsigned_V mask = (mod_bits == sizeof(V) * 8) ? static_cast<unsigned_V>(-1)
                                                      : ((unsigned_V)1 << (unsigned_V)mod_bits) - (unsigned_V)1;

        std::vector<std::vector<unsigned_V>> temp_mat(n);
        for (size_t i = 0 ; i < n; ++i){
            temp_mat[i].resize(k);
            for (size_t j = 0; j < k; ++j){
                unsigned_V temp = 0;
                for (size_t l = 0; l < m; ++l){
                    temp += (unsigned_V)matX[i][l] * (unsigned_V)matY[l][j];
                }
                temp_mat[i][j] = temp & mask;
            }
        }

        // resize
        result.resize(n);
        if (std::is_unsigned_v<V>){
            for (size_t i = 0; i < n; ++i){
                result[i].resize(k);
                std::copy_n(temp_mat[i].data(), k, result[i].data());
            }
        }else{
            for (size_t i = 0; i < n; ++i){
                result[i].resize(k);
                UintVecToInt(temp_mat[i].data(), result[i].data(), k, mod_bits);
            }
        }
    }

    // find max diff over the ring Z_{2^mod_bits}
    template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
    size_t CompVector(const T *vec0, const T *vec1, size_t vec_len, uint64_t &max_diff, uint32_t mod_bits)
    {
        if (mod_bits > sizeof(T) * 8)
            throw std::invalid_argument("mod bits too large");
        max_diff = 0;
        using unsigned_T = typename std::make_unsigned_t<T>;
        unsigned_T cur_max_diff = 0;
        unsigned_T half_th = (unsigned_T)1 << (mod_bits - 1);
        size_t max_diff_index = 0;
        unsigned_T mask = (mod_bits == sizeof(T) * 8) ? static_cast<unsigned_T>(-1) : (unsigned_T(1) << mod_bits) - unsigned_T(1);
        for (size_t i = 0; i < vec_len; ++i)
        {
            cur_max_diff = unsigned_T(vec0[i] - vec1[i]) & mask;
            if (cur_max_diff > half_th){
                if (mod_bits == sizeof(T) * 8){
                    cur_max_diff = -cur_max_diff;
                }else{
                    cur_max_diff = (unsigned_T(1) << mod_bits) - cur_max_diff;
                }
            }
            if (cur_max_diff > max_diff){
                max_diff = cur_max_diff;
                max_diff_index = i;
            }
        }
        return max_diff_index;
    }

    template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
    std::tuple<size_t, size_t> CompMatrix(const std::vector<std::vector<T>> &mat0, const std::vector<std::vector<T>> &mat1,
                                          size_t nrows, size_t ncols, uint64_t &max_diff, uint32_t mod_bits)
    {
        if (mod_bits > sizeof(T) * 8)
            throw std::invalid_argument("mod bits too large");
        if (mat0.size() != nrows || mat1.size() != nrows)
            throw std::invalid_argument("dimension mismatch");

        max_diff = 0;
        using unsigned_T = typename std::make_unsigned_t<T>;
        unsigned_T cur_max_diff = 0;
        unsigned_T half_th = (unsigned_T)1 << (mod_bits - 1);
        size_t max_diff_row_index = 0;
        size_t max_diff_col_index = 0;
        unsigned_T mask = (mod_bits == sizeof(T) * 8) ? static_cast<unsigned_T>(-1) : (unsigned_T(1) << mod_bits) - unsigned_T(1);

        for (size_t i = 0; i < nrows; ++i){
            for (size_t j = 0; j < ncols; ++j){
                cur_max_diff = unsigned_T(mat0[i][j] - mat1[i][j]) & mask;
                if (cur_max_diff > half_th){
                    if (mod_bits == sizeof(T) * 8){
                        cur_max_diff = -cur_max_diff;
                    }else{
                        cur_max_diff = (unsigned_T(1) << mod_bits) - cur_max_diff;
                    }
                }
                if (cur_max_diff > max_diff){
                    max_diff = cur_max_diff;
                    max_diff_row_index = i;
                    max_diff_col_index = j;
                }
            }
        }
        return std::make_tuple(max_diff_row_index, max_diff_col_index);
    }

    // shift over coefficients mod X^N+1.
    // Input: coefficient vector of a(X)
    // Output: coefficient vector of a(X) * X^{exp} mod (X^N+1)
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void shift_vector(const std::vector<T> &vec, std::vector<T> &output, uint32_t exp, uint32_t poly_degree)
    {
        if (exp >= 2 * poly_degree)
            throw std::invalid_argument("exp should be in [0, 2N)");
        bool sign = (exp >= poly_degree);
        uint32_t n_mask = poly_degree - 1;
        uint32_t pow = exp & n_mask;

        output.resize(poly_degree);
        if (sign){
            std::transform(vec.cbegin(), vec.cbegin() + poly_degree - pow, output.begin() + pow, [&](auto elt){
                return -elt;});
            std::copy_n(vec.cbegin() + poly_degree - pow, pow, output.begin());
        }else{
            std::copy_n(vec.cbegin(), poly_degree - pow, output.begin() + pow);
            std::transform(vec.cbegin() + poly_degree - pow, vec.cend(), output.begin(), [&](auto elt){
                return -elt;});
        }
    }

    // Input: coefficient vector of a(X)
    // Output: coefficient vector of a(X^ge) mod (X^N+1)
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void aut_permute(const std::vector<T> &vec, std::vector<T> &output, uint32_t ge, uint32_t poly_degree)
    {
        uint32_t m_mask = (poly_degree << 1) - 1;
        uint32_t n_mask = poly_degree - 1;

        output.resize(poly_degree);
        uint32_t index_raw = 0;
        for (uint32_t i = 0; i < poly_degree; ++i){
            uint32_t index = index_raw & m_mask;
            if (index >= poly_degree)
                output[index & n_mask] = -vec[i];
            else
                output[index] = vec[i];
            index_raw += ge;
        }
    }


    template <typename T, typename std::enable_if<std::is_unsigned_v<T>, T>::type * = nullptr>
    uint32_t BitLength(T n){
        if (n == 0)
            return 0;
        uint32_t count = 0;
        while (n){
            n >>= T(1);
            ++count;
        }
        return count;
    }


    // a^{-1} in Z_{mod}
    // If a don't have an inverse, return 0.
    template <typename T, typename std::enable_if<std::is_unsigned_v<T>, T>::type * = nullptr>
    T ModInv(T a, T mod){
        std::vector<std::vector<T>> Euc_table(4);
        for (size_t i = 0; i < 4; ++i)
            Euc_table[i].resize(2);
        Euc_table[0][0] = mod;
        Euc_table[0][1] = a;
        Euc_table[1][1] = mod / a;
        Euc_table[2][0] = T(1);
        Euc_table[2][1] = 0;
        Euc_table[3][0] = 0;
        Euc_table[3][1] = T(1);

        size_t i = 1;
        while (Euc_table[0][i] != (T)1){
            ++i;
            Euc_table[0].push_back(Euc_table[0][i-2] - Euc_table[0][i-1] * Euc_table[1][i-1]);
            Euc_table[2].push_back(Euc_table[2][i-2] - Euc_table[2][i-1] * Euc_table[1][i-1]);
            Euc_table[3].push_back(Euc_table[3][i-2] - Euc_table[3][i-1] * Euc_table[1][i-1]);
            Euc_table[1].push_back(Euc_table[0][i-1] / Euc_table[0][i]);
            if (Euc_table[0][i] == 0)
                return 0;
        }
        return Euc_table[3][i] % mod;
    }


    template <typename T, typename std::enable_if<std::is_arithmetic_v<T>, T>::type * = nullptr>
    void print_vector(const T * vec, size_t vec_len, size_t output_num = 3)
    {
        std::cout << "[";
        if (vec_len <= 2 * output_num){
            for (size_t i = 0; i < vec_len - 1; ++i){
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec_len - 1];
        }else{
            for (size_t i = 0; i < output_num; ++i)
                std::cout << vec[i] << ", ";
            std::cout << "..., ";
            for (size_t i = vec_len - output_num; i < vec_len - 1; ++i)
                std::cout << vec[i] << ", ";
            std::cout << vec[vec_len - 1];
        }
        std::cout << "]" << std::endl;
    }

    // print first output_rows rows of mat
    template <typename T, typename std::enable_if<std::is_arithmetic_v<T>, T>::type * = nullptr>
    void print_matrix(const T *mat, size_t rows, size_t cols, size_t output_rows = 3, size_t output_cols = 3)
    {
        if (rows <= output_rows + 2){
            for (size_t i = 0; i < rows; ++i){
                print_vector(mat + i * cols, cols, output_cols);
            }
        }else{
            for (size_t i = 0; i < output_rows; ++i){
                print_vector(mat + i * cols, cols, output_cols);
            }
            std::cout << std::endl;
        }
    }


    template <typename T, typename std::enable_if<std::is_arithmetic_v<T>, T>::type * = nullptr>
    void print_matrix(const std::vector<const T*> &mat, size_t cols, size_t output_rows = 3, size_t output_cols = 3)
    {
        size_t rows = mat.size();
        if (rows == 0)
            throw std::invalid_argument("empty matrix");
        if (cols == 0)
            throw std::invalid_argument("empty matrix");
        if (rows <= output_rows + 2){
            for (size_t i = 0; i < rows; ++i){
                print_vector(mat[i], cols, output_cols);
            }
        }else{
            for (size_t i = 0; i < output_rows; ++i){
                print_vector(mat[i], cols, output_cols);
            }
            std::cout << std::endl;
        }
    }

    template <typename T, typename std::enable_if<std::is_arithmetic_v<T>, T>::type * = nullptr>
    void print_matrix(const std::vector<std::vector<T>> &mat, size_t output_rows = 3, size_t output_cols = 3)
    {
        size_t rows = mat.size();
        if (rows == 0)
            throw std::invalid_argument("empty matrix");
        size_t cols = mat[0].size();
        if (cols == 0)
            throw std::invalid_argument("empty matrix");

        if (rows <= output_rows + 2){
            for (size_t i = 0; i < rows; ++i){
                print_vector(mat[i].data(), cols, output_cols);
            }
        }else{
            for (size_t i = 0; i < output_rows; ++i){
                print_vector(mat[i].data(), cols, output_cols);
            }
            std::cout << std::endl;
        }
    }

    template <typename T, typename std::enable_if<std::is_arithmetic_v<T>, T>::type * = nullptr>
    void print_sub_matrix(const T *mat, size_t rows, size_t cols,
                          size_t row_bgn, size_t row_num, size_t col_bgn, size_t col_num,
                          size_t output_rows = 3, size_t output_cols = 3)
    {
        if (row_bgn + row_num > rows || col_bgn + col_num > cols)
            throw std::invalid_argument("out of index");

        // construct sub-matrix
        std::vector<const T *> sub_mat(row_num);
        for (size_t i = 0; i < row_num; ++i){
            sub_mat[i] = mat + (row_bgn + i) * cols + col_bgn;
        }

        print_matrix(sub_mat, col_num, output_rows, output_cols);
    }

}

#endif 
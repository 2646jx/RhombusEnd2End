
#include "seal_api.h"
#include "seal/ciphertext.h"
#include "seal/galoiskeys.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/uintarithmod.h"
#include "seal/util/rlwe.h"
#include <thread>

namespace antchain
{

    void gen_aux_params(AuxParms &aux_parms, uint32_t mod_bits, const std::vector<uint64_t> &AllMods, uint32_t kNumMods)
    {

        if (mod_bits < 1 || mod_bits > 64)
            throw std::invalid_argument("the modulus bits should be in [1, 64]");
        aux_parms.pow_2_t_mod_qi_.resize(kNumMods - 1);
        aux_parms.inv_2_t_mod_qi_.resize(kNumMods - 1);
        aux_parms.inv_qi_mod_2_t_.resize(kNumMods - 1);
        aux_parms.pow_2_64_mod_qi_.resize(kNumMods - 1);
        if (mod_bits == 64){
            for (size_t i = 0; i < kNumMods - 1; ++i)
                aux_parms.pow_2_t_mod_qi_[i] = seal::util::multiply_uint_mod((1ULL << 32), (1ULL << 32), AllMods[i]);
            uint64_t pow_2_t[]{0, 1};
            for (size_t i = 0; i < kNumMods - 1; ++i){
                aux_parms.pow_2_64_mod_qi_[i] = seal::util::multiply_uint_mod((1ULL << 32), (1ULL << 32), AllMods[i]);
                seal::util::try_invert_uint_mod(aux_parms.pow_2_t_mod_qi_[i], AllMods[i], aux_parms.inv_2_t_mod_qi_[i]);
                uint64_t vec_qi[]{AllMods[i], 0};
                uint64_t vec_inv_qi_mod_2_t[]{0, 0};

                seal::util::try_invert_uint_mod(vec_qi, pow_2_t, 2, vec_inv_qi_mod_2_t, seal::MemoryManager::GetPool());
                aux_parms.inv_qi_mod_2_t_[i] = vec_inv_qi_mod_2_t[0];
            }
        }else{
            uint64_t pow_2_t{1ULL << mod_bits};
            for (size_t i = 0; i < kNumMods - 1; ++i)
                aux_parms.pow_2_t_mod_qi_[i] = seal::util::barrett_reduce_64(pow_2_t, AllMods[i]);

            for (size_t i = 0; i < kNumMods - 1; ++i){
                aux_parms.pow_2_64_mod_qi_[i] = seal::util::multiply_uint_mod((1ULL << 32), (1ULL << 32), AllMods[i]);
                seal::util::try_invert_uint_mod(aux_parms.pow_2_t_mod_qi_[i], AllMods[i], aux_parms.inv_2_t_mod_qi_[i]);
                seal::util::try_invert_uint_mod(AllMods[i], pow_2_t, aux_parms.inv_qi_mod_2_t_[i]);
            }
        }

        if (kNumMods >= 3)
            seal::util::try_invert_uint_mod(AllMods[1], AllMods[0], aux_parms.inv_q1_mod_q0_);
        if (kNumMods == 4){
            seal::util::try_invert_uint_mod(AllMods[2], AllMods[0], aux_parms.inv_q2_mod_q0_);
            seal::util::try_invert_uint_mod(AllMods[2], AllMods[1], aux_parms.inv_q2_mod_q1_);
        }
    }

    static void divide_and_round_q_last_ntt_inplace(seal::util::RNSIter input, size_t input_base_size,
                                                    seal::util::ConstNTTTablesIter rns_ntt_tables, const seal::util::RNSTool *rns_tool, seal::MemoryPoolHandle pool)
    {
        using namespace seal;
        using namespace seal::util;
        if (!input or !rns_tool or !rns_ntt_tables or !pool)
        {
            throw std::invalid_argument("null ptr");
        }

        const size_t coeff_count = input.poly_modulus_degree();
        auto base_q = rns_tool->base_q();
        const size_t base_q_size = base_q->size();
        if (input_base_size >= base_q_size || input_base_size < 1)
        {
            throw std::invalid_argument("invalid size");
        }

        CoeffIter last_input = input[input_base_size];
        // Convert to non-NTT form
        inverse_ntt_negacyclic_harvey(last_input, rns_ntt_tables[base_q_size - 1]);
        // Add (qi-1)/2 to change from flooring to rounding
        const Modulus &last_modulus = (*base_q)[base_q_size - 1];
        uint64_t half = last_modulus.value() >> 1;
        add_poly_scalar_coeffmod(last_input, coeff_count, half, last_modulus,
                                 last_input);

        SEAL_ALLOCATE_GET_COEFF_ITER(temp, coeff_count, pool);
        SEAL_ITERATE(
            iter(input, rns_tool->inv_q_last_mod_q(), base_q->base(), rns_ntt_tables),
            input_base_size, [&](auto I)
            {
                    // (ct mod qk) mod qi
                    if (get<2>(I).value() < last_modulus.value()) {
                        modulo_poly_coeffs(last_input, coeff_count, get<2>(I), temp);
                    } else {
                        set_uint(last_input, coeff_count, temp);
                    }

                    // Lazy subtraction here. ntt_negacyclic_harvey_lazy can take 0 < x <
                    // 4*qi input.
                    uint64_t neg_half_mod =
                            get<2>(I).value() - barrett_reduce_64(half, get<2>(I));

                    // Note: lambda function parameter must be passed by reference here
                    SEAL_ITERATE(temp, coeff_count, [&](auto &J) { J += neg_half_mod; });

#if SEAL_USER_MOD_BIT_COUNT_MAX <= 60
                    // Since SEAL uses at most 60-bit moduli, 8*qi < 2^63.
                    // This ntt_negacyclic_harvey_lazy results in [0, 4*qi).
                    uint64_t qi_lazy = get<2>(I).value() << 2;
                    ntt_negacyclic_harvey_lazy(temp, get<3>(I));
#else
                    // 2^60 < pi < 2^62, then 4*pi < 2^64, we perfrom one reduction from [0, 4*qi) to [0, 2*qi) after ntt.
        uint64_t qi_lazy = get<2>(I).value() << 1;
        ntt_negacyclic_harvey_lazy(temp, get<3>(I));

        // Note: lambda function parameter must be passed by reference here
        SEAL_ITERATE(temp, coeff_count, [&](auto &J) {
            J -= (qi_lazy & static_cast<uint64_t>(-static_cast<int64_t>(J >= qi_lazy)));
        });
#endif
                    // Lazy subtraction again, results in [0, 2*qi_lazy),
                    // The reduction [0, 2*qi_lazy) -> [0, qi) is done implicitly in
                    // multiply_poly_scalar_coeffmod.
                    SEAL_ITERATE(iter(get<0>(I), temp), coeff_count,
                                 [&](auto J) { get<0>(J) += qi_lazy - get<1>(J); });

                    // qk^(-1) * ((ct mod qi) - (ct mod qk)) mod qi
                    multiply_poly_scalar_coeffmod(get<0>(I), coeff_count, get<1>(I),
                                                  get<2>(I), get<0>(I)); });
    }

    static void switch_key_inplace(seal::Ciphertext &encrypted, seal::util::ConstRNSIter target_iter,
                                   const seal::KSwitchKeys &kswitch_keys, size_t kswitch_keys_index,
                                   seal::MemoryPoolHandle pool, const seal::SEALContext &context)
    {
        using namespace ::seal;
        using namespace ::seal::util;
        auto parms_id = encrypted.parms_id();
        auto &context_data = *context.get_context_data(parms_id);
        auto &parms = context_data.parms();
        auto &key_context_data = *context.key_context_data();
        auto &key_parms = key_context_data.parms();

        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context) ||
            !is_buffer_valid(encrypted))
        {
            throw std::invalid_argument("invalid encrypted buffer");
        }
        if (!target_iter)
        {
            throw std::invalid_argument("invalid iter");
        }
        if (!context.using_keyswitching())
        {
            throw std::invalid_argument("unsupport for keyswitching");
        }
        if (kswitch_keys.parms_id() != context.key_parms_id())
        {
            throw std::invalid_argument("invalid parameter id for ksk");
        }
        if (kswitch_keys_index >= kswitch_keys.data().size())
        {
            throw std::invalid_argument("ksk index out of range");
        }

        // Extract encryption parameters.
        size_t coeff_count = parms.poly_modulus_degree();
        size_t decomp_modulus_size = parms.coeff_modulus().size();
        auto &key_modulus = key_parms.coeff_modulus();
        size_t key_modulus_size = key_modulus.size();
        size_t rns_modulus_size = decomp_modulus_size + 1;
        auto key_ntt_tables = iter(key_context_data.small_ntt_tables());
        bool is_ntt_form = encrypted.is_ntt_form();

        if (!product_fits_in(coeff_count, rns_modulus_size, size_t(2)))
        {
            throw std::invalid_argument("dim mismatch");
        }

        auto &key_vector = kswitch_keys.data()[kswitch_keys_index];
        size_t key_component_count = key_vector[0].data().size();

        for (auto &each_key : key_vector)
        {
            if (!is_metadata_valid_for(each_key, context) ||
                !is_buffer_valid(each_key))
            {
                throw std::invalid_argument("invalid buffer");
            }
        }

        SEAL_ALLOCATE_GET_RNS_ITER(t_target, coeff_count, decomp_modulus_size, pool);
        set_uint(target_iter, decomp_modulus_size * coeff_count, t_target);

        if (is_ntt_form)
        {
            inverse_ntt_negacyclic_harvey(t_target, decomp_modulus_size,
                                          key_ntt_tables);
        }
        auto t_poly_prod(allocate_zero_poly_array(key_component_count, coeff_count,
                                                  rns_modulus_size, pool));

        SEAL_ITERATE(iter(size_t(0)), rns_modulus_size, [&](auto I)
                     {
            size_t key_index = (I == decomp_modulus_size ? key_modulus_size - 1 : I);
            size_t lazy_reduction_summand_bound =
                    size_t(SEAL_MULTIPLY_ACCUMULATE_USER_MOD_MAX);
            size_t lazy_reduction_counter = lazy_reduction_summand_bound;

            auto t_poly_lazy(
                    allocate_zero_poly_array(key_component_count, coeff_count, 2, pool));
            PolyIter accumulator_iter(t_poly_lazy.get(), 2, coeff_count);

            SEAL_ITERATE(iter(size_t(0)), decomp_modulus_size, [&](auto J) {
                SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);
                ConstCoeffIter t_operand;

                // RNS-NTT form exists in input
                if (is_ntt_form && (I == J)) {
                    t_operand = target_iter[J];
                }
                    // Perform RNS-NTT conversion
                else {
                    // No need to perform RNS conversion (modular reduction)
                    if (key_modulus[J] <= key_modulus[key_index]) {
                        set_uint(t_target[J], coeff_count, t_ntt);
                    }
                        // Perform RNS conversion (modular reduction)
                    else {
                        modulo_poly_coeffs(t_target[J], coeff_count, key_modulus[key_index],
                                           t_ntt);
                    }
                    // NTT conversion lazy outputs in [0, 4q)
                    ntt_negacyclic_harvey_lazy(t_ntt, key_ntt_tables[key_index]);
                    t_operand = t_ntt;
                }

                // Multiply with keys and modular accumulate products in a lazy fashion
                SEAL_ITERATE(
                        iter(key_vector[J].data(), accumulator_iter), key_component_count,
                        [&](auto K) {
                            if (!lazy_reduction_counter) {
                                SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)),
                                             coeff_count, [&](auto L) {
                                            unsigned long long qword[2]{0, 0};
                                            multiply_uint64(get<0>(L), get<1>(L), qword);

                                            // Accumulate product of t_operand and t_key_acc to
                                            // t_poly_lazy and reduce
                                            add_uint128(qword, get<2>(L).ptr(), qword);
                                            get<2>(L)[0] = barrett_reduce_128(
                                                    qword, key_modulus[key_index]);
                                            get<2>(L)[1] = 0;
                                        });
                            } else {
                                // Same as above but no reduction
                                SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)),
                                             coeff_count, [&](auto L) {
                                            unsigned long long qword[2]{0, 0};
                                            multiply_uint64(get<0>(L), get<1>(L), qword);
                                            add_uint128(qword, get<2>(L).ptr(), qword);
                                            get<2>(L)[0] = qword[0];
                                            get<2>(L)[1] = qword[1];
                                        });
                            }
                        });

                if (!--lazy_reduction_counter) {
                    lazy_reduction_counter = lazy_reduction_summand_bound;
                }
            });

            // PolyIter pointing to the destination t_poly_prod, shifted to the
            // appropriate modulus
            PolyIter t_poly_prod_iter(t_poly_prod.get() + (I * coeff_count),
                                      coeff_count, rns_modulus_size);

            // Final modular reduction
            SEAL_ITERATE(
                    iter(accumulator_iter, t_poly_prod_iter), key_component_count,
                    [&](auto K) {
                        if (lazy_reduction_counter == lazy_reduction_summand_bound) {
                            SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                                get<1>(L) = static_cast<uint64_t>(*get<0>(L));
                            });
                        } else {
                            // Same as above except need to still do reduction
                            SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                                get<1>(L) =
                                        barrett_reduce_128(get<0>(L).ptr(), key_modulus[key_index]);
                            });
                        }
                    }); });

        // Accumulated products are now stored in t_poly_prod
        // Perform modulus switching with scaling
        PolyIter t_poly_prod_iter(t_poly_prod.get(), coeff_count, rns_modulus_size);
        auto rns_tool = key_context_data.rns_tool();
        SEAL_ITERATE(
            iter(encrypted, t_poly_prod_iter), key_component_count, [&](auto I)
            {
                    divide_and_round_q_last_ntt_inplace(get<1>(I), decomp_modulus_size,
                                                        key_ntt_tables, rns_tool, pool);
                    if (!is_ntt_form) {
                        inverse_ntt_negacyclic_harvey(get<1>(I), decomp_modulus_size,
                                                      key_ntt_tables);
                    }
                    add_poly_coeffmod(get<0>(I), get<1>(I), decomp_modulus_size,
                                      key_modulus, get<0>(I)); });
    }

    static void switch_key_inplace_alternative(seal::Ciphertext &encrypted, seal::util::ConstRNSIter target_iter,
                                               const seal::KSwitchKeys &kswitch_keys, size_t kswitch_keys_index,
                                               seal::MemoryPoolHandle pool, const seal::SEALContext &context)
    {
        using namespace ::seal;
        using namespace ::seal::util;
        auto parms_id = encrypted.parms_id();
        auto &context_data = *context.get_context_data(parms_id);
        auto &parms = context_data.parms();
        auto &key_context_data = *context.key_context_data();
        auto &key_parms = key_context_data.parms();

        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context) ||
            !is_buffer_valid(encrypted))
        {
            throw std::invalid_argument("invalid encrypted buffer");
            // throw invalid_argument("encrypted is not valid for encryption
            // parameters");
        }
        if (!target_iter)
        {
            throw std::invalid_argument("invalid iter");
            // throw invalid_argument("target_iter");
        }
        if (!context.using_keyswitching())
        {
            throw std::invalid_argument("unsupport for keyswitching");
            // throw logic_error("keyswitching is not supported by the context");
        }

        // Don't validate all of kswitch_keys but just check the parms_id.
        if (kswitch_keys.parms_id() != context.key_parms_id())
        {
            throw std::invalid_argument("invalid parameter id for ksk");
            // throw invalid_argument("parameter mismatch");
        }

        if (kswitch_keys_index >= kswitch_keys.data().size())
        {
            throw std::invalid_argument("ksk index out of range");
            // throw out_of_range("kswitch_keys_index");
        }

        // Extract encryption parameters.
        size_t coeff_count = parms.poly_modulus_degree();
        size_t decomp_modulus_size = parms.coeff_modulus().size();
        auto &key_modulus = key_parms.coeff_modulus();
        size_t key_modulus_size = key_modulus.size();
        size_t rns_modulus_size = decomp_modulus_size + 1;
        auto key_ntt_tables = iter(key_context_data.small_ntt_tables());
        bool is_ntt_form = encrypted.is_ntt_form();
        auto modswitch_factors = key_context_data.rns_tool()->inv_q_last_mod_q();

        // Size check
        if (!product_fits_in(coeff_count, rns_modulus_size, size_t(2)))
        {
            throw std::invalid_argument("dim mismatch");
            // throw logic_error("invalid parameters");
        }

        // Prepare input
        auto &key_vector = kswitch_keys.data()[kswitch_keys_index];
        size_t key_component_count = key_vector[0].data().size();

        // Check only the used component in KSwitchKeys.
        for (auto &each_key : key_vector)
        {
            if (!is_metadata_valid_for(each_key, context) ||
                !is_buffer_valid(each_key))
            {
                throw std::invalid_argument("invalid buffer");
                // throw invalid_argument(
                //     "kswitch_keys is not valid for encryption parameters");
            }
        }

        // Create a copy of target_iter
        SEAL_ALLOCATE_GET_RNS_ITER(t_target, coeff_count, decomp_modulus_size, pool);
        set_uint(target_iter, decomp_modulus_size * coeff_count, t_target);

        // In CKKS t_target is in NTT form; switch back to normal form
        if (is_ntt_form)
        {
            inverse_ntt_negacyclic_harvey(t_target, decomp_modulus_size,
                                          key_ntt_tables);
        }

        // Temporary result
        auto t_poly_prod(allocate_zero_poly_array(key_component_count, coeff_count,
                                                  rns_modulus_size, pool));

        SEAL_ITERATE(iter(size_t(0)), rns_modulus_size, [&](auto I)
                     {
            size_t key_index = (I == decomp_modulus_size ? key_modulus_size - 1 : I);

            // Product of two numbers is up to 60 + 60 = 120 bits, so we can sum up to
            // 256 of them without reduction.
            size_t lazy_reduction_summand_bound =
                    size_t(SEAL_MULTIPLY_ACCUMULATE_USER_MOD_MAX);
            size_t lazy_reduction_counter = lazy_reduction_summand_bound;

            // Allocate memory for a lazy accumulator (128-bit coefficients)
            auto t_poly_lazy(
                    allocate_zero_poly_array(key_component_count, coeff_count, 2, pool));

            // Semantic misuse of PolyIter; this is really pointing to the data for a
            // single RNS factor
            PolyIter accumulator_iter(t_poly_lazy.get(), 2, coeff_count);

            // Multiply with keys and perform lazy reduction on product's coefficients
            SEAL_ITERATE(iter(size_t(0)), decomp_modulus_size, [&](auto J) {
                SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);
                ConstCoeffIter t_operand;

                // RNS-NTT form exists in input
                if (is_ntt_form && (I == J)) {
                    t_operand = target_iter[J];
                }
                    // Perform RNS-NTT conversion
                else {
                    // No need to perform RNS conversion (modular reduction)
                    if (key_modulus[J] <= key_modulus[key_index]) {
                        set_uint(t_target[J], coeff_count, t_ntt);
                    }
                        // Perform RNS conversion (modular reduction)
                    else {
                        modulo_poly_coeffs(t_target[J], coeff_count, key_modulus[key_index],
                                           t_ntt);
                    }
                    // NTT conversion lazy outputs in [0, 4q)
                    ntt_negacyclic_harvey_lazy(t_ntt, key_ntt_tables[key_index]);
                    t_operand = t_ntt;
                }

                // Multiply with keys and modular accumulate products in a lazy fashion
                SEAL_ITERATE(
                        iter(key_vector[J].data(), accumulator_iter), key_component_count,
                        [&](auto K) {
                            if (!lazy_reduction_counter) {
                                SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)),
                                             coeff_count, [&](auto L) {
                                            unsigned long long qword[2]{0, 0};
                                            multiply_uint64(get<0>(L), get<1>(L), qword);

                                            // Accumulate product of t_operand and t_key_acc to
                                            // t_poly_lazy and reduce
                                            add_uint128(qword, get<2>(L).ptr(), qword);
                                            get<2>(L)[0] = barrett_reduce_128(
                                                    qword, key_modulus[key_index]);
                                            get<2>(L)[1] = 0;
                                        });
                            } else {
                                // Same as above but no reduction
                                SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)),
                                             coeff_count, [&](auto L) {
                                            unsigned long long qword[2]{0, 0};
                                            multiply_uint64(get<0>(L), get<1>(L), qword);
                                            add_uint128(qword, get<2>(L).ptr(), qword);
                                            get<2>(L)[0] = qword[0];
                                            get<2>(L)[1] = qword[1];
                                        });
                            }
                        });

                if (!--lazy_reduction_counter) {
                    lazy_reduction_counter = lazy_reduction_summand_bound;
                }
            });

            // PolyIter pointing to the destination t_poly_prod, shifted to the
            // appropriate modulus
            PolyIter t_poly_prod_iter(t_poly_prod.get() + (I * coeff_count),
                                      coeff_count, rns_modulus_size);

            // Final modular reduction
            SEAL_ITERATE(
                    iter(accumulator_iter, t_poly_prod_iter), key_component_count,
                    [&](auto K) {
                        if (lazy_reduction_counter == lazy_reduction_summand_bound) {
                            SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                                get<1>(L) = static_cast<uint64_t>(*get<0>(L));
                            });
                        } else {
                            // Same as above except need to still do reduction
                            SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                                get<1>(L) =
                                        barrett_reduce_128(get<0>(L).ptr(), key_modulus[key_index]);
                            });
                        }
                    }); });

        // Accumulated products are now stored in t_poly_prod
        // Perform modulus switching with scaling
        PolyIter t_poly_prod_iter(t_poly_prod.get(), coeff_count, rns_modulus_size);
        SEAL_ITERATE(iter(encrypted, t_poly_prod_iter), key_component_count, [&](auto I)
                     {
            // Lazy reduction; this needs to be then reduced mod qi
            CoeffIter t_last(get<1>(I)[decomp_modulus_size]);
            inverse_ntt_negacyclic_harvey(t_last, key_ntt_tables[key_modulus_size - 1]);

            // Add (p-1)/2 to change from flooring to rounding.
            uint64_t qk = key_modulus[key_modulus_size - 1].value();
            uint64_t qk_half = qk >> 1;
            SEAL_ITERATE(t_last, coeff_count, [&](auto &J) {
                J = barrett_reduce_64(J + qk_half, key_modulus[key_modulus_size - 1]);
            });

            SEAL_ITERATE(iter(I, key_modulus, key_ntt_tables, modswitch_factors), decomp_modulus_size, [&](auto J) {
                SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);

                // (ct mod 4qk) mod qi
                uint64_t qi = get<1>(J).value();
                if (qk > qi)
                {
                    // This cannot be spared. NTT only tolerates input that is less than 4*modulus (i.e. qk <=4*qi).
                    modulo_poly_coeffs(t_last, coeff_count, get<1>(J), t_ntt);
                }
                else
                {
                    set_uint(t_last, coeff_count, t_ntt);
                }

                // Lazy substraction, results in [0, 2*qi), since fix is in [0, qi].
                uint64_t fix = qi - barrett_reduce_64(qk_half, get<1>(J));
                SEAL_ITERATE(t_ntt, coeff_count, [fix](auto &K) { K += fix; });

                uint64_t qi_lazy = qi << 1; // some multiples of qi
                if (is_ntt_form)
                {
                    // This ntt_negacyclic_harvey_lazy results in [0, 4*qi).
                    ntt_negacyclic_harvey_lazy(t_ntt, get<2>(J));
#if SEAL_USER_MOD_BIT_COUNT_MAX > 60
                    // Reduce from [0, 4qi) to [0, 2qi)
                    SEAL_ITERATE(t_ntt, coeff_count, [&](auto &K) { K -= SEAL_COND_SELECT(K >= qi_lazy, qi_lazy, 0); });
#else
                    // Since SEAL uses at most 60bit moduli, 8*qi < 2^63.
                    qi_lazy = qi << 2;
#endif
                }
                else
                {
                    inverse_ntt_negacyclic_harvey_lazy(get<0, 1>(J), get<2>(J));
                }

                // ((ct mod qi) - (ct mod qk)) mod qi
                SEAL_ITERATE(iter(get<0, 1>(J), t_ntt), coeff_count, [&](auto K) { get<0>(K) += qi_lazy - get<1>(K); });

                // qk^(-1) * ((ct mod qi) - (ct mod qk)) mod qi
                multiply_poly_scalar_coeffmod(get<0, 1>(J), coeff_count, get<3>(J), get<1>(J), get<0, 1>(J));
                add_poly_coeffmod(get<0, 1>(J), get<0, 0>(J), coeff_count, get<1>(J), get<0, 0>(J));
            }); });
    }

    void apply_galois_inplace(seal::Ciphertext &encrypted, uint32_t galois_elt,
                              const seal::GaloisKeys &galois_keys, const seal::SEALContext &context)
    {
        using namespace seal;
        using namespace seal::util;
        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context) ||
            !is_buffer_valid(encrypted))
        {
            std::cerr << "encrypted is not valid for encryption parameters"
                      << std::endl;
            throw std::invalid_argument("encrypted buffer invalid");
        }

        // Don't validate all of galois_keys but just check the parms_id.
        if (galois_keys.parms_id() != context.key_parms_id())
        {
            std::cerr << "galois_keys is not valid for encryption parameters"
                      << std::endl;
            throw std::invalid_argument("parameter id mismatch");
        }

        auto &context_data = *context.get_context_data(encrypted.parms_id());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_size = encrypted.size();
        // Use key_context_data where permutation tables exist since previous runs.
        auto galois_tool = context.key_context_data()->galois_tool();
        auto pool = MemoryManager::GetPool(mm_prof_opt::mm_force_thread_local);

        // Size check
        if (!product_fits_in(coeff_count, coeff_modulus_size))
        {
            // throw logic_error("invalid parameters");
            throw std::invalid_argument("dim mismatch");
        }

        // Check if Galois key is generated or not.
        if (!galois_keys.has_key(galois_elt))
        {
            // throw invalid_argument("Galois key not present");
            throw std::invalid_argument("key is not exist");
        }

        uint64_t m = mul_safe(static_cast<uint64_t>(coeff_count), uint64_t(2));

        // Verify parameters
        if (!(galois_elt & 1) || unsigned_geq(galois_elt, m))
        {
            // throw invalid_argument("Galois element is not valid");
            throw std::invalid_argument("galois elt must be an odd");
        }
        if (encrypted_size > 2)
        {
            // throw invalid_argument("encrypted size must be 2");
            throw std::invalid_argument("encrypted size must be 2");
        }

        SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, coeff_modulus_size, pool);

        // DO NOT CHANGE EXECUTION ORDER OF FOLLOWING SECTION
        // BEGIN: Apply Galois for each ciphertext
        // Execution order is sensitive, since apply_galois is not inplace!
        if (!encrypted.is_ntt_form())
        {
            // !!! DO NOT CHANGE EXECUTION ORDER!!!

            // First transform encrypted.data(0)
            auto encrypted_iter = iter(encrypted);
            galois_tool->apply_galois(encrypted_iter[0], coeff_modulus_size, galois_elt,
                                      coeff_modulus, temp);

            // Copy result to encrypted.data(0)
            set_poly(temp, coeff_count, coeff_modulus_size, encrypted.data(0));

            // Next transform encrypted.data(1)
            galois_tool->apply_galois(encrypted_iter[1], coeff_modulus_size, galois_elt,
                                      coeff_modulus, temp);
        }
        else
        {
            // !!! DO NOT CHANGE EXECUTION ORDER!!!

            // First transform encrypted.data(0)
            auto encrypted_iter = iter(encrypted);
            galois_tool->apply_galois_ntt(encrypted_iter[0], coeff_modulus_size,
                                          galois_elt, temp);

            // Copy result to encrypted.data(0)
            set_poly(temp, coeff_count, coeff_modulus_size, encrypted.data(0));

            // Next transform encrypted.data(1)
            galois_tool->apply_galois_ntt(encrypted_iter[1], coeff_modulus_size,
                                          galois_elt, temp);
        }

        // Wipe encrypted.data(1)
        set_zero_poly(coeff_count, coeff_modulus_size, encrypted.data(1));

        // END: Apply Galois for each ciphertext
        // REORDERING IS SAFE NOW

        // Calculate (temp * galois_key[0], temp * galois_key[1]) + (ct[0], 0)

        if (std::any_of(*temp, *temp + coeff_modulus_size * coeff_count,
                        [](uint64_t u)
                        { return u > 0; }))
        {
            switch_key_inplace(encrypted, temp,
                               static_cast<const KSwitchKeys &>(galois_keys),
                               GaloisKeys::get_index(galois_elt), pool, context);
        }
    }

    void set_zero_ct(seal::Ciphertext &ct, const seal::SEALContext &context, seal::parms_id_type parms_id, bool is_ntt_form)
    {
        ct.resize(context, parms_id, 2);
        auto context_data = context.get_context_data(parms_id);
        auto &parms = context_data->parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        size_t coeff_modulus_size = parms.coeff_modulus().size();
        size_t total_coeff_count = poly_modulus_degree * coeff_modulus_size * 2;
        std::fill_n(ct.data(), total_coeff_count, 0);
        ct.scale() = 1.;
        ct.is_ntt_form() = is_ntt_form;
    }

    void add(const seal::Ciphertext &ctxt1, const seal::Ciphertext &ctxt2, seal::Ciphertext &destination,
             const seal::SEALContext &context)
    {
        if (!seal::util::are_close(ctxt1.scale(), ctxt2.scale()))
            throw std::invalid_argument("add: ciphertexts scale mismatch");
        if (ctxt1.parms_id() != ctxt2.parms_id())
            throw std::invalid_argument("ciphertexts parms id must be same");
        if (ctxt1.is_ntt_form() != ctxt2.is_ntt_form())
            throw std::invalid_argument("ct1 and ct2 must be simultaneous in NTT or INTT form");
        auto context_data_ptr = context.get_context_data(ctxt1.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();

        // restrict the ctxt size to 2 (not necessary)
        if (ctxt1.size() != 2 || ctxt2.size() != 2)
        {
            throw std::invalid_argument("ciphertext size must be 2 here");
        }

        destination.resize(context, context_data.parms_id(), 2);
        seal::util::add_poly_coeffmod(ctxt1, ctxt2, 2, coeff_modulus, destination);
        destination.scale() = ctxt1.scale();
        destination.is_ntt_form() = ctxt1.is_ntt_form();
    }

    void add_inplace(seal::Ciphertext &destination, const seal::Ciphertext &ctxt, const seal::SEALContext &context)
    {
        using namespace seal::util;
        if (!are_close(destination.scale(), ctxt.scale()))
            throw std::invalid_argument("add_inplace: ciphertexts scale mismatch");
        if (destination.parms_id() != ctxt.parms_id())
            throw std::invalid_argument("ciphertexts parms_id must be same");
        if (destination.is_ntt_form() != ctxt.is_ntt_form())
            throw std::invalid_argument("NTT form mismatch");

        auto context_data = context.get_context_data(destination.parms_id());
        auto &parms = context_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();

        size_t dest_size = destination.size();
        size_t ctxt_size = ctxt.size();

        // not necessary
        if (dest_size != 2 || ctxt_size != 2)
            throw std::invalid_argument("ciphertext size must be 2 here");
        // not necessary here
        destination.resize(context, context_data->parms_id(), dest_size);
        add_poly_coeffmod(destination, ctxt, 2, coeff_modulus, destination);
    }

    void add_plain(const seal::Ciphertext &ctxt, const seal::Plaintext &plain, seal::Ciphertext &destination,
                   const seal::SEALContext &context)
    {
        using namespace seal::util;
        if (!seal::util::are_close(ctxt.scale(), plain.scale()))
            throw std::invalid_argument("add_plain: ctxt, plain scale mismatch");
        destination = ctxt;
        size_t N = ctxt.poly_modulus_degree();
        size_t L = ctxt.coeff_modulus_size();
        if (plain.coeff_count() != N * L)
            throw std::invalid_argument("add_plain: plain coeff_count mismatch");
        auto context_data = context.get_context_data(plain.parms_id());
        auto &coeff_modulus = context_data->parms().coeff_modulus();
        ConstRNSIter plain_iter(plain.data(), N);
        RNSIter dest_iter = *iter(destination);
        add_poly_coeffmod(dest_iter, plain_iter, L, coeff_modulus, dest_iter);
    }

    void add_plain_inplace(seal::Ciphertext &destination, const seal::Plaintext &plain, const seal::SEALContext &context)
    {
        using namespace seal::util;
        size_t N = destination.poly_modulus_degree();
        size_t L = destination.coeff_modulus_size();
        if (plain.coeff_count() != N * L)
            throw std::invalid_argument("add_plain_inplace: plain coeff_count mismatch");
        auto context_data = context.get_context_data(plain.parms_id());
        auto &coeff_modulus = context_data->parms().coeff_modulus();
        ConstRNSIter plain_iter(plain.data(), N);
        RNSIter dest_iter = *iter(destination);
        add_poly_coeffmod(dest_iter, plain_iter, L, coeff_modulus, dest_iter);
    }

    void sub(const seal::Ciphertext &ctxt1, const seal::Ciphertext &ctxt2, seal::Ciphertext &destination,
             const seal::SEALContext &context)
    {
        if (!seal::util::are_close(ctxt1.scale(), ctxt2.scale()))
            throw std::invalid_argument("ciphertexts scale mismatch");
        if (ctxt1.parms_id() != ctxt2.parms_id())
            throw std::invalid_argument("ciphertexts parms id must be the same");
        auto context_data_ptr = context.get_context_data(ctxt1.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();

        // restrict the ctxt size to 2 (not necessary)
        if (ctxt1.size() != 2 || ctxt2.size() != 2)
            throw std::invalid_argument("ciphertext size must be 2 here");

        destination.resize(context, context_data.parms_id(), 2);
        seal::util::sub_poly_coeffmod(ctxt1, ctxt2, 2, coeff_modulus, destination);
        destination.scale() = ctxt1.scale();
    }

    void sub_plain(const seal::Ciphertext &ctxt, const seal::Plaintext &plain, seal::Ciphertext &destination,
                   const seal::SEALContext &context)
    {
        using namespace seal::util;
        if (!seal::util::are_close(ctxt.scale(), plain.scale()))
            throw std::invalid_argument("sub_plain: ctxt, plain scale mismatch");
        destination = ctxt;
        size_t N = ctxt.poly_modulus_degree();
        size_t L = ctxt.coeff_modulus_size();
        if (plain.coeff_count() != N * L)
            throw std::invalid_argument("sub_plain: plain coeff_count mismatch");
        auto context_data = context.get_context_data(plain.parms_id());
        auto &coeff_modulus = context_data->parms().coeff_modulus();
        ConstRNSIter plain_iter(plain.data(), N);
        RNSIter dest_iter = *iter(destination);
        sub_poly_coeffmod(dest_iter, plain_iter, L, coeff_modulus, dest_iter);
    }

    void sub_plain_inplace(seal::Ciphertext &destination, const seal::Plaintext &plain, const seal::SEALContext &context)
    {
        using namespace seal::util;
        if (!seal::util::are_close(destination.scale(), plain.scale()))
            throw std::invalid_argument("sub_plain_inplace: destination, plain scale mismatch");
        size_t N = destination.poly_modulus_degree();
        size_t L = destination.coeff_modulus_size();
        if (plain.coeff_count() != N * L)
            throw std::invalid_argument("sub_plain_inplace: plain coeff_count mismatch");
        auto context_data = context.get_context_data(plain.parms_id());
        auto &coeff_modulus = context_data->parms().coeff_modulus();
        ConstRNSIter plain_iter(plain.data(), N);
        RNSIter dest_iter = *iter(destination);
        sub_poly_coeffmod(dest_iter, plain_iter, L, coeff_modulus, dest_iter);
    }

    void multiply_power_of_x(const seal::Ciphertext &encrypted, seal::Ciphertext &destination,
                             uint64_t exponent_of_x, const seal::SEALContext &context)
    {
        if (encrypted.is_ntt_form())
            throw std::invalid_argument("Encrypted cannot be in NTT form");
        auto parms_id = encrypted.parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();

        uint64_t coeff_mod_count = parms.coeff_modulus().size();
        uint64_t coeff_count = parms.poly_modulus_degree();

        uint64_t encrypted_count = encrypted.size();
        destination = encrypted;

        for (size_t i = 0; i < encrypted_count; ++i)
        {
            for (size_t j = 0; j < coeff_mod_count; ++j)
            {
                seal::util::negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count),
                                                           coeff_count, exponent_of_x, parms.coeff_modulus()[j],
                                                           destination.data(i) + (j * coeff_count));
            }
        }
    }

    void pack_rlwes(const std::vector<seal::Ciphertext> &rlwe_ctxts, uint64_t rlwe_num, uint64_t init_exp_of_x,
                    const seal::GaloisKeys &galois_keys, seal::Ciphertext &destination,
                    const seal::SEALContext &context, uint64_t multi_thread_count)
    {
        auto context_data_ptr = context.get_context_data(rlwe_ctxts[0].parms_id());
        int log_rlwe_ctxt_num = seal::util::get_power_of_two(rlwe_num);
        if (log_rlwe_ctxt_num == -1)
            throw std::invalid_argument("The number of ciphertexts must be a power of 2");
        if (rlwe_num == 1)
        {
            destination = rlwe_ctxts[0];
            return;
        }

        std::vector<seal::Ciphertext> temp(rlwe_num >> 1);
        uint32_t iter_num = log_rlwe_ctxt_num;

        auto pack_two_rlwes = [&](const seal::Ciphertext &rlwe1, const seal::Ciphertext &rlwe2,
                                  seal::Ciphertext &destination, uint32_t exp_of_x)
        {
            seal::Ciphertext temp, temp_add, temp_sub;
            uint32_t N = rlwe1.poly_modulus_degree();
            multiply_power_of_x(rlwe2, temp, exp_of_x, context);
            uint32_t galois_elt = N / exp_of_x + 1;
            add(rlwe1, temp, temp_add, context);
            sub(rlwe1, temp, temp_sub, context);
            apply_galois_inplace(temp_sub, galois_elt, galois_keys, context);
            add(temp_add, temp_sub, destination, context);
        };

        size_t thread_count{0};
        for (uint32_t i = 0; i < iter_num; ++i)
        {
            size_t steps_j = rlwe_num >> (i + 1);
            for (size_t j = 0; j < steps_j; j += thread_count)
            {
                size_t step_last = steps_j - j;
                thread_count = (step_last < multi_thread_count) ? step_last : multi_thread_count;
                std::vector<std::thread> threads(thread_count);
                for (size_t k = 0; k < thread_count; ++k)
                {
                    if (i == 0)
                        threads[k] = std::thread(pack_two_rlwes, std::ref(rlwe_ctxts[(j + k) << 1]),
                                                 std::ref(rlwe_ctxts[((j + k) << 1) + 1]), std::ref(temp[j + k]), init_exp_of_x);
                    else
                        threads[k] = std::thread(pack_two_rlwes, std::ref(temp[(j + k) << 1]), std::ref(temp[((j + k) << 1) + 1]),
                                                 std::ref(temp[j + k]), init_exp_of_x >> i);
                }
                std::for_each(threads.begin(), threads.end(), [](std::thread &t)
                              { t.join(); });
            }
        }
        destination = temp[0];
    }

    void pack_rlwes(const std::vector<seal::Ciphertext> &rlwe_ctxts, uint32_t h, const seal::GaloisKeys &galois_keys,
                    seal::Ciphertext &destination, const seal::SEALContext &context, uint32_t threads_count)
    {
        uint32_t N = rlwe_ctxts[0].poly_modulus_degree();
        if (N < (h << 1))
            throw std::invalid_argument("invalid input h");
        uint64_t init_exp_of_x = N >> (h + 1);
        pack_rlwes(rlwe_ctxts, rlwe_ctxts.size(), init_exp_of_x, galois_keys, destination, context, threads_count);
    }

    void PackRLWEs(const std::vector<seal::Ciphertext> &rlwe_ctxts, uint32_t h, const seal::GaloisKeys &galois_keys,
                   seal::Ciphertext &destination, const seal::SEALContext &context, uint32_t threads)
    {
        uint32_t ct_num = rlwe_ctxts.size();
        PackRLWEs(rlwe_ctxts.data(), ct_num, h, galois_keys, destination, context, threads);
    }

    void PackRLWEs(const seal::Ciphertext *rlwe_ctxts, uint32_t ct_num, uint32_t h, const seal::GaloisKeys &galois_keys,
                   seal::Ciphertext &destination, const seal::SEALContext &context, uint32_t threads_count)
    {
        auto context_data_ptr = context.get_context_data(rlwe_ctxts[0].parms_id());
        uint32_t N = rlwe_ctxts[0].poly_modulus_degree();
        if (ct_num * (1 << h) > N)
            throw std::invalid_argument("do not have enough coefficients");
        if (ct_num == 1)
        {
            destination = rlwe_ctxts[0];
            return;
        }

        uint32_t log_pad_rlwe_num = seal::util::get_significant_bit_count(ct_num - 1);
        uint32_t pad_rlwe_num = 1UL << log_pad_rlwe_num;
        std::vector<seal::Ciphertext> temp(pad_rlwe_num >> 1);

        // If rlwe2 is zero, it will be omitted.
        auto pack_two_rlwes = [&](const seal::Ciphertext &rlwe1, const seal::Ciphertext &rlwe2,
                                  bool rlwe2_is_zero, seal::Ciphertext &destination, uint32_t exp_of_x){
            uint32_t N = rlwe1.poly_modulus_degree();
            uint32_t galois_elt = N / exp_of_x + 1;
            seal::Ciphertext temp;
            if (rlwe2_is_zero){
                temp = rlwe1;
                apply_galois_inplace(temp, galois_elt, galois_keys, context);
                add(rlwe1, temp, destination, context);
            }else{
                seal::Ciphertext temp_add, temp_sub;
                multiply_power_of_x(rlwe2, temp, exp_of_x, context);
                add(rlwe1, temp, temp_add, context);
                sub(rlwe1, temp, temp_sub, context);
                apply_galois_inplace(temp_sub, galois_elt, galois_keys, context);
                add(temp_add, temp_sub, destination, context);
            }
        };

        uint32_t thread_count = 0;
        uint32_t init_exp_of_x = N >> (h + 1);
        for (uint32_t i = 0; i < log_pad_rlwe_num; ++i)
        {
            uint32_t steps_j = pad_rlwe_num >> (i + 1);
            for (uint32_t j = 0; j < steps_j; j += thread_count)
            {
                uint32_t step_last = steps_j - j;
                thread_count = (step_last < threads_count) ? step_last : threads_count;
                std::vector<std::thread> threads_pool(thread_count);
                for (size_t k = 0; k < thread_count; ++k){
                    if (i == 0){
                        uint32_t rlwe1_index = j + k;
//                        std::cout << "rlwe1_index: " << rlwe1_index << std::endl;
                        uint32_t bit_rev_r1 = seal::util::reverse_bits(rlwe1_index, log_pad_rlwe_num);
                        uint32_t rlwe2_index = seal::util::reverse_bits(bit_rev_r1 + 1, log_pad_rlwe_num);
//                        std::cout << "rlwe2_index: " << rlwe2_index << std::endl;
                        uint32_t dest_index = bit_rev_r1 >> 1;
                        bool is_zero = (rlwe2_index >= ct_num);
                        if (is_zero){
                            threads_pool[k] = std::thread(pack_two_rlwes, std::cref(rlwe_ctxts[rlwe1_index]),
                                                          std::cref(rlwe_ctxts[rlwe1_index]) /* no use */, is_zero, std::ref(temp[dest_index]), init_exp_of_x);
                        }else {
                            threads_pool[k] = std::thread(pack_two_rlwes, std::cref(rlwe_ctxts[rlwe1_index]),
                                                          std::cref(rlwe_ctxts[rlwe2_index]),
                                                          is_zero, std::ref(temp[dest_index]), init_exp_of_x);
                        }
                    }else{
                        threads_pool[k] = std::thread(pack_two_rlwes, std::cref(temp[(j + k) << 1]), std::cref(temp[((j + k) << 1) + 1]),
                                                      false, std::ref(temp[j + k]), init_exp_of_x >> i);
                    }
                }
                std::for_each(threads_pool.begin(), threads_pool.end(), [](std::thread &t){t.join();});
            }
        }
        destination = temp[0];
    }

    void expand(const seal::Ciphertext &ctxt, std::vector<seal::Ciphertext> &destination, int ell,
                const seal::GaloisKeys &galois_keys, const seal::SEALContext &context, int threads)
    {
        uint32_t N = ctxt.poly_modulus_degree();
        int logN = seal::util::get_power_of_two(N);
        if (ell < 0 || ell > logN)
            throw std::invalid_argument("invalid ell, ell should be in range [0, logN]");

        uint32_t ct_num = 1UL << ell;
        destination.resize(ct_num);
        destination[0] = ctxt;

        auto expand_to_two_rlwes = [&](const seal::Ciphertext &input_ct, seal::Ciphertext &dest1,
                seal::Ciphertext &dest2, uint32_t j){
            seal::Ciphertext temp0, temp1;
            temp0 = input_ct;
            apply_galois_inplace(temp0, (N >> j) + 1, galois_keys, context);
            sub(input_ct, temp0, temp1, context);
            add(temp0, input_ct, dest1, context);
            multiply_power_of_x(temp1, dest2, 2 * N - (1 << j), context);
        };

        for (uint32_t j = 0; j < ell; ++j){
            for (uint32_t k = 0; k < (1UL << j); k += threads){
                uint32_t step_last = (1UL << j) - k;
                uint32_t thread_count = (step_last < threads) ? step_last : threads;
                std::vector<std::thread> threads_pool(thread_count);
                for (uint32_t h = 0; h < thread_count; ++h){
                    threads_pool[h] = std::thread(expand_to_two_rlwes, std::cref(destination[k + h]),
                                                  std::ref(destination[k + h]), std::ref(destination[(k + h) + (1 << j)]), j);
                }
                std::for_each(threads_pool.begin(), threads_pool.end(), [](std::thread &t){t.join();});
            }
        }
    }

    void Expand(const seal::Ciphertext &ctxt, std::vector<seal::Ciphertext> &destination, uint32_t ell, uint64_t expand_factor,
                const seal::GaloisKeys &galois_keys, const seal::SEALContext &context, uint32_t threads)
    {
        uint32_t N = ctxt.poly_modulus_degree();
        uint32_t logN = seal::util::get_power_of_two(N);

        if (ell > logN)
            throw std::invalid_argument("invalid parameter ell, it should be in [0, logN]");
        if (ctxt.is_ntt_form()){
            seal::Ciphertext dup_ctxt = ctxt;
            mul_inv_pow2_inplace(dup_ctxt, context, expand_factor);
            transform_from_ntt_inplace(dup_ctxt, context);
            expand(dup_ctxt, destination, ell, galois_keys, context, threads);
            /// TODO: add multi-thread
            for (auto &ct : destination)
                transform_to_ntt_inplace(ct, context);
        }
        else{
            seal::Ciphertext dup_ct = ctxt;
            mul_inv_pow2_inplace(dup_ct, context, expand_factor);
            expand(dup_ct, destination, ell, galois_keys, context, threads);
        }
    }

    void rescale_to_next_inplace(seal::Ciphertext &encrypted, const seal::SEALContext &context)
    {
        using namespace seal;
        using namespace seal::util;
        auto context_data_ptr = context.get_context_data(encrypted.parms_id());
        auto pool = MemoryManager::GetPool(mm_prof_opt::mm_force_thread_local);

        auto &context_data = *context_data_ptr;
        auto &next_context_data = *context_data.next_context_data();
        auto &next_parms = next_context_data.parms();
        auto rns_tool = context_data.rns_tool();

        size_t encrypted_size = encrypted.size();
        size_t coeff_count = next_parms.poly_modulus_degree();
        size_t next_coeff_modulus_size = next_parms.coeff_modulus().size();

        Ciphertext encrypted_copy(pool);
        encrypted_copy = encrypted;
        bool is_ntt_form = encrypted.is_ntt_form();

        if (is_ntt_form)
        {
            SEAL_ITERATE(iter(encrypted_copy), encrypted_size, [&](auto I)
                         { rns_tool->divide_and_round_q_last_ntt_inplace(I, context_data.small_ntt_tables(), pool); });
        }
        else
        {
            SEAL_ITERATE(iter(encrypted_copy), encrypted_size, [&](auto I)
                         { rns_tool->divide_and_round_q_last_inplace(I, pool); });
        }

        encrypted.resize(context, next_context_data.parms_id(), encrypted_size);
        SEAL_ITERATE(iter(encrypted_copy, encrypted), encrypted_size, [&](auto I)
                     { set_poly(get<0>(I), coeff_count, next_coeff_modulus_size, get<1>(I)); });
    }

    void multiply_plain_ntt(const seal::Ciphertext &encrypted_ntt, const seal::Plaintext &plain_ntt,
                            seal::Ciphertext &destination, const seal::SEALContext &context)
    {
        using namespace seal;
        using namespace seal::util;
        if (!plain_ntt.is_ntt_form())
            throw std::invalid_argument("plain_ntt is not in NTT form");
        if (encrypted_ntt.parms_id() != plain_ntt.parms_id())
            throw std::invalid_argument("encrypted_ntt and plain_ntt parameter mismatch");

        auto context_data_ptr = context.get_context_data(encrypted_ntt.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_ntt_size = encrypted_ntt.size();

        ConstRNSIter plain_ntt_iter(plain_ntt.data(), coeff_count);
        destination = encrypted_ntt;
        SEAL_ITERATE(iter(destination), encrypted_ntt_size, [&](auto I)
                     { dyadic_product_coeffmod(I, plain_ntt_iter, coeff_modulus_size, coeff_modulus, I); });
        destination.scale() = encrypted_ntt.scale() * plain_ntt.scale();
    }

    // customized functionality, be careful when use this function
    void mul_inv_pow2_inplace(seal::Ciphertext &encrypted,
                              const seal::SEALContext &context, uint64_t pack_coeff)
    {
        using namespace seal::util;
        int bits = get_power_of_two(pack_coeff);
        if (bits == -1)
            throw std::invalid_argument("pack coeff must be a power of two");
        if (bits == 0)
            return;

        auto context_data_ptr = context.get_context_data(encrypted.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = coeff_modulus.size();
        size_t N = parms.poly_modulus_degree();

        std::vector<uint64_t> inv_2_bits_mod_q(coeff_modulus_size);
        for (size_t i = 0; i < coeff_modulus_size; ++i)
            inv_2_bits_mod_q[i] = coeff_modulus[i].value() - (coeff_modulus[i].value() >> bits);

        SEAL_ITERATE(iter(encrypted), encrypted.size(), [&](auto I)
                     { SEAL_ITERATE(iter(I, size_t(0), coeff_modulus), coeff_modulus.size(), [&](auto J)
                                    { multiply_poly_scalar_coeffmod(get<0>(J), N, inv_2_bits_mod_q[get<1>(J)], get<2>(J), get<0>(J)); }); });
    }

    void multiply_scalar(const seal::Ciphertext &encrypted, int64_t scalar, seal::Ciphertext &destination,
                         const seal::SEALContext &context)
    {
        using namespace seal;
        using namespace seal::util;

        auto context_data_ptr = context.get_context_data(encrypted.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = coeff_modulus.size();

        std::vector<uint64_t> scalar_in_Zq(coeff_modulus_size);
        for (size_t i = 0; i < coeff_modulus_size; ++i)
        {
            if (scalar >= 0)
                scalar_in_Zq[i] = scalar;
            else
                scalar_in_Zq[i] = int64_t(coeff_modulus[i].value()) + scalar;
        }
        destination = encrypted;
        size_t N = parms.poly_modulus_degree();
        SEAL_ITERATE(iter(encrypted, destination), encrypted.size(), [&](auto I)
                     { SEAL_ITERATE(iter(I, size_t(0), coeff_modulus), coeff_modulus_size, [&](auto J)
                                    { multiply_poly_scalar_coeffmod(get<0, 0>(J), N, scalar_in_Zq[get<1>(J)], get<2>(J), get<0, 1>(J)); }); });
    }

    // customized functionality, be careful when use this function
    void mul_uint_inplace(seal::Ciphertext &encrypted, uint64_t factor, const seal::SEALContext &context)
    {
        using namespace seal::util;
        auto context_data = context.get_context_data(encrypted.parms_id());
        auto &parms = context_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t N = parms.poly_modulus_degree();

        SEAL_ITERATE(iter(encrypted), encrypted.size(), [&](auto I)
                     { SEAL_ITERATE(iter(I, coeff_modulus), coeff_modulus.size(), [&](auto J)
                                    { multiply_poly_scalar_coeffmod(get<0>(J), N, factor, get<1>(J), get<0>(J)); }); });
    }

    void transform_from_ntt_inplace(seal::Ciphertext &ct, const seal::SEALContext &context)
    {
        seal::parms_id_type parms_id = ct.parms_id();
        auto context_data = context.get_context_data(parms_id);
        auto ntt_tables = context_data->small_ntt_tables();
        seal::util::inverse_ntt_negacyclic_harvey(ct, ct.size(), ntt_tables);
        ct.is_ntt_form() = false;
    }

    void transform_to_ntt_inplace(seal::Ciphertext &ct, const seal::SEALContext &context)
    {
        seal::parms_id_type parms_id = ct.parms_id();
        auto context_data = context.get_context_data(parms_id);
        auto ntt_tables = context_data->small_ntt_tables();
        seal::util::ntt_negacyclic_harvey(ct, ct.size(), ntt_tables);
        ct.is_ntt_form() = true;
    }

    void drop_unrelated_coeffs(seal::Ciphertext &ctxt, uint64_t vec_size, const seal::SEALContext &context)
    {
        using namespace seal::util;
        auto context_data = context.get_context_data(ctxt.parms_id());
        auto poly_modulus_degree = context_data->parms().poly_modulus_degree();
        auto coeff_modulus_size = context_data->parms().coeff_modulus().size();
        if (vec_size == 0 || vec_size > poly_modulus_degree)
            throw std::invalid_argument("vector size out of range");
        bool is_ntt_form = ctxt.is_ntt_form();
        if (is_ntt_form)
            transform_from_ntt_inplace(ctxt, context);
        uint64_t vec_size_extend = 1UL << get_significant_bit_count(vec_size - 1);
        uint64_t stride_interval = poly_modulus_degree / vec_size_extend;
        for (size_t mod = 0; mod < coeff_modulus_size; ++mod)
        {
            for (size_t i = 0; i < vec_size_extend; i++)
                std::fill_n(ctxt.data(0) + mod * poly_modulus_degree + i * stride_interval + 1, stride_interval - 1, 0);
        }
        if (is_ntt_form)
            transform_to_ntt_inplace(ctxt, context);
    }

    void encode_to_coeff_double(seal::Plaintext &plain, const double *real_vec, size_t real_size, double scale,
                         const seal::SEALContext &context, seal::parms_id_type parms_id, seal::MemoryPoolHandle pool)
    {
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &parms = context_data_ptr->parms();
        auto N = parms.poly_modulus_degree();
        if (!context_data_ptr)
        {
            throw std::invalid_argument("parms_id is not valid for encryption parameters");
        }
        if (!real_vec || real_size == 0)
        {
            throw std::invalid_argument("vector cannot be null");
        }
        if (real_size > N)
        {
            throw std::invalid_argument("vector size is too long (> N) ");
        }
        if (!pool)
        {
            throw std::invalid_argument("pool is uninitialized");
        }
        auto &context_data = *context_data_ptr;
        auto &coeff_modulus = parms.coeff_modulus();
        std::size_t coeff_modulus_size = coeff_modulus.size();
        std::size_t coeff_count = parms.poly_modulus_degree();

        // Quick sanity check
        if (!seal::util::product_fits_in(coeff_modulus_size, coeff_count))
        {
            throw std::logic_error("invalid parameters");
        }

        // Check that scale is positive and not too large
        if (scale <= 0 || (static_cast<int>(log2(scale)) + 1 >= context_data.total_coeff_modulus_bit_count()))
        {
            throw std::invalid_argument("scale out of bounds");
        }

        auto ntt_tables = context_data.small_ntt_tables();
        auto real_values = seal::util::allocate<double>(N, pool, 0);
        for (std::size_t i = 0; i < real_size; ++i)
        {
            real_values[i] = real_vec[i] * scale;
        }

        double max_coeff = 0;
        for (std::size_t i = 0; i < real_size; ++i)
        {
            max_coeff = std::max<>(max_coeff, std::fabs(real_vec[i]));
        }

        int max_coeff_bit_count = static_cast<int>(std::ceil(std::log2(std::max<>(max_coeff, 1.0)))) + 1;
        if (max_coeff_bit_count >= context_data.total_coeff_modulus_bit_count())
        {
            throw std::invalid_argument("encoded values are too large");
        }

        double two_pow_64 = std::pow(2.0, 64);

        plain.parms_id() = seal::parms_id_zero;
        plain.resize(seal::util::mul_safe(coeff_count, coeff_modulus_size));

        if (max_coeff_bit_count <= 64)
        {
            for (std::size_t i = 0; i < N; i++)
            {
                double coeffd = std::round(real_values[i]);
                bool is_negative = std::signbit(coeffd);

                std::uint64_t coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));

                if (is_negative)
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = seal::util::negate_uint_mod(
                            seal::util::barrett_reduce_64(coeffu, coeff_modulus[j]), coeff_modulus[j]);
                    }
                }
                else
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = seal::util::barrett_reduce_64(coeffu, coeff_modulus[j]);
                    }
                }
            }
        }
        else if (max_coeff_bit_count <= 128)
        {
            for (std::size_t i = 0; i < N; i++)
            {
                double coeffd = std::round(real_values[i]);
                bool is_negative = std::signbit(coeffd);
                coeffd = std::fabs(coeffd);

                std::uint64_t coeffu[2]{static_cast<std::uint64_t>(std::fmod(coeffd, two_pow_64)),
                                        static_cast<std::uint64_t>(coeffd / two_pow_64)};

                if (is_negative)
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = seal::util::negate_uint_mod(
                            seal::util::barrett_reduce_128(coeffu, coeff_modulus[j]), coeff_modulus[j]);
                    }
                }
                else
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = seal::util::barrett_reduce_128(coeffu, coeff_modulus[j]);
                    }
                }
            }
        }
        else
        {
            auto coeffu(seal::util::allocate_uint(coeff_modulus_size, pool));
            for (std::size_t i = 0; i < N; i++)
            {
                double coeffd = std::round(real_values[i]);
                bool is_negative = std::signbit(coeffd);
                coeffd = std::fabs(coeffd);

                // We are at this point guaranteed to fit in the allocated space
                seal::util::set_zero_uint(coeff_modulus_size, coeffu.get());
                auto coeffu_ptr = coeffu.get();
                while (coeffd >= 1)
                {
                    *coeffu_ptr++ = static_cast<std::uint64_t>(std::fmod(coeffd, two_pow_64));
                    coeffd /= two_pow_64;
                }

                // Next decompose this coefficient
                context_data.rns_tool()->base_q()->decompose(coeffu.get(), pool);

                // Finally replace the sign if necessary
                if (is_negative)
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = seal::util::negate_uint_mod(coeffu[j], coeff_modulus[j]);
                    }
                }
                else
                {
                    for (std::size_t j = 0; j < coeff_modulus_size; j++)
                    {
                        plain[i + (j * coeff_count)] = coeffu[j];
                    }
                }
            }
        }

        // transform to ntt
        for (std::size_t i = 0; i < coeff_modulus_size; ++i)
        {
            seal::util::ntt_negacyclic_harvey(plain.data(i * coeff_count), ntt_tables[i]);
        }
        plain.parms_id() = context_data.parms_id();
        plain.scale() = scale;
    }

    void encode_to_coeff128(seal::Plaintext &plain, const uint64_t *vec, size_t vec_len, const AuxParms &aux_parms,
                            seal::parms_id_type parms_id, const seal::SEALContext &context)
    {
        auto context_data = context.get_context_data(parms_id);
        auto &parms = context_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = coeff_modulus.size();
        auto ntt_table = context_data->small_ntt_tables();

        uint32_t poly_mod_degree = parms.poly_modulus_degree();

        plain.parms_id() = seal::parms_id_zero;
        size_t buffer_size = poly_mod_degree * coeff_modulus_size;
        plain.resize(buffer_size);
        std::fill_n(plain.data(), buffer_size, 0);

        for (size_t i = 0; i < coeff_modulus_size; ++i){
            auto offset = i * poly_mod_degree;
            uint64_t pow_2_128_mod_qi = seal::util::multiply_uint_mod(aux_parms.pow_2_64_mod_qi_[i],
                                                                      aux_parms.pow_2_64_mod_qi_[i], coeff_modulus[i]);
            for (size_t j = 0; j < poly_mod_degree; ++j){

                // Note:
                // [0, 2^128) = [0, 2^127) \cup [2^127, 2^128)
                // Regard [0, 2^127) as positive numbers, [2^127, 2^128) as negative numbers.
                // For a \in [2^127, 2^128), it represents the negative number a - 2^128.

                // positive number
                if (vec[(j << 1) + 1] < (1ULL << 63)){
                    plain[offset + j] = seal::util::barrett_reduce_128(vec + (j << 1), coeff_modulus[i]);
                }
                else{
                    // A negative number a will be represented by 2^128 + a
                    // 2^128 + a mod qi
                    uint64_t v = seal::util::barrett_reduce_128(vec + (j << 1), coeff_modulus[i]);
                    // subtract the 2^128
                    plain[offset + j] = seal::util::sub_uint_mod(v, pow_2_128_mod_qi, coeff_modulus[i]);
                }
            }

            // transform to ntt form, 'lazy' version could be ok in some scenarios.
            seal::util::ntt_negacyclic_harvey(plain.data(offset), ntt_table[i]);
        }
        plain.parms_id() = parms_id;
        plain.scale() = 1.;
    }

    void decode_from_coeff_double(const seal::Plaintext &plain, double *destination, size_t vec_size,
                           const seal::SEALContext &context, seal::MemoryPoolHandle pool)
    {
        // verify parameters
        if (!seal::is_valid_for(plain, context))
        {
            throw std::invalid_argument("plain is not valid for encryption parameters");
        }
        if (!plain.is_ntt_form())
        {
            throw std::invalid_argument("plain is not in NTT form");
        }
        if (!destination || vec_size == 0)
        {
            throw std::invalid_argument("destination cannot be null or vec_size == 0");
        }
        if (!pool)
        {
            throw std::invalid_argument("pool is uninitialized");
        }

        auto &context_data = *context.get_context_data(plain.parms_id());
        auto &parms = context_data.parms();
        std::size_t coeff_modulus_size = parms.coeff_modulus().size();
        std::size_t coeff_count = parms.poly_modulus_degree();
        std::size_t rns_poly_uint64_count = seal::util::mul_safe(coeff_count, coeff_modulus_size);
        auto ntt_tables = context_data.small_ntt_tables();

        // Check that scale is positive and not too large
        if (plain.scale() <= 0 ||
            (static_cast<int>(log2(plain.scale())) >= context_data.total_coeff_modulus_bit_count()))
        {
            throw std::invalid_argument("scale out of bounds");
        }

        auto decryption_modulus = context_data.total_coeff_modulus();
        auto upper_half_threshold = context_data.upper_half_threshold();
        int logn = seal::util::get_power_of_two(coeff_count);

        // Quick sanity check
        if ((logn < 0) || (coeff_count < SEAL_POLY_MOD_DEGREE_MIN) || (coeff_count > SEAL_POLY_MOD_DEGREE_MAX))
        {
            throw std::logic_error("invalid parameters");
        }

        double inv_scale = double(1.0) / plain.scale();
        // Create mutable copy of input
        auto plain_copy(seal::util::allocate_uint(rns_poly_uint64_count, pool));
        seal::util::set_uint(plain.data(), rns_poly_uint64_count, plain_copy.get());

        // Transform each polynomial from NTT domain
        for (std::size_t i = 0; i < coeff_modulus_size; i++)
        {
            seal::util::inverse_ntt_negacyclic_harvey(plain_copy.get() + (i * coeff_count), ntt_tables[i]);
        }

        // CRT-compose the polynomial
        context_data.rns_tool()->base_q()->compose_array(plain_copy.get(), coeff_count, pool);

        // Create floating-point representations of the multi-precision integer coefficients
        double two_pow_64 = std::pow(2.0, 64);
        auto res(seal::util::allocate<double>(coeff_count, pool));
        for (std::size_t i = 0; i < coeff_count; i++)
        {
            res[i] = 0.0;
            if (seal::util::is_greater_than_or_equal_uint(
                    plain_copy.get() + (i * coeff_modulus_size), upper_half_threshold, coeff_modulus_size))
            {
                double scaled_two_pow_64 = inv_scale;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
                {
                    if (plain_copy[i * coeff_modulus_size + j] > decryption_modulus[j])
                    {
                        auto diff = plain_copy[i * coeff_modulus_size + j] - decryption_modulus[j];
                        res[i] += diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    }
                    else
                    {
                        auto diff = decryption_modulus[j] - plain_copy[i * coeff_modulus_size + j];
                        res[i] -= diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    }
                }
            }
            else
            {
                double scaled_two_pow_64 = inv_scale;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
                {
                    auto curr_coeff = plain_copy[i * coeff_modulus_size + j];
                    res[i] += curr_coeff ? static_cast<double>(curr_coeff) * scaled_two_pow_64 : 0.0;
                }
            }

            // Scaling instead incorporated above; this can help in cases
            // where otherwise pow(two_pow_64, j) would overflow due to very
            // large coeff_modulus_size and very large scale
            // res[i] = res_accum * inv_scale;
        }
        std::copy_n(res.get(), vec_size, destination);
    }

    void decode_from_coeff_double(const seal::Plaintext &plain, double *destination, size_t vec_size, uint32_t stride,
                           const seal::SEALContext &context, seal::MemoryPoolHandle pool)
    {
        // verify parameters
        if (!seal::is_valid_for(plain, context))
        {
            throw std::invalid_argument("plain is not valid for encryption parameters");
        }
        if (!plain.is_ntt_form())
        {
            throw std::invalid_argument("plain is not in NTT form");
        }
        if (!destination || stride == 0)
        {
            throw std::invalid_argument("destination cannot be null");
        }
        if (!pool)
        {
            throw std::invalid_argument("pool is uninitialized");
        }

        auto &context_data = *context.get_context_data(plain.parms_id());
        auto &parms = context_data.parms();
        uint32_t coeff_modulus_size = parms.coeff_modulus().size();
        uint32_t coeff_count = parms.poly_modulus_degree();
        uint32_t rns_poly_uint64_count = seal::util::mul_safe(coeff_count, coeff_modulus_size);

        if ((stride * vec_size) > coeff_count)
        {
            throw std::invalid_argument("stride, vec_len mismatches");
        }
        auto ntt_tables = context_data.small_ntt_tables();

        // Check that scale is positive and not too large
        if (plain.scale() <= 0 ||
            (static_cast<int>(log2(plain.scale())) >= context_data.total_coeff_modulus_bit_count()))
        {
            throw std::invalid_argument("scale out of bounds");
        }

        auto decryption_modulus = context_data.total_coeff_modulus();
        auto upper_half_threshold = context_data.upper_half_threshold();
        int logn = seal::util::get_power_of_two(coeff_count);

        // Quick sanity check
        if ((logn < 0) || (coeff_count < SEAL_POLY_MOD_DEGREE_MIN) || (coeff_count > SEAL_POLY_MOD_DEGREE_MAX))
        {
            throw std::logic_error("invalid parameters");
        }

        double inv_scale = double(1.0) / plain.scale();
        // Create mutable copy of input
        auto plain_copy(seal::util::allocate_uint(rns_poly_uint64_count, pool));
        seal::util::set_uint(plain.data(), rns_poly_uint64_count, plain_copy.get());

        // Transform each polynomial from NTT domain
        for (std::size_t i = 0; i < coeff_modulus_size; i++)
        {
            seal::util::inverse_ntt_negacyclic_harvey(plain_copy.get() + (i * coeff_count), ntt_tables[i]);
        }

        // CRT-compose the polynomial
        context_data.rns_tool()->base_q()->compose_array(plain_copy.get(), coeff_count, pool);

        // Create floating-point representations of the multi-precision integer coefficients
        double two_pow_64 = std::pow(2.0, 64);
        auto res(seal::util::allocate<double>(vec_size, pool));

        for (uint32_t i = 0; i < vec_size * stride; i += stride){
            uint32_t res_index = i / stride;
            res[res_index] = 0.0;
            if (seal::util::is_greater_than_or_equal_uint(
                    plain_copy.get() + (i * coeff_modulus_size), upper_half_threshold, coeff_modulus_size))
            {
                double scaled_two_pow_64 = inv_scale;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
                {
                    if (plain_copy[i * coeff_modulus_size + j] > decryption_modulus[j])
                    {
                        auto diff = plain_copy[i * coeff_modulus_size + j] - decryption_modulus[j];
                        res[res_index] += diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    }
                    else
                    {
                        auto diff = decryption_modulus[j] - plain_copy[i * coeff_modulus_size + j];
                        res[res_index] -= diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    }
                }
            }
            else
            {
                double scaled_two_pow_64 = inv_scale;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
                {
                    auto curr_coeff = plain_copy[i * coeff_modulus_size + j];
                    res[res_index] += curr_coeff ? static_cast<double>(curr_coeff) * scaled_two_pow_64 : 0.0;
                }
            }
        }

        std::copy_n(res.get(), vec_size, destination);
    }

    void encrypt(const seal::Plaintext &plain, const seal::PublicKey &pk, bool is_ntt_form,
                 std::string &ctxt, const seal::SEALContext &context)
    {
        seal::Ciphertext destination;
        seal::util::encrypt_zero_asymmetric(pk, context, plain.parms_id(), is_ntt_form, destination);
        add_plain_inplace(destination, plain, context);
        destination.scale() = plain.scale();
        std::ostringstream ostr;
        destination.save(ostr);
        ostr.str().swap(ctxt);
    }

    void encrypt(const seal::Plaintext &plain, const seal::PublicKey &pk, bool is_ntt_form,
                 seal::Ciphertext &ctxt, const seal::SEALContext &context){
        seal::util::encrypt_zero_asymmetric(pk, context, plain.parms_id(), is_ntt_form, ctxt);
        add_plain_inplace(ctxt, plain, context);
        ctxt.scale() = plain.scale();
    }

    void encrypt(const seal::Plaintext &plain, const seal::SecretKey &sk, bool is_ntt_form,
                 std::string &ctxt, const seal::SEALContext &context)
    {
        seal::Ciphertext destination;
        seal::util::encrypt_zero_symmetric(sk, context, plain.parms_id(), is_ntt_form, /* save_seed= */ true, destination);
        add_plain_inplace(destination, plain, context);
        destination.scale() = plain.scale();
        std::ostringstream ostr;
        destination.save(ostr);
        ostr.str().swap(ctxt);
    }

    void encrypt(const seal::Plaintext &plain, const seal::SecretKey &sk, bool is_ntt_form,
                 seal::Ciphertext &ctxt, const seal::SEALContext &context){

        seal::util::encrypt_zero_symmetric(sk, context, plain.parms_id(), is_ntt_form, false, ctxt);
        add_plain_inplace(ctxt, plain, context);
        ctxt.scale() = plain.scale();
    }

    void decrypt(const std::string &ctxt, const seal::SecretKey &sk, seal::Plaintext &plain, const seal::SEALContext &context)
    {
        seal::Ciphertext ct;
        ct.load(context, (seal::seal_byte *)ctxt.data(), ctxt.size());
        seal::Decryptor decryptor(context, sk);
        decryptor.decrypt(ct, plain);
    }

    void decrypt(const seal::Ciphertext &ctxt, const seal::SecretKey &sk, seal::Plaintext &plain, const seal::SEALContext &context)
    {
        seal::Decryptor decryptor(context, sk);
        decryptor.decrypt(ctxt, plain);
    }

    void scale_down(const seal::Plaintext &plain, uint64_t *result, size_t vec_len, DcdRole dcd_role, const AuxParms &aux_parms,
                    const seal::SEALContext &context, uint32_t mod_bits){
        if (result == nullptr || vec_len == 0)
            throw std::invalid_argument("invalid vector");

        auto &parms = context.get_context_data(plain.parms_id())->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = coeff_modulus.size();
        uint32_t poly_mod_degree = parms.poly_modulus_degree();
        int vec_len_bit_count = seal::util::get_significant_bit_count(vec_len - 1);
        size_t strides = poly_mod_degree >> vec_len_bit_count;

        uint64_t mask = (mod_bits == 64) ? static_cast<uint64_t>(-1) : (1ULL << mod_bits) - 1;
        uint64_t half_q0 = coeff_modulus[0].value() >> 1;
        uint64_t half_q1 = 0;
        if (coeff_modulus_size > 1)
            half_q1 = coeff_modulus[1].value() >> 1;
        if (coeff_modulus_size > 3)
            throw std::invalid_argument("DO NOT support now");

        auto rescale_two_mod = [&](uint64_t a_t, uint64_t a_0, uint64_t a_1) -> uint64_t
        {
            uint64_t rescale_q1_2_t, rescale_q1_q0, temp0, temp1;
            // rescale (a_t, a_0, a_1) to (b_t, b_0)
            temp0 = seal::util::add_uint_mod(a_1, half_q1, coeff_modulus[1]);
            temp1 = (a_t + half_q1) - temp0;
            // b_t = [(a_t + q1/2) - (a_1 + q1/2) mod q1] * q1^{-1} mod 2^t
            rescale_q1_2_t = (temp1 * aux_parms.inv_qi_mod_2_t_[1]) & mask;

            temp1 = seal::util::add_uint_mod(a_0, half_q1, coeff_modulus[0]);
            temp1 = seal::util::sub_uint_mod(temp1, temp0, coeff_modulus[0]);
            // b_0 = [(a_0 + q1/2) - (a_1 + q1/2) mod q1] * q1^{-1} mod q0
            rescale_q1_q0 = seal::util::multiply_uint_mod(temp1, aux_parms.inv_q1_mod_q0_, coeff_modulus[0]);

            // rescale (b_t, b_0) to get the final result
            // r = [(b_t + q0/2) - (b_0 + q0/2) mod q0] * q0^{-1} mod 2^t
            temp0 = seal::util::add_uint_mod(rescale_q1_q0, half_q0, coeff_modulus[0]);
            temp1 = rescale_q1_2_t + half_q0 - temp0;
            return (temp1 * aux_parms.inv_qi_mod_2_t_[0]) & mask;
        };

        if (dcd_role == Dcd_SCALED_STRIDE){
            // directly use float point number to round
            if (coeff_modulus_size == 1)
            {
                for (size_t i = 0; i < vec_len; ++i)
                {
                    size_t index = i * strides;
                    result[i] = double(plain[index]) * pow(2, mod_bits) / coeff_modulus[0].value();
                }
            }
            else if (coeff_modulus_size == 2)
            {
                for (size_t i = 0; i < vec_len; ++i)
                {
                    // represent the number as CRT form: (mod_2_t, mod_q0, mod_q1) = (a_t, a_0, a_1)
                    // a_t = (elt * 2^t) = 0 mod 2^t
                    size_t index = i * strides;
                    uint64_t a_0, a_1;
                    a_0 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[0], plain[index], coeff_modulus[0]);
                    a_1 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[1], plain[poly_mod_degree + index], coeff_modulus[1]);
                    result[i] = rescale_two_mod(0, a_0, a_1);
                }
            }
            else if (coeff_modulus_size == 3)
            {
                uint64_t half_q2 = coeff_modulus[2].value() >> 1;
                for (size_t i = 0; i < vec_len; ++i)
                {
                    // represent the number as CRT form (b_t, b_0, b_1, b_2), b_t = elt * 2^t = 0 mod 2^t
                    size_t index = i * strides;
                    uint64_t a_t, a_0, a_1, b_0, b_1, b_2, temp0, temp1;
                    b_0 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[0], plain[index], coeff_modulus[0]);
                    b_1 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[1], plain[index + poly_mod_degree], coeff_modulus[1]);
                    b_2 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[2], plain[index + poly_mod_degree * 2], coeff_modulus[2]);

                    // (b_2 + q2/2) mod q2
                    temp0 = seal::util::add_uint_mod(b_2, half_q2, coeff_modulus[2]);
                    a_t = ((half_q2 - temp0) * aux_parms.inv_qi_mod_2_t_[2]) & mask;

                    temp1 = seal::util::add_uint_mod(b_0, half_q2, coeff_modulus[0]);
                    temp1 = seal::util::sub_uint_mod(temp1, temp0, coeff_modulus[0]);
                    a_0 = seal::util::multiply_uint_mod(temp1, aux_parms.inv_q2_mod_q0_, coeff_modulus[0]);

                    temp1 = seal::util::add_uint_mod(b_1, half_q2, coeff_modulus[1]);
                    temp1 = seal::util::sub_uint_mod(temp1, temp0, coeff_modulus[1]);
                    a_1 = seal::util::multiply_uint_mod(temp1, aux_parms.inv_q2_mod_q1_, coeff_modulus[1]);
                    result[i] = rescale_two_mod(a_t, a_0, a_1);
                }
            }
        }else if (dcd_role == Dcd_SCALED_IN_ORDER) {
            // directly use float point number to round
            if (coeff_modulus_size == 1) {
                for (size_t i = 0; i < vec_len; ++i) {
                    result[i] = double(plain[i]) * pow(2, mod_bits) / coeff_modulus[0].value();
                }
            } else if (coeff_modulus_size == 2) {
                for (size_t i = 0; i < vec_len; ++i) {
                    // represent the number as CRT form: (mod_2_t, mod_q0, mod_q1) = (a_t, a_0, a_1)
                    // a_t = (elt * 2^t) = 0 mod 2^t
                    uint64_t a_0, a_1;
                    a_0 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[0], plain[i], coeff_modulus[0]);
                    a_1 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[1], plain[poly_mod_degree + i],
                                                        coeff_modulus[1]);
                    result[i] = rescale_two_mod(0, a_0, a_1);
                }
            } else // coeff_modulus_size == 3
            {
                uint64_t half_q2 = coeff_modulus[2].value() >> 1;
                for (size_t i = 0; i < vec_len; ++i) {
                    // represent the number as CRT form (b_t, b_0, b_1, b_2), b_t = elt * 2^t = 0 mod 2^t
                    size_t index = i;
                    uint64_t a_t, a_0, a_1, b_0, b_1, b_2, temp0, temp1;
                    b_0 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[0], plain[index], coeff_modulus[0]);
                    b_1 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[1],
                                                        plain[index + poly_mod_degree], coeff_modulus[1]);
                    b_2 = seal::util::multiply_uint_mod(aux_parms.pow_2_t_mod_qi_[2],
                                                        plain[index + poly_mod_degree * 2], coeff_modulus[2]);

                    // (b_2 + q2/2) mod q2
                    temp0 = seal::util::add_uint_mod(b_2, half_q2, coeff_modulus[2]);
                    a_t = ((half_q2 - temp0) * aux_parms.inv_qi_mod_2_t_[2]) & mask;

                    temp1 = seal::util::add_uint_mod(b_0, half_q2, coeff_modulus[0]);
                    temp1 = seal::util::sub_uint_mod(temp1, temp0, coeff_modulus[0]);
                    a_0 = seal::util::multiply_uint_mod(temp1, aux_parms.inv_q2_mod_q0_, coeff_modulus[0]);

                    temp1 = seal::util::add_uint_mod(b_1, half_q2, coeff_modulus[1]);
                    temp1 = seal::util::sub_uint_mod(temp1, temp0, coeff_modulus[1]);
                    a_1 = seal::util::multiply_uint_mod(temp1, aux_parms.inv_q2_mod_q1_, coeff_modulus[1]);
                    result[i] = rescale_two_mod(a_t, a_0, a_1);
                }
            }
        }
    }

    void transform_batched_to_ntt_inplace(std::vector<seal::Ciphertext> &ct, const seal::SEALContext &context, uint32_t threads)
    {
        uint32_t ct_num = ct.size();
        if (ct_num == 0)
            throw std::invalid_argument("empty ciphertext vector");

        // transform to ntt inplace
        auto tran2ntt = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                transform_to_ntt_inplace(ct[i], context);
            }
        };
        // multi-thread
        uint32_t thread_block = (ct_num + threads - 1) / threads;
        std::vector<std::thread> thread_pool(threads);
        for (size_t i = 0; i < threads; ++i){
            size_t bgn = i * thread_block;
            size_t end = std::min<size_t>(bgn + thread_block, ct_num);
            thread_pool[i] = std::thread(tran2ntt, bgn, end);
        }
        std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
    }
}
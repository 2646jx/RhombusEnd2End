#include "Rhombus/matmul_protocol.h"

using namespace std;
using namespace antchain;

int party, port;
shared_ptr<RhombusLinear> rhombus_;

int main(int argc, char **argv)
{
    sci::parse_party_and_port(argv, &party, &port);
    auto *io = new sci::NetIO(party == sci::ALICE ? nullptr : "127.0.0.1", port);

    uint32_t nthreads = 4;
    uint32_t ell = 37;
    uint32_t N = 8192;
    vector<int> coeff_mod{50, 50, 60};

    cout << nthreads << " threads" << endl;
    cout << "Ring bits: " << ell << endl;
#ifdef SEAL_USE_INTEL_HEXL
    cout << "AVX512 = ON" << endl;
#else
    cout << "AVX512 = OFF" << endl;
#endif    
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    rhombus_ = make_shared<RhombusLinear>(io, party, N, ell, coeff_mod);

    // the number of rows, columns
    uint32_t nr, nc;
    nr = 3479;
    nc = 246;

    cout << "Matrix: (" << nr << ", " << nc << ")" << endl;

    auto *mat = new uint64_t[nr * nc];
    auto *vec = new uint64_t[nc];
    auto *mv = new uint64_t[nr];

    if (party == 1){
        GenRandUintVector(mat, nr * nc, 15);
        GenRandUintVector(vec, nc, ell);
    }else{
        GenRandUintVector(vec, nc, ell);
    }

    cout << "StartComputation: " << endl;
    time_start = chrono::high_resolution_clock::now();
    uint64_t before_sent = rhombus_->IO()->counter;

    // matrix is held by server, vector is secret shared
    if (party == 1)
    {
        rhombus_->MatVecMul(nr, nc, mat, vec, mv);
    }
    else
    {
        rhombus_->MatVecMul(nr, nc, nullptr, vec, mv);
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    
    cout << "Elapsed time: " << time_diff.count() / 1000. << " ms" << endl;
    cout << "Sent: " << (rhombus_->IO()->counter - before_sent) / pow(2, 20) << " MB" << endl;

    // check
    if (party == 1) // SERVER
    {
        vector<uint64_t> vec_c(nc);
        vector<uint64_t> mv_c(nr);
        vector<uint64_t> vec_clear(nc);
        vector<uint64_t> mv_clear(nr);
        vector<uint64_t> mv_HE(nr);

        io->recv_data(vec_c.data(), sizeof(uint64_t) * nc);
        io->recv_data(mv_c.data(), sizeof(uint64_t) * nr);

        AddVecMod(vec, vec_c.data(), vec_clear.data(), nc, ell);
        MatVecMulMod(mat, vec_clear.data(), mv_clear.data(), nr, nc, ell);

        AddVecMod(mv, mv_c.data(), mv_HE.data(), nr, ell);
        uint64_t max_diff = 0;
        CompVector(mv_clear.data(), mv_HE.data(), nr, max_diff, ell);

        cout << "max error: " << max_diff << endl;
    }
    else // CLIENT
    {
        io->send_data(vec, sizeof(uint64_t) * nc);
        io->send_data(mv, sizeof(uint64_t) * nr);
    }

    delete[] mat;
    delete[] vec;
    delete[] mv;

    delete io;

}
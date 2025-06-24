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

    // set method
    // Note: only support PackRLWEs + V1 or Expand + V2 now, we will opensource other methods later.
    rhombus_->set_method(true);

    uint32_t n, m, k;
    n = 256;
    m = 256;
    k = 128;
    cout << "Matrix: (" << n << ", " << m << ", " << k << ")" << endl;

    auto *inputMatX = new uint64_t[n * m];
    auto *inputMatY = new uint64_t[m * k];
    auto *outMatXY = new uint64_t[n * k];

    if (party == 1) // SERVER / ALICE
    {
        GenRandUintVector(inputMatX, n * m, 15);
    }
    else
    {
        GenRandUintVector(inputMatY, m * k, 12);
    }

    cout << "StartComputation: " << endl;
    time_start = chrono::high_resolution_clock::now();
    uint64_t before_sent = rhombus_->IO()->counter;

    if (party == 1)
    {
        rhombus_->MatMul(n, m, k, inputMatX, nullptr, outMatXY);
    }
    else
    {
        rhombus_->MatMul(n, m, k, nullptr, inputMatY, outMatXY);
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    
    cout << "Elapsed time: " << time_diff.count() / 1000. << " ms" << endl;
    cout << "Sent: " << (rhombus_->IO()->counter - before_sent) / pow(2, 20) << " MB" << endl;

    // check correctness
    if (party == 1)
    {
        vector<uint64_t> MatY_c(m * k);
        vector<uint64_t> MatXY_c(n * k);

        io->recv_data(MatY_c.data(), sizeof(uint64_t) * m * k);
        io->recv_data(MatXY_c.data(), sizeof(uint64_t) * n * k);

        vector<uint64_t> MatXY(n * k);
        vector<uint64_t> MatXY_HE(n * k);
        MatMulMod(inputMatX, MatY_c.data(), MatXY.data(), n, m, k, ell);
        AddVecMod(outMatXY, MatXY_c.data(), MatXY_HE.data(), n * k, ell);

        uint64_t max_diff = 0;
        CompVector(MatXY.data(), MatXY_HE.data(), n * k, max_diff, ell);

        cout << "max error: " << max_diff << endl;

    }
    else
    {
        io->send_data(inputMatY, sizeof(uint64_t) * m * k);
        io->send_data(outMatXY, sizeof(uint64_t) * n * k);
    }

    delete[] inputMatX;
    delete[] inputMatY;
    delete[] outMatXY;

    delete io;

    return 0;
}
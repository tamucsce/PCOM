#include "KZG_verifier.h"
#include "enviroment.h"
#include "trusted_setup.h"
#include "KZG_commit.h"
#include "KZG_improved_prover.h"
#include "fft.h"
#include <algorithm>
#include <chrono>

mpz_class omega;
extern bn::Ec1 g1;
extern bn::Ec2 g2;

int main(int argc, char *argv[])
{
    int lg_degree, lg_players;
    sscanf(argv[1], "%d", &lg_degree);
    sscanf(argv[2], "%d", &lg_players);
    init_scratch_pad(std::max(1 << lg_players, 1 << lg_degree) * 2);
    init_enviroment();
    auto pp = trusted_setup_ec1(1 << lg_degree);
    auto commitment = commit_proof(lg_degree, pp.first);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto proofs = generate_proof(lg_degree, lg_players, pp.first);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    printf("Proof generation time: %f\n", time_span.count());
    omega = fastpow(root_of_unity, 1LL << (max_lg_order - lg_players));
    mpz_class x = 1;
    for(int i = 0; i < (1 << lg_players); ++i)
    {
        assert(verify(x, proofs[i].ans, proofs[i].pi, commitment, g1, g2, pp.second));
        x = x * omega % p;
    }
    return 0;
}
//
// Created by Tiancheng Xie on 3/13/2021.
//

#include "bn.h"
#include "gmpxx.h"
#include "gmp.h"
#include "enviroment.h"
#include <chrono>

bn::Ec1 fastmul(bn::Ec1 x, mpz_class y)
{
    bn::Ec1 x_arr[] = {g1 * 0, x, x + x, x + x + x};
    bn::Ec1 res = g1 * 0;
    int y_bits[256];
    int max_bit_pairs = 0;
    mpz_class four = 4;
    while(y != 0)
    {
        mpz_class mod4 = (y % four);
        y_bits[max_bit_pairs++] = mod4.get_si();
        y /= 4;
    }
    for(int i = max_bit_pairs - 1; i >= 0; --i)
    {
        res = res + x_arr[y_bits[i]];
        res = res + res;
        res = res + res;
    }
    return res;
}

const int maxn = 10000;

mpz_class rand_coef[maxn];
mie::Vuint rand_coef_vuint[maxn];

int main()
{
    init_enviroment();

    bn::Ec1 a = g1 * mie::Vuint(rand());
    bn::Ec1 res1 = g1 * 0;
    bn::Ec1 res2 = g1 * 0;

    for(int i = 0; i < maxn; ++i)
    {
        mpz_urandomm(rand_coef[i].get_mpz_t(), r_state, p.get_mpz_t());
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < maxn; ++i)
    {
        res1 = res1 + fastmul(a, rand_coef[i]);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    printf("res1 generation time: %f\n", time_span.count());

    t0 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < maxn; ++i)
    {
        res2 = res2 + a * rand_coef[i];
    }
    t1 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    printf("res2 generation time: %f\n", time_span.count());

    assert(res1 == res2);

    return 0;
}

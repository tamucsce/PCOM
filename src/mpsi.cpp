/*
    a tiny sample of optimal ate pairing
*/

#include <array>
#include <cstdlib>
#include <fstream>
#include <future>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <time.h>
#include <vector>
#include <optional>
#include <ligero/util/timer.hpp>

#include "bn.h"
#include "test_point.hpp"
// /#include "fft.h"
//#include "enviroment.h"
#include <chrono>
#include <omp.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <ligero/tcp/server.hpp>
#include <ligero/tcp/client.hpp>
#include <ligero/network/serialize.hpp>

// ============================================================
/*  Include this to enable vector serialization */
#include <boost/serialization/vector.hpp>
// ============================================================

#include "utils/csv.hpp"

#define G 16
int num_bits = 254;
#define MIE_ATE_USE_GMP 1


using namespace ligero::net;
using namespace ligero::util;
using namespace std::chrono_literals;
using namespace std::string_literals;

using serializer_t = serializer<boost_portable_binary>;
using tcp_server = basic_tcp_server<serializer_t>;
using tcp_client = basic_tcp_client<serializer_t>;

using namespace std;
using namespace bn;

namespace boost::serialization {


//serialization for mpz_class



//serialization for Ec1
template <typename Archive>
void serialize(Archive& ar, mie::Fp& fp, const unsigned int version) {
    for (auto i = 0; i < Fp::N; i++) {
        ar & fp[i];
    }
}


template <typename Archive>
void serialize(Archive& ar, Ec1& ec, const unsigned int version) {
    ar & ec.p[0];
    ar & ec.p[1];
    ar & ec.p[2];
}

//serialization for Ec2
template <typename Archive>
void serialize(Archive& ar, bn::Fp2& fp2, const unsigned int version) {
    //for (auto i = 0; i < Fp::N; i++) {
        ar & fp2.a_;
        ar & fp2.b_;

    //}
}


template <typename Archive>
void serialize(Archive& ar, Ec2& ec, const unsigned int version) {
    ar & ec.p[0];
    ar & ec.p[1];
    ar & ec.p[2];
}


template <typename Archive>
void serialize(Archive& ar, bn::Fp6& fp, const unsigned int version) {
    ar & fp.a_;
    ar & fp.b_;
    ar & fp.c_;
}

template <typename Archive>
void serialize(Archive& ar, bn::Fp12& fp, const unsigned int version) {
    ar & fp.a_;
    ar & fp.b_;
}

}

mpz_class zero(0);
mpz_class exponetiate(mpz_class base, mpz_class pow);
void printCipherPoly(Ec1* poly, int size);

// ------------------------------------------------------------
struct pubkey {
    pubkey() { }
    pubkey(const Ec1& u, const Ec1& v, const Ec1& h) : data_({u, v, h}) { }

    const Ec1& operator[](size_t idx) const { return data_[idx]; }
    Ec1& operator[](size_t idx) { return data_[idx]; }
    Ec1 get(int idx){return data_[idx]; }
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        ar & data_[0];
        ar & data_[1];
        ar & data_[2];
    }

    pubkey& operator*=(const mpz_class& exponent) {
        data_[0] = data_[0] * exponent;
        data_[1] = data_[1] * exponent;
        data_[2] = data_[2] * exponent;
        return *this;
    }


    pubkey& operator+=(const pubkey& other) {
        data_[0] += other.data_[0];
        data_[1] += other.data_[1];
        data_[2] += other.data_[2];
        return *this;
    }

public:
    std::array<Ec1, 3> data_;
};

pubkey operator+(pubkey k, const pubkey& k1) {
    return k += k1;
}

using ciphertext = pubkey;
// ------------------------------------------------------------

ciphertext encrypt(const pubkey& key, const Ec1& g, mpz_class message) {
    mpz_class mpr = mpz_class((Param::r).toString());
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));

    mpz_class r = rand.get_z_bits(400) % mpr;
    mpz_class s = rand.get_z_bits(400) % mpr;
    //mpz_class r(0);
    //mpz_class s(0);
    mpz_class rs = r + s;

    ciphertext cipher{
        key[0] * r,
        key[1] * s,
        g * rs + key[2] * message
    };
    DEBUG << "Encryption is u = " << cipher[0] << ", v = " << cipher[1] << ", w = " << cipher[2] << std::endl;
    return cipher;
}

std::vector<ciphertext> rerandomize(const pubkey& key, const Ec1& g, const std::vector<ciphertext>& ctxt_vec ) {

    mpz_class mpr = mpz_class((Param::r).toString());
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));

    std::vector<ciphertext> reenc;
    for (auto i=0; i < ctxt_vec.size(); i++) {
      mpz_class temp = rand.get_z_bits(400) % mpr;
      ciphertext ctemp(ctxt_vec[i]);
      ctemp *= temp;
      ctemp += encrypt(key,g,mpz_class(0));
      reenc.emplace_back(ctemp);
    }
    return reenc;
}
std::vector<ciphertext> homorphic_polyeval(const pubkey& key, const Ec1& g, const std::vector<mpz_class>& input,std::vector<ciphertext>& enc_poly ) {
   

    std::vector<ciphertext> enc_eval;
    #pragma omp parallel for
    //wenxuan NEW VERSION
    int num_eval = input.size();
    Ec1* eval1 = new Ec1[num_eval];
    mpz_class idx;
    for (int j = 0; j < num_eval; ++j)
    {
        mpz_class t_mpz = input[j];
        Ec1 eval0_ = enc_poly[0][0]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
        Ec1 eval1_ = enc_poly[0][1]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
        Ec1 eval2_ = enc_poly[0][2]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
        for (int i = 1; i < enc_poly.size(); i++)
        {
            idx = mpz_class(i);
            eval0_ += enc_poly[i][0]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));
            eval1_ += enc_poly[i][1]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));
            eval2_ += enc_poly[i][2]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));

        }
        ciphertext eval(eval0_,eval1_,eval2_);
        enc_eval.emplace_back(eval);

    }
    return enc_eval;
    //wenxuan


    //OLD version
    for (auto i=0; i < input.size(); i++) {
      ciphertext eval_sum(enc_poly[0]);
      mpz_class point(input[i]);

      #pragma omp parallel for
      
      for( auto j=1; j < enc_poly.size(); j++) {
       //std::cout << "point[" << j << "] = " << point << std::endl;
        ciphertext temp(enc_poly[j]);
        temp *= point;
        eval_sum += temp;
        point *= input[i];
            
      }
      enc_eval.emplace_back(eval_sum);
    }
    return enc_eval;
}

ciphertext homorphic_polyeval(const pubkey& key, const Ec1& g, const mpz_class input,std::vector<ciphertext>& enc_poly ) {
    //cout<<"inputs: "<<input<<endl;

    // std::vector<ciphertext> enc_eval;
    
    //Wenxuan's NEW Version
    mpz_class t_mpz = input;
    mpz_class idx;

    Ec1 eval0_ = enc_poly[0][0]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
    Ec1 eval1_ = enc_poly[0][1]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
    Ec1 eval2_ = enc_poly[0][2]*(mie::Vuint(exponetiate(t_mpz,zero).get_str().c_str()));
    for (int i = 1; i < enc_poly.size(); i++)
    {
        idx = mpz_class(i);
        eval0_ += enc_poly[i][0]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));
        eval1_ += enc_poly[i][1]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));
        eval2_ += enc_poly[i][2]*(mie::Vuint(exponetiate(t_mpz,idx).get_str().c_str()));

    }
    ciphertext eval(eval0_,eval1_,eval2_);
    return eval;


    //OLD Version
    ciphertext eval_sum(enc_poly[0]);
    mpz_class point(input);

    //#pragma omp parallel for
    for( auto j=1; j < enc_poly.size(); j++) {
      ciphertext temp(enc_poly[j]);
      temp *= point;
      eval_sum += temp;
      point *= input;
    }
    return eval_sum;
}

static int errNum = 0;
unsigned long int seed;
mpz_class mod("21888242871839275222246405745257275088548364400416034343698204186575808495617");
mpz_class root_of_unity;

gmp_randstate_t r_state; //included in environment.h
mpz_class p;             //included in environment.h
Ec2* v;          //pk
int m,n,f_len,g_len;
Ec1 g1;
Ec2 g2;

Ec1* w;         //ck2
Ec2 h1;         //ck3
Ec2 h2;         //ck3

mpz_class* b;    //batched eval vectors

//initialize all the proofs
Ec1* gy;
vector<Fp12> proof_L;
vector<Fp12> proof_R;
vector<Ec1*> proof_g_fL;
vector<Ec1*> proof_g_fR;

vector<Ec1> p_;
vector<Ec1> r_;

Ec1 final_proof_1;
mpz_class final_proof_2;

//debug proof
vector<Ec1*> proof_A;
vector<mpz_class*> proof_b;
Ec1 final_proof;
//challenges x
vector<mpz_class> chals;
mpz_class* r; //for batching

//selector matrices & eval matrices
mpz_class** T;//used to compute b
mpz_class** matrix_T;//for phase 2
mpz_class** matrix_A;
mpz_class** matrix_B;
mpz_class** matrix_C;
mpz_class** matrix_D;

//Prover's broadcast msgs


//Prover's broadcast msgs
Ec1* g_fL;
Ec1* g_fR;
Ec1* g_fR_cy;
int newsize;
Fp12 com;

mpz_class* eval_points;
//tools to calculate runtime
float eval_time = 0;
float pairing_time = 0;
float folding_time = 0;

float folding_time_veri = 0;
float eval_time_veri = 0;
int test_round;


float Ec1_fft_time = 0;
//calculate f0/f1
Ec1* f0_;
Ec1* f1_;
int num_eval_pt; // = deg(f)
//Lookup Table
Ec1** table;

//some repetitive values
//mpz_class zero = mod;
mpz_class one = mpz_class(1);
mie::Vuint neg("21888242871839275222246405745257275088548364400416034343698204186575808495616");
mie::Vuint ling("21888242871839275222246405745257275088548364400416034343698204186575808495617");
mie::Vuint er("2");
mpz_class neg_mpz = mpz_class(-1);

Ec1 g1_p;

//================Phase 2==========================//
//phase 2 pk
Ec1* v_2;          //ck1
Ec1 cku;         //ck3
//Phase 2 initialize all proofs
vector<Ec1> proof_L_phase2;
vector<Ec1> proof_R_phase2;
//debug proof
mpz_class final_proof_1_phase2;
mpz_class final_proof_2_phase2;
//matrix commitment Ti, Xi
Ec1* com_Ti;

Ec1* com_Xi;
Ec1* com_Yi;
Ec1* com_Zi;
Ec1* com_Fi;


//Phase2 2 randomness
mpz_class* s;
mpz_class* u;

//Phase 2 X,Y,Z
mpz_class** matrix_X;
mpz_class** matrix_Y;
mpz_class** matrix_Z;
mpz_class** matrix_F;

//Phase2 Sa, Tu, Xu
mpz_class* Tu;

mpz_class* Sa;
mpz_class* Xu;

mpz_class* Sb;
mpz_class* Yu;

mpz_class* Sc;
mpz_class* Zu;

mpz_class* Sd;
mpz_class* Fu;
//phase2 com_l, com_r
Ec1 com_l_X;
Ec1 com_r_X;

Ec1 com_l_Y;
Ec1 com_r_Y;

Ec1 com_l_Z;
Ec1 com_r_Z;

Ec1 com_l_F;
Ec1 com_r_F;

Ec1 com_phase2_X;
Ec1 com_phase2_Y;
Ec1 com_phase2_Z;
Ec1 com_phase2_F;
//================Phase 3==========================//
//phase 3 randomness
mpz_class *R;
mpz_class W;
mpz_class x_phase3;

//phase3 Li, Ri
mpz_class** Li;
mpz_class** Ri;
//phase3 com_Li, com_Ri
Ec1* com_Li;
Ec1* com_Ri;
//phase3 commitment
Ec1 com_d;

//phase3 proof
vector<Ec1> proof_coml;
vector<Ec1> proof_comr;
int proof_size = 0;
int commitment_size = 0 ;
int total_communication = 0;
float comTi_time; 
// Client / Server (Party P1) code
//
//


//evaluation results
vector<Ec1> eval_results;



mie::Vuint invert(mpz_class x){ // calculate x^-1 mod p

    mpz_class x_inverse;
    mpz_invert (x_inverse.get_mpz_t(),x.get_mpz_t(), mod.get_mpz_t());
    mie::Vuint x_inv(x_inverse.get_str().c_str());//x_negOne = x^(-1) mod p

    return x_inv;
}
mpz_class invert_mpz_class(mpz_class x){ // calculate x^-1 mod p

    mpz_class x_inverse;
    mpz_invert (x_inverse.get_mpz_t(),x.get_mpz_t(), mod.get_mpz_t());
    mpz_class x_inv(x_inverse.get_str().c_str());//x_negOne = x^(-1) mod p

    return x_inv;
}

bool isPowerOfTwo(int n)
{
   if(n==0)
   return false;
   return (ceil(log2(n)) == floor(log2(n)));
}

vector<mpz_class> generate_phase1_random_challenge(int rounds){
    vector<mpz_class> x(rounds);

    mpz_class x_each;
    for(int i=0;i<rounds;i++){
        mpz_urandomm(x_each.get_mpz_t(),r_state,p.get_mpz_t());
        x[i] = x_each;
    }
    return x;
}
void init_proof(int rounds){
     //phase 1
    proof_L.resize(rounds);
    proof_R.resize(rounds);
    proof_g_fL.resize(rounds);
    proof_g_fR.resize(rounds);

    proof_A.resize(rounds);
    proof_b.resize(rounds);

    //phase2
    proof_L_phase2.resize(rounds);
    proof_R_phase2.resize(rounds);

    //phase3
    proof_coml.resize(rounds+1);
    proof_comr.resize(rounds+1);
}
mpz_class exponetiate(mpz_class base, mpz_class pow){
        mpz_class result;
        //mpz_class mod("16798108731015832284940804142231733909759579603404752749028378864165570215949");
        mpz_powm(result.get_mpz_t(),base.get_mpz_t(),pow.get_mpz_t(),mod.get_mpz_t());
        return result;  
}
Ec1* ec1AddEc1(Ec1* a, Ec1* b, int size){
    Ec1* ret = new Ec1[size];
    for (int i = 0; i < size; ++i)
    {
        ret[i] = (a[i] + b[i]);
    }
    return ret;
}
Ec1* mpzTimesEc1Vector(mpz_class a, Ec1* b, int size){
    Ec1* ret = new Ec1[size];
    for (int i = 0; i < size; ++i)
    {
        ret[i] = b[i] * (mie::Vuint(a.get_str().c_str()));
    }
    return ret;
}
mpz_class* mpzhadamardmpz(mpz_class* a, mpz_class* b, int size){
    mpz_class* ret = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        ret[i] = (a[i] * b[i])%mod;
    }
    return ret;
}

mpz_class* mpzAddmpz(mpz_class* a, mpz_class* b, int size){
    mpz_class* ret = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        ret[i] = (a[i] + b[i])%mod;
    }
    return ret;
}
mpz_class* mpzTimesMpzVector(mpz_class a, mpz_class* b, int size){
    mpz_class* ret = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        ret[i] = (a * b[i])%mod+mod;
    }
    return ret;
}
void init_phase2_randomness(int size, int num_eval){
    //init s
    int s_size = size - 1;

    mpz_class s_;
    s = new mpz_class[s_size];
    mpz_urandomm(s_.get_mpz_t(),r_state,p.get_mpz_t());
    for (int i = 0; i < s_size; ++i)
    {
        s[i] = exponetiate(s_,i);
    }

    //init u
    int u_size = num_eval;
    u = new mpz_class[u_size];
    for (int i = 0; i < u_size; ++i)
    {
        mpz_class u_each;
        mpz_urandomm(u_each.get_mpz_t(),r_state,p.get_mpz_t());
        u[i] = u_each;
    }
}
void init_phase3_randomness(int size){
    //init R
    int R_size = size - 1;
    mpz_class r_;
    R = new mpz_class[R_size];
    mpz_urandomm(r_.get_mpz_t(),r_state,p.get_mpz_t());
    for (int i = 0; i < R_size; ++i)
    {
        R[i] = exponetiate(r_,i);
    }

    mpz_urandomm(W.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(x_phase3.get_mpz_t(),r_state,p.get_mpz_t());

}

mpz_class * pad_zero(mpz_class *poly, int size_own, int &size){
        // cout<<"size "<<size<<endl;
        // cout<<"size_own "<<size_own<<endl;

    if(isPowerOfTwo(size)==0){
        int n = ceil(log2(size));
        int newsize = pow(2,n);

        mpz_class *poly_new = new mpz_class[newsize];
        for (int i = 0; i < newsize; ++i)
        {

            if(i<size_own){
                poly_new[i] = poly[i];

            }
            else{
                poly_new[i] = p;

            }
        }
        size = newsize;
        //delete [] poly;
        return poly_new;
    }
    else{
        mpz_class *poly_new = new mpz_class[size];
        for (int i = 0; i < size; ++i)
        {
            if(i<size_own){
                poly_new[i] = poly[i];
            }
            else{
                poly_new[i] = p;
            }
        }
        //delete [] poly;
        return poly_new;
    }
    return poly;
}
Ec1 * pad_zero_Ec1(Ec1 *poly, int size_own, int &size,Ec1 g1){


    if(isPowerOfTwo(size)==0){
        int n = ceil(log2(size));
        int newsize = pow(2,n);
        Ec1 *poly_new = new Ec1[newsize];
        for (int i = 0; i < newsize; ++i)
        {
            if(i<size_own){
                poly_new[i] = poly[i];
            }
            else{
                poly_new[i] = g1*mie::Vuint(0);
            }
        }
        size = newsize;
        //delete[] poly;
        return poly_new;
    }
    else{
        Ec1 *poly_new = new Ec1[size];
        for (int i = 0; i < size; ++i)
        {
            if(i<size_own){
                poly_new[i] = poly[i];
            }
            else{
                poly_new[i] = g1*mie::Vuint(0);
            }
        }
        //delete[] poly;
        return poly_new;
    }
    return poly;
}
///new things 22/03/2021
mpz_class* __dst[3];
bn::Ec1* __dst_ec[3];
mpz_class* twiddle_factor;
bn::Ec2* twiddle_factor_ec2;

void init_scratch_pad(int order)
{
    __dst[0] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[1] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[2] = new mpz_class[order];
    __dst_ec[0] = new bn::Ec1[order];
    __dst_ec[1] = new bn::Ec1[order];
    __dst_ec[2] = new bn::Ec1[order];
    // twiddle_factor = new mpz_class[order];
    // twiddle_factor_ec2 = new bn::Ec2[order];
}
template<class exp_T>
mpz_class fastpow(mpz_class x, exp_T y)
{

    mpz_class ret = 1;
    exp_T two = 2;
    while(y != 0)
    {
        if(y % 2 == 1)
        {
            ret = ret * x % p;
        }
        x = x * x % p;
        y = y / two;
    }
    return ret;
}
void FFT(const mpz_class *coefficients, int coef_len, int order, mpz_class root_of_unity, mpz_class *result)
{

    mpz_class rot_mul[62];
    //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
    //In sake of both memory and time efficiency, use the non-recursive version
    int lg_order = -1;
    rot_mul[0] = root_of_unity;
    for(int i = 0; i < 62; ++i)
    {
        if(i > 0)
            rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1] % p;
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    //we can merge both cases, but I just don't want to do so since it's easy to make mistake
    if(lg_coef > lg_order)
    {
        assert(false);
    }
    else
    {
        //initialize leaves
        int blk_sz = (order / coef_len);
        for(int j = 0; j < blk_sz; ++j)
        {
            for(int i = 0; i < coef_len; ++i)
            {       
                __dst[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
            }
        }

        mpz_class *x_arr = new mpz_class[1 << lg_order];

        {
            for(int dep = lg_coef - 1; dep >= 0; --dep)
            {
                int blk_size = 1 << (lg_order - dep);
                int half_blk_size = blk_size >> 1;

                int cur = dep & 1;
                int pre = cur ^ 1;
                mpz_class x = 1;
                x_arr[0] = 1;
                for(int j = 1; j < blk_size; ++j)
                    x_arr[j] = x_arr[j - 1] * rot_mul[dep] % p;
                for(int k = 0; k < blk_size / 2; ++k)
                {
                    int double_k = (k) & (half_blk_size - 1);
                    for(int j = 0; j < (1 << dep); ++j)
                    {
                        mpz_class l_value = __dst[pre][double_k << (dep + 1) | j], r_value = x_arr[k] * __dst[pre][double_k << (dep + 1) | (1 << dep) | j];
                        __dst[cur][k << dep | j] = (l_value + r_value) % p;
                        __dst[cur][(k + blk_size / 2) << dep | j] = ((l_value - r_value) % p + p) % p;
                    }
                }
            }
        }
        delete[] x_arr;
    }

    for(int i = 0; i < order; ++i)
        result[i] = __dst[0][i]%mod;
}
void iFFT(mpz_class *evaluations, int coef_len, int order, mpz_class root_of_unity, mpz_class *dst)
{
    if(coef_len > order)
    {
        //more coefficient than evaluation
        fprintf(stderr, "Warning, Request do inverse fft with inefficient number of evaluations.");
        fprintf(stderr, "Will construct a polynomial with less order than required.");
        coef_len = order;
    }

    //assume coef_len <= order

    //subsample evalutions

    mpz_class *sub_eval;
    bool need_free = false;
    if(coef_len != order)
    {
        need_free = true;
        sub_eval = new mpz_class[coef_len];
        for(int i = 0; i < coef_len; ++i)
        {
            sub_eval[i] = evaluations[i * (order / coef_len)];
        }
    }
    else
        sub_eval = evaluations;

    mpz_class new_rou = 1;
    for(int i = 0; i < order / coef_len; ++i)
        new_rou = new_rou * root_of_unity % p;
    order = coef_len;

    mpz_class inv_rou = 1, tmp = new_rou;
    int lg_order = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    for(int i = 0; i < lg_order; ++i)
    {
        inv_rou = inv_rou * tmp % p;
        tmp = tmp * tmp % p;
    }

    FFT(sub_eval, order, coef_len, inv_rou, dst);

    if(need_free)
        delete[] sub_eval;

    mpz_class inv_n = fastpow<mpz_class>(mpz_class(order), p - 2);

    for(int i = 0; i < coef_len; ++i)
    {
        dst[i] = dst[i] * inv_n;
    }

}
void FFT_cipher(const Ec1 *coefficients, int coef_len, int order, mpz_class root_of_unity, bn::Ec1 *result)
{
    mpz_class rot_mul[62];
    //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
    //In sake of both memory and time efficiency, use the non-recursive version
    int lg_order = -1;
    rot_mul[0] = root_of_unity;
    for(int i = 0; i < 62; ++i)
    {
        if(i > 0)
            rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1] % p;
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    //we can merge both cases, but I just don't want to do so since it's easy to make mistake
    if(lg_coef > lg_order)
    {
        assert(false);
    }
    else
    {
        //initialize leaves
        int blk_sz = (order / coef_len);
        for(int j = 0; j < blk_sz; ++j)
        {
            for(int i = 0; i < coef_len; ++i)
            {
                __dst_ec[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
            }
        }

        mpz_class *x_arr = new mpz_class[1 << lg_order];
        {
            for(int dep = lg_coef - 1; dep >= 0; --dep)
            {
                int blk_size = 1 << (lg_order - dep);
                int half_blk_size = blk_size >> 1;
                int cur = dep & 1;
                int pre = cur ^ 1;

                mpz_class x = 1;
                x_arr[0] = 1;
                for(int j = 1; j < blk_size; ++j)
                    x_arr[j] = x_arr[j - 1] * rot_mul[dep] % p;
                for(int k = 0; k < blk_size / 2; ++k)
                {
                    int double_k = (k) & (half_blk_size - 1);
                    for(int j = 0; j < (1 << dep); ++j)
                    {

                        Ec1 l_value = __dst_ec[pre][double_k << (dep + 1) | j];
                        Ec1 r_value = __dst_ec[pre][double_k << (dep + 1) | (1 << dep) | j];
                        r_value = r_value * mie::Vuint(x_arr[k].get_str().c_str());               //i modify this line
                        __dst_ec[cur][k << dep | j] = (l_value + r_value);
                        __dst_ec[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                    }
                }
            }
        }
        delete[] x_arr;
    }

    for(int i = 0; i < order; ++i)
        result[i] = __dst_ec[0][i];
}

void iFFT_cipher(Ec1 *evaluations, int coef_len, int order, mpz_class root_of_unity, bn::Ec1 *dst)
{
    if(coef_len > order)
    {
        //more coefficient than evaluation
        fprintf(stderr, "Warning, Request do inverse fft with inefficient number of evaluations.");
        fprintf(stderr, "Will construct a polynomial with less order than required.");
        coef_len = order;
    }

    //assume coef_len <= order

    //subsample evalutions

    Ec1 *sub_eval;
    bool need_free = false;
    if(coef_len != order)
    {
        need_free = true;
        sub_eval = new Ec1[coef_len];
        for(int i = 0; i < coef_len; ++i)
        {
            sub_eval[i] = evaluations[i * (order / coef_len)];
        }
    }
    else
        sub_eval = evaluations;

    mpz_class new_rou = 1;
    for(int i = 0; i < order / coef_len; ++i)
        new_rou = new_rou * root_of_unity % p;
    order = coef_len;

    mpz_class inv_rou = 1, tmp = new_rou;
    int lg_order = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    for(int i = 0; i < lg_order; ++i)
    {
        inv_rou = inv_rou * tmp % p;
        tmp = tmp * tmp % p;
    }

    FFT_cipher(sub_eval, order, coef_len, inv_rou, dst);
    if(need_free)
        delete [] sub_eval;

    mpz_class inv_n = fastpow<mpz_class>(mpz_class(order), p - 2);

    for(int i = 0; i < coef_len; ++i)
    {
        dst[i] = dst[i] * mie::Vuint(inv_n.get_str().c_str());
    }

    //wenxuan delete

}

void printPoly(mpz_class* poly, int size){
    for (int i = 0; i < size; ++i)
    {
        cout<<poly[i]%mod<<endl;
    }
}
void printCipherPoly(Ec1* poly, int size){
    for (int i = 0; i < size; i++)
    {
        cout<<poly[i]<<endl;
    }
}
void printCipher2Poly(Ec2* poly, int size){
    for (int i = 0; i < size; i++)
    {
        cout<<poly[i]<<endl;
    }
}
void printCipherPoly_vec(vector<Ec1> poly, int size){
    for (int i = 0; i < size; i++)
    {
        //cout<<poly[i]<<endl;
    }
}

void init_selector_matrices(int num_eval){
    int matrix_size = num_eval -1; //only select 1 less of the total 

    //initialize Matrix A
    matrix_A = new mpz_class*[matrix_size];
    for (int i = 0; i < matrix_size; ++i)
    {
        matrix_A[i] = new mpz_class[num_eval];
        for (int j = 0; j < num_eval; ++j)
        {
            if(i == j ){
                matrix_A[i][j] = 1;
            }
            else{
                matrix_A[i][j] = 0;
            }   

        }
    }
    //initialize Matrix B
    matrix_B = new mpz_class*[matrix_size];
    for (int i = 0; i < matrix_size; ++i)
    {
        matrix_B[i] = new mpz_class[num_eval];
        for (int j = 0; j < num_eval; ++j)
        {
            if(j==1){
                matrix_B[i][j] = 1;
            }
            else{
                matrix_B[i][j] = 0;
            }   
            //cout<<matrix_B[i][j]<<"";   
        }
            //cout<<""<<endl;   

    }
             //cout<<"====================="<<endl;   

    //initialize Matrix C
    matrix_C = new mpz_class*[matrix_size];
    for (int i = 0; i < matrix_size; ++i)
    {
        matrix_C[i] = new mpz_class[num_eval];
        for (int j = 0; j < num_eval; ++j)
        {
            if(j == i+1 ){
                matrix_C[i][j] = 1;
            }
            else{
                matrix_C[i][j] = 0;
            }
            //cout<<matrix_C[i][j]<<"";   
        }
            //cout<<""<<endl;   
    }

    //initialize Matrix D
    matrix_D = new mpz_class*[matrix_size];
    for (int i = 0; i < matrix_size; ++i)
    {
        matrix_D[i] = new mpz_class[num_eval];
        for (int j = 0; j < num_eval; ++j)
        {
            if(j==0){
                matrix_D[i][j] = 1;
            }
            else{
                matrix_D[i][j] = 0;
            }   
            //cout<<matrix_B[i][j]<<"";   
        }
            //cout<<""<<endl;   

    }

}


Ec1 ped_enc(mpz_class* a, mpz_class* b, mpz_class d, Ec1* ck1, Ec1* ck2, int size){
    Ec1 com_ped;
    com_ped = ck1[0] * (mie::Vuint(a[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_ped = com_ped + ck1[i] * (mie::Vuint(a[i].get_str().c_str()));
    }

    for (int i = 0; i < size; ++i)
    {
        com_ped = com_ped + ck2[i] * (mie::Vuint(b[i].get_str().c_str()));
    }
    d = d%mod;
    com_ped = com_ped + cku*(mie::Vuint(d.get_str().c_str()));

    return com_ped;
}

Fp12 bbm_enc(Ec1* A, mpz_class* b, Ec1 d, Ec2* ck1, Ec1* ck2, int size){
    Fp12 com_bbm;
    opt_atePairing(com_bbm,  ck1[0],A[0]);
    for (int i = 1; i < size; i++)
    {
        Fp12 e;
        opt_atePairing(e, ck1[i],A[i]);
        com_bbm*=e;
    }
    Ec1 com_b = ck2[0]*(mie::Vuint(b[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_b = com_b + ck2[i]*(mie::Vuint(b[i].get_str().c_str()));
    }
    // Ec1 com_b = com_Ti[0];
    // for (int i = 1; i < size; ++i)
    // {
    //     com_b = com_b + com_Ti[i];
    // }


    Fp12 e_;
    opt_atePairing(e_, h1,com_b);
    com_bbm *=e_;
    opt_atePairing(e_, h2, d);
    com_bbm *=e_;
    return com_bbm;

}
Fp12 compute_com_phase1(Ec1* A, mpz_class* b, Ec1 d, Ec2* ck1, Ec1* ck2, mpz_class* r, int size, int num_eval_pts){
    Fp12 com_bbm;
    opt_atePairing(com_bbm,  ck1[0],A[0]);
    for (int i = 1; i < size; i++)
    {
        Fp12 e;
        opt_atePairing(e, ck1[i],A[i]);
        com_bbm*=e;
    }
    // Ec1 com_b = ck2[0]*(mie::Vuint(b[0].get_str().c_str()));
    // for (int i = 1; i < size; ++i)
    // {
    //     com_b = com_b + ck2[i]*(mie::Vuint(b[i].get_str().c_str()));
    // }
    Ec1 com_b = com_Ti[0]*(mie::Vuint(r[0].get_str().c_str()));
    for (int i = 1; i < num_eval_pts; ++i)
    {
        com_b = com_b + com_Ti[i]*(mie::Vuint(r[i].get_str().c_str()));
    }


    Fp12 e_;
    opt_atePairing(e_, h1,com_b);
    com_bbm *=e_;
    opt_atePairing(e_, h2, d);
    com_bbm *=e_;
    return com_bbm;

}
Ec1* multiplyEc1_mpz(Ec1* poly1,int size_poly1, mpz_class* poly2, int size_poly2, Ec1 g1){

    //int size_poly1 = get_size(poly1);
    //int size_poly2 = get_size(poly2);


    //cout<<"size_poly1: "<<size_poly1<<endl;
    //cout<<"size_poly2: "<<size_poly2<<endl;

    int size = max(size_poly1,size_poly2);
    poly1 = pad_zero_Ec1(poly1, size_poly1,size,g1);
    poly2 = pad_zero(poly2, size_poly2,size);
    
    //cout<<size<<endl;
    //cout<<"size_poly1: "<<size_poly1<<endl;
    //cout<<"size_poly2: "<<size_poly2<<endl;
    

    //int size = max(size_poly1,size_poly2);
    //FFT on poly1
    init_scratch_pad(size*2);
    mpz_class rou = fastpow(root_of_unity, (1LL << (28)) / (size * 2));
    //cout<<"here"<<endl;
    Ec1* eval1 = new Ec1[size*2];
    FFT_cipher(poly1, size, size*2, rou,eval1);
    //cout<<"here"<<endl;



    //FFT on poly2
    mpz_class* eval2 = new mpz_class[size*2];
    FFT(poly2, size, size*2, rou,eval2);
    //multiply 2 evals
    Ec1 *eval = new Ec1[size*2];

    #pragma omp parallel for
    for(int i=0;i<size*2;i++){
        // eval[i] = eval1[i] * (mie::Vuint(eval2[i].get_str().c_str()));
        eval[i] = eval1[i] * eval2[i];
    }


    //iFFT on eval
    Ec1 *coeff = new Ec1[size*2];
    iFFT_cipher(eval,size*2,size*2, rou, coeff);

    //delete is necessary
    delete[] poly1;
    delete[] poly2;
    delete[] eval1;
    delete[] eval2;
    delete[] eval;

    delete[] __dst[0];
    delete[] __dst[1];
    delete[] __dst[2];
    delete[] __dst_ec[0];
    delete[] __dst_ec[1];
    delete[] __dst_ec[2];
    return coeff;
}
mpz_class* multiply2Poly(mpz_class* poly1,int size_poly1, mpz_class* poly2,int size_poly2){

    // int size_poly1 = get_size(poly1);
    // int size_poly2 = get_size(poly2);


    // cout<<"size_poly1: "<<size_poly1<<endl;
    // cout<<"size_poly2: "<<size_poly2<<endl;

    int size = max(size_poly1,size_poly2);

    poly1 = pad_zero(poly1, size_poly1,size);
    poly2 = pad_zero(poly2, size_poly2,size);
    // cout<<"size_poly1_new: "<<size<<endl;
    // cout<<"size_poly2_new: "<<size<<endl;  
    
    //int size = max(size_poly1,size_poly2);
    //FFT on poly1
    init_scratch_pad(size*2);
    mpz_class rou = fastpow(root_of_unity, (1LL << (28)) / (size * 2));

    mpz_class* eval1 = new mpz_class[size*2];
    FFT(poly1, size, size*2, rou,eval1);

    //FFT on poly2
    mpz_class* eval2 = new mpz_class[size*2];
    FFT(poly2, size, size*2, rou,eval2);

    //multiply 2 evals

    mpz_class *eval = new mpz_class[size*2];
    for(int i=0;i<size*2;i++){
        eval[i] = eval1[i] * eval2[i];
    }
    //iFFT on eval
    mpz_class *coeff = new mpz_class[size*2];
    iFFT(eval,size*2,size*2, rou, coeff);
    //delete is necessary
    delete[] poly1;
    delete[] poly2;
    delete[] eval1;
    delete[] eval2;
    delete[] eval;
    delete[] __dst[0];
    delete[] __dst[1];
    delete[] __dst[2];
    delete[] __dst_ec[0];
    delete[] __dst_ec[1];
    delete[] __dst_ec[2];
    return coeff;
}

mpz_class* multiplyPoly(mpz_class** polys,int num_poly){

    if(num_poly==1){
        return polys[0];
    }
    if(num_poly == 2){
        return multiply2Poly(polys[0],2,polys[1],2);
    }
    else{
        int size_left = num_poly/2;
        int size_right = num_poly - (num_poly/2);
        mpz_class* l_half[size_left];
        mpz_class* r_half[size_right];


        for(int i=0;i<size_left;i++){
            l_half[i] = polys[i];
        }
        for(int i=0;i<size_right;i++){
            r_half[i] = polys[i+size_left];
        }
        
        return multiply2Poly(multiplyPoly(l_half,size_left),size_left+1, multiplyPoly(r_half,size_right),size_right+1);
    }
}

mpz_class* getCoeff(const mpz_class* roots, int m){
    mpz_class **polys = new mpz_class*[m];
    for (int i = 0; i < m; ++i)
    {   
        polys[i] = new mpz_class[2];
        polys[i][0] = p-roots[i];
        polys[i][1] = mpz_class(1);

    }
    mpz_class *coeffs = multiplyPoly(polys,m);
    // for (int i = 0; i < m; ++i)
    // {
    //     delete[] polys[i];
    // }
    
    delete [] polys;
    return coeffs;
}

Ec1 innerEc1Mpz(Ec1* vec1, mpz_class* vec2, int size){
    //cout<<vec2[1]<<endl;
    Ec1 innerProduct = vec1[0] * (mie::Vuint(vec2[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        innerProduct = innerProduct + vec1[i] * (mie::Vuint(vec2[i].get_str().c_str()));
    }
    return innerProduct;
}
mpz_class innerMpzMpz(mpz_class* vec1, mpz_class* vec2, int size){

    mpz_class innerProduct = vec1[0] * vec2[0];
    for (int i = 1; i < size; ++i)
    {
        innerProduct = (innerProduct + vec1[i] * vec2[i])%mod;
    }
    return innerProduct;
}

mpz_class* subtract2Poly(mpz_class* poly1, int size1, mpz_class* poly2, int size2){


    int size_big = max(size1,size2);
    int size_small = min(size1,size2);

    mpz_class *poly = new mpz_class[size_big];
    for (int i = 0; i < size_big; ++i)
    {
        if(i<size_small){
            poly[i] = poly1[i] - poly2[i];
        }
        else{
            if(size1>size2){
                poly[i] = poly1[i];
            }
            else{
                poly[i] = p-poly2[i];
            }
        }
    }
    delete [] poly1;
    delete [] poly2;
    return poly;

}

mpz_class* poly_mod_x(mpz_class* poly, int power){
    //int size = get_size(poly)+20;
    mpz_class* poly_new = new mpz_class[power];
    for (int i = 0; i < power; ++i)
    {
        poly_new[i] = poly[i]%mod;

    }
    delete[]poly;
    return poly_new;


}
mpz_class* flip_Poly_Coeff(mpz_class* poly, int size){//it takes a original polynomial and the size of the polynomial

    mpz_class *poly_flip = new mpz_class[size];
    int idx = size;
    for(int i = 0;i<size;i++){
        poly_flip[i] = poly[idx-1];
        idx = idx-1;
    }
    //delete[] poly;
    return poly_flip;
}
Ec1* flip_Poly_Ec1_Coeff(Ec1* poly, int size_to_flip, int size_total){//it takes a original polynomial and the size of the polynomial

    Ec1 *poly_flip = new Ec1[size_total];
    int idx = size_to_flip;
    for(int i = 0;i<size_total;i++){
        if(i<size_to_flip){
            poly_flip[i] = poly[idx-1];
            idx = idx-1;
        }else{
            poly_flip[i] = poly[i];

        }
    }
    //delete[] poly;
    return poly_flip;
}

mpz_class* compute_h(mpz_class* g_flip, int g_len,int n, int m){// deg(f) = n, deg(g) = m
    int l = ceil(log2(n-m+1));
    //int h_size = (n-m+1);
    int h_size = pow(2,l);

    //cout<<"l: "<<l<<endl;
    //cout<<"n-m+1: "<<n-m+1<<endl;

    mpz_class*h = new mpz_class[1];
    h[0] = mpz_class(1);

    mpz_class* two = new mpz_class[1];
    two[0] = mpz_class(2);
    

    mpz_class* h_next;

    int h_size_i = 1;
    for(int  i = 0;i<l;i++){
    
        int h_square_size = h_size_i+h_size_i-1;
        h_next =  subtract2Poly(multiply2Poly(two,1,h,h_size_i),h_size_i, multiply2Poly(g_flip,g_len, multiply2Poly(h,h_size_i,h,h_size_i),h_square_size), (g_len+h_square_size-1));
        h_size_i = pow(2,i+1);
        h_next = poly_mod_x(h_next,pow(2,i+1));
        h = h_next;
    }



    for (int i = 0; i < h_size; i++)
    {
        h[i] = h[i]%p;
    }
    delete [] two;
    return h;

}


void keygen_public2(int size, Ec1 g1, Ec2 g2){
    //construct v

    mpz_class rand;
    //mpz_urandomm(rand.get_mpz_t(), r_state,rand.get_mpz_t());
    mpz_class alpha_mpz;
    mpz_class beta_mpz;

    //generating alpha and beta for public keys
    mpz_urandomm(alpha_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(beta_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mie::Vuint alpha(alpha_mpz.get_str().c_str());

    v = new Ec2[size];

    for(int i=0;i<size;i++){
        mpz_class beta_each = exponetiate(beta_mpz,2*i);
        mie::Vuint beta(beta_each.get_str().c_str());

        v[i] = g2 *beta;    
    }
}
void keygen_public(int size, Ec1 g1, Ec2 g2){
    //construct v,w,h1,h2

    mpz_class rand;
    //mpz_urandomm(rand.get_mpz_t(), r_state,rand.get_mpz_t());
    mpz_class alpha_mpz;
    mpz_class beta_mpz;
    mpz_class omega_mpz;
    mpz_class c1_mpz;
    mpz_class c2_mpz;


    //generating alpha and beta for public keys
    mpz_urandomm(alpha_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(beta_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(omega_mpz.get_mpz_t(),r_state,p.get_mpz_t());

    mpz_urandomm(c1_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(c2_mpz.get_mpz_t(),r_state,p.get_mpz_t());

    //mie::Vuint alpha(alpha_mpz.get_str().c_str());

    v = new Ec2[size];
    v_2 = new Ec1[size];
    w = new Ec1[size];
    for(int i=0;i<size;i++){
        mpz_class beta_each = exponetiate(alpha_mpz,2*i);
        mpz_class alpha_each = exponetiate(beta_mpz,2*i);
 
        mie::Vuint beta(beta_each.get_str().c_str());
        mie::Vuint alpha(alpha_each.get_str().c_str());

        v[i] = g2 *beta;
        v_2[i] = g1 *beta;

        w[i] = g1 *alpha;    
    
    }
    mie::Vuint c1(c1_mpz.get_str().c_str());
    mie::Vuint c2(c2_mpz.get_str().c_str()); 
    h1 = g2*c1;
    h2 = g2*c2;


    //phase 2
    mie::Vuint omega(omega_mpz.get_str().c_str());
    cku = g1*omega;

}
Ec1* dividePoly(Ec1* A, mpz_class *roots,int A_len, int deg_A,int num_roots, string side, Ec1 g1){
    Ec1* A_flip = flip_Poly_Ec1_Coeff(A,deg_A+1,A_len);
    clock_t f0 = clock();
    clock_t t = clock();


    g_len = num_roots+1;
    mpz_class*g = new mpz_class[g_len];
    g = getCoeff(roots, num_roots);


    #pragma omp parallel for
    for (int i = 0; i < g_len; ++i)
    {
       g[i] = g[i]%p;
    }
      

    //compute g_flip = g flip coeff
    mpz_class *g_flip = new mpz_class[g_len];
    g_flip = flip_Poly_Coeff(g, g_len);
 

    m = num_roots;//deg(g)
    n = 2*m;//deg(f)
    n = deg_A;

    //try to use the same notation used from lec 5 by yupeng
    t = clock();
    //printf("%s\n","here" );
    int l = ceil(log2(n-m+1));
    int h_size = pow(2,l);
    //int h_size = n-m+1;
    mpz_class * h = new mpz_class[h_size];
    h = compute_h(g_flip,g_len,n,m);
   


    //=========================================//
    //==========calcualte r (f0 or f1) ========//
    //=========================================//
    
    //compute q_flip = A_flip * h;   A_flip is cipher, but h is plain text
    int q_size = n-m+1;




    t = clock();
    Ec1* q_flip =  multiplyEc1_mpz(A_flip,A_len,h,h_size,g1);
    Ec1_fft_time += clock()-t;


    float q_flip_time = clock() - t;


    //compute q = flip(q_flip);
    Ec1 *q =flip_Poly_Ec1_Coeff(q_flip,q_size,q_size);

   

   //compute r = f-qg
    t=clock();
    int qg_len = (q_size-1)+(g_len-1)+1;
    // /cout<<qg_len<<endl;
    t = clock();
    Ec1 * qg = multiplyEc1_mpz(q,q_size,g, g_len,g1);
    Ec1_fft_time += clock()-t;

    //compute remainder r = A - qg
    int r_size = A_len/2+1;

    Ec1* r = new Ec1[r_size];
    for (int i = 0; i < r_size; ++i)
    {
        r[i] = A[i] + (qg[i]*mie::Vuint(neg));

    }

    //cleaning up 
    delete [] g;
    delete [] g_flip;
    delete [] h;
    delete [] qg;
    delete [] q;

    return r;

   

}
Ec1 multi_scalar2(Ec1* g, mpz_class* b, int size){

 
    // preprocess lookup table
    int num_table = size/G;
    Ec1 com_mul_all = g1_p;
    //compute_init_time += clock() - t_init;

    for (int i = 0; i < num_table; ++i)
    {
        int idx = G*i;
        Ec1 com_mul = g1_p;
    
        for (int j = num_bits; j >=0; j--)
        {
            //clock_t t_compute_selector = clock();
            int selector = 0;
            for (int k = 0; k < G; ++k)
            {
                selector+=mpz_tstbit(b[idx+k].get_mpz_t(),j);
                if (k!=G-1)
                {
                     selector*=2;
                }
            }
            //compute_selector_time += clock() - t_compute_selector;



            //clock_t t_look_table = clock();
            if (j == 0)
            {
                com_mul =  com_mul + table[i][selector];
            }else{
                com_mul =  (com_mul + table[i][selector]) * er;
            }
            //look_table_time += clock() - t_look_table;

        }
        com_mul_all = com_mul_all +com_mul;
    }
    return com_mul_all;

}
void compute_com_Ti(int num_eval, int size){

    // ===================== testing =================
    float newtime;
    clock_t t_multiscalar = clock();
    com_Ti = new Ec1[num_eval];
    for (int i = 0; i < num_eval; ++i)
    {
        mpz_class* T_each = T[i];
        com_Ti[i] = multi_scalar2(w,T_each,size);
    }
    newtime = clock() - t_multiscalar;

    // =======================ending testing==========
    cout<<"multiscalar exp time: "<<(double)(newtime)/CLOCKS_PER_SEC<<"s\n";
}
void compute_com_Xi_com_Yi_com_Zi_com_Fi(int size, int num_eval_pts){
     //commit Xi
    com_Xi = new Ec1[num_eval_pts];
    for (int i = 0; i < num_eval_pts; ++i)
    {
        mpz_class* X_each = new mpz_class[size-1];
        for (int j = 0; j < size-1; ++j)
        {
            X_each[j] = matrix_X[j][i];
        }
        com_Xi[i] = multi_scalar2(w,X_each,size);
    }
    //commit Yi
    com_Yi = new Ec1[num_eval_pts];
    for (int i = 0; i < num_eval_pts; ++i)
    {
        mpz_class* Y_each = new mpz_class[size-1];
        for (int j = 0; j < size-1; ++j)
        {
            Y_each[j] = matrix_Y[j][i];

        }
        com_Yi[i] = multi_scalar2(w,Y_each,size);
    }

    //commit Zi
    com_Zi = new Ec1[num_eval_pts];
    for (int i = 0; i < num_eval_pts; ++i)
    {
        mpz_class* Z_each = new mpz_class[size-1];
        for (int j = 0; j < size-1; ++j)
        {
            Z_each[j] = matrix_Z[j][i]+mod%mod;

        }
        com_Zi[i] = multi_scalar2(w,Z_each,size);
    }

    //commit Fi
    com_Fi = new Ec1[num_eval_pts];
    for (int i = 0; i < num_eval_pts; ++i)
    {
        mpz_class* F_each = new mpz_class[size-1];
        for (int j = 0; j < size-1; ++j)
        {
            F_each[j] = matrix_F[j][i]+mod%mod;

        }
        com_Fi[i] = multi_scalar2(w,F_each,size);
    }
}
void compute_L_and_R(mpz_class* L,mpz_class* R,mpz_class chal_x,int size,int idx){
    if(idx == 0){//L_X, R_X
        //compute L_X, R_X
        for (int i = 0; i < size; ++i)
        {
            L[i] = (Sa[i] + invert_mpz_class(chal_x)*s[i])%mod;
        }

        for (int i = 0; i < size; ++i)
        {
            R[i] = (Tu[i] - chal_x*Xu[i])%mod;
        }

    }
    else if(idx == 1){
         //compute L_Y, R_Y
        for (int i = 0; i < size; ++i)
        {
            L[i] = (Sb[i] + invert_mpz_class(chal_x)*s[i])%mod;
        }

        for (int i = 0; i < size; ++i)
        {
            R[i] = (Tu[i] - chal_x*Yu[i])%mod;
        }

    }
    else if(idx == 2){
        
         //compute L_Z, R_Z
        for (int i = 0; i < size; ++i)
        {
            L[i] = (Sc[i] + invert_mpz_class(chal_x)*s[i])%mod;
        }
        for (int i = 0; i < size; ++i)
        {
            R[i] = (Tu[i] - chal_x*Zu[i])%mod;
        }

    }
    else{
        //compute L_F, R_F
        for (int i = 0; i < size; ++i)
        {
            L[i] = (Sd[i] + invert_mpz_class(chal_x)*s[i])%mod;
        }
        for (int i = 0; i < size; ++i)
        {
            R[i] = (Tu[i] - chal_x*Fu[i])%mod;
        }
    }
    
}

void compute_l_r_commitments(Ec1 cku, int size, int num_eval,mpz_class x){
    //compute com_l_X, com_r_X
    mpz_class l,r;
    size = size-1;
    int size_s = size + 1;
    mpz_class* Xu_minus = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        Xu_minus[i] = (neg_mpz*Xu[i])%mod;
    }

    l = innerMpzMpz(Sa, Xu_minus,size)+mod;
    r = innerMpzMpz(s, Tu,size) + mod;

    com_l_X = cku * (mie::Vuint(l.get_str().c_str()));
    com_r_X = cku * (mie::Vuint(r.get_str().c_str()));

    //compute com_l_Y, com_r_Y
    mpz_class* Yu_minus = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        Yu_minus[i] = (neg_mpz*Yu[i])%mod;
    }
    l = innerMpzMpz(Sb, Yu_minus,size)+mod;
    //r = innerMpzMpz(s, Tu,size) + mod;

    com_l_Y = cku * (mie::Vuint(l.get_str().c_str()));
    com_r_Y = cku * (mie::Vuint(r.get_str().c_str()));


    //compute com_l_Z, com_r_Z
    mpz_class* Zu_minus = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        Zu_minus[i] = (neg_mpz*Zu[i])%mod;
    }
    l = innerMpzMpz(Sc, Zu_minus,size)+mod;
    //r = innerMpzMpz(s, Tu,size) + mod;

    com_l_Z = cku * (mie::Vuint(l.get_str().c_str()));
    com_r_Z = cku * (mie::Vuint(r.get_str().c_str()));


    //compute com_l_F, com_r_F
    mpz_class* Fu_minus = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        Fu_minus[i] = (neg_mpz*Fu[i])%mod;
    }
    l = innerMpzMpz(Sd, Fu_minus,size)+mod;
    //r = innerMpzMpz(s, Tu,size) + mod;

    com_l_F = cku * (mie::Vuint(l.get_str().c_str()));
    com_r_F = cku * (mie::Vuint(r.get_str().c_str()));
    // mpz_class lr = l*x + invert_mpz_class(x)*r;
    // Ec1 test =  cku*(mie::Vuint(lr.get_str().c_str()));
    // cout<<"test"<<endl;
    // cout<<test<<endl;
}
void compute_Sa_Xu_Tu_Sb_Yu_Sc_Zu_Sd_Fu(int size, int num_eval){

    //compute Tu
    Tu = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        mpz_class Tu_each = innerMpzMpz(matrix_T[i],u,num_eval);
        Tu[i] = Tu_each;
    }
    
    //compute Sa
    int size_s = size - 1;
    Sa = new mpz_class[size];

    for (int i = 0; i < size; ++i)
    {
        mpz_class Sa_each = 0;
        for (int j = 0; j < size_s; ++j)
        {
            Sa_each = Sa_each + s[j] * matrix_A[j][i];
        }
        Sa[i] = Sa_each;
    }
 

    //compute Xu
    int t_size = size;
    int Xu_size = size-1;
    Xu = new mpz_class[Xu_size];
    for (int i = 0; i < Xu_size; ++i)
    {
        Xu[i] = innerMpzMpz(matrix_X[i],u,num_eval);
    }

    //compute Sb
    Sb = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        mpz_class Sb_each = 0;
        for (int j = 0; j < size_s; ++j)
        {
            Sb_each = Sb_each + s[j] * matrix_B[j][i];
        }
        Sb[i] = Sb_each;

    }

    //compute Yu
    int Yu_size = size-1;
    Yu = new mpz_class[Yu_size];
    for (int i = 0; i < Yu_size; ++i)
    {
        Yu[i] = innerMpzMpz(matrix_Y[i],u,num_eval);
    }

    //compute Sc
    Sc = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        mpz_class Sc_each = 0;
        for (int j = 0; j < size_s; ++j)
        {
            Sc_each = (Sc_each + s[j] * matrix_C[j][i]);
        }
        Sc[i] = Sc_each;
    }

    //compute Zu
    int Zu_size = size-1;
    Zu = new mpz_class[Zu_size];
    for (int i = 0; i < Zu_size; ++i)
    {
        Zu[i] = innerMpzMpz(matrix_Z[i],u,num_eval);
    }

    //compute Sd
    Sd = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        mpz_class Sd_each = 0;
        for (int j = 0; j < size_s; ++j)
        {
            Sd_each = (Sd_each + s[j] * matrix_D[j][i]);
        }
        Sd[i] = Sd_each;
    }
    //compute Fu
    int Fu_size = size-1;
    Fu = new mpz_class[Fu_size];
    for (int i = 0; i < Fu_size; ++i)
    {
        Fu[i] = innerMpzMpz(matrix_F[i],u,num_eval);
    }

}


Ec1 verifier_compute_com_phase2_X(int size,int num_eval,mie::Vuint x,  mie::Vuint x_inv , mpz_class x_inv_mpz){
    Ec1 com_X;
    Ec1 com_L, com_R;
    //compute com_L
    mpz_class * L_phase2 = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = Sa[i] + x_inv_mpz*s[i];
    }


    com_L = v_2[0] * (mie::Vuint(L_phase2[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_L = com_L + v_2[i] * (mie::Vuint(L_phase2[i].get_str().c_str()));
    }

    // cout<<"com_L_X"<<endl;
    // cout<<com_L<<endl;
    //1 com_Ti * u
    Ec1 com_Tu = (com_Ti[0] * (mie::Vuint(u[0].get_str().c_str())));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Tu = com_Tu + (com_Ti[i] * (mie::Vuint(u[i].get_str().c_str())));
    }

    //test Xu commitment
    //1 com_Xi *u
    Ec1 com_Xu = com_Xi[0] * (mie::Vuint(u[0].get_str().c_str()));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Xu = com_Xu + (com_Xi[i] * (mie::Vuint(u[i].get_str().c_str())));
    }


    com_R = com_Tu - (com_Xu * x);
    
    Ec1 com_u2 = (com_l_X*x) + (com_r_X*x_inv);//verifier compute u2's commitment

    com_X = com_L + com_R + com_u2;

    return com_X;
}
Ec1 verifier_compute_com_phase2_Y(int size,int num_eval,mie::Vuint x,  mie::Vuint x_inv , mpz_class x_inv_mpz){
    Ec1 com_Y;
    Ec1 com_L, com_R;
    //compute com_L
    mpz_class * L_phase2 = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = Sb[i] + x_inv_mpz*s[i];
    }

    com_L = v_2[0] * (mie::Vuint(L_phase2[0].get_str().c_str()));
    for (int i = 1; i < size-1; ++i)
    {
        com_L = com_L + v_2[i] * (mie::Vuint(L_phase2[i].get_str().c_str()));
    }

    //1 com_Ti * u
    Ec1 com_Tu = (com_Ti[0] * (mie::Vuint(u[0].get_str().c_str())));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Tu = com_Tu + (com_Ti[i] * (mie::Vuint(u[i].get_str().c_str())));
    }
    //test Yu commitment
    //1 com_Yi *u
    Ec1 com_Yu = com_Yi[0] * (mie::Vuint(u[0].get_str().c_str()));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Yu = com_Yu + (com_Yi[i] * (mie::Vuint(u[i].get_str().c_str())));
    }

    com_R = com_Tu - (com_Yu * x);
  
    Ec1 com_u2 = (com_l_Y*x) + (com_r_Y*x_inv);//verifier compute u2's commitment

    com_Y = com_L + com_R + com_u2;

    return com_Y;
}
Ec1 verifier_compute_com_phase2_Z(int size,int num_eval, mie::Vuint x,  mie::Vuint x_inv , mpz_class x_inv_mpz){

    Ec1 com_Z;
    Ec1 com_L, com_R;
    //compute com_L
    mpz_class * L_phase2 = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = (Sc[i] + x_inv_mpz*s[i])%mod;
    }
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = L_phase2[i]+mod;
        
    }
 

    com_L = v_2[0] * (mie::Vuint(L_phase2[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_L = com_L + v_2[i] * (mie::Vuint(L_phase2[i].get_str().c_str()));
    }

    //1 com_Ti * u
    Ec1 com_Tu = (com_Ti[0] * (mie::Vuint(u[0].get_str().c_str())));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Tu = com_Tu + (com_Ti[i] * (mie::Vuint(u[i].get_str().c_str())));
    }
    //test Zu commitment
    //1 com_Zi *u
    Ec1 com_Zu = com_Zi[0] * (mie::Vuint(u[0].get_str().c_str()));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Zu = com_Zu + (com_Zi[i] * (mie::Vuint(u[i].get_str().c_str())));
    }

    com_R = com_Tu - (com_Zu * x);

    Ec1 com_u2 = (com_l_Z*x) + (com_r_Z*x_inv);//verifier compute u2's commitment

    com_Z = com_L + com_R + com_u2;

    return com_Z;
}
Ec1 verifier_compute_com_phase2_F(int size,int num_eval, mie::Vuint x,  mie::Vuint x_inv , mpz_class x_inv_mpz){

    Ec1 com_F;
    Ec1 com_L, com_R;
    //compute com_L
    mpz_class * L_phase2 = new mpz_class[size];
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = (Sd[i] + x_inv_mpz*s[i])%mod;
    }
    for (int i = 0; i < size; ++i)
    {
        L_phase2[i] = L_phase2[i]+mod;
        
    }
 

    com_L = v_2[0] * (mie::Vuint(L_phase2[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_L = com_L + v_2[i] * (mie::Vuint(L_phase2[i].get_str().c_str()));
    }

    //1 com_Ti * u
    Ec1 com_Tu = (com_Ti[0] * (mie::Vuint(u[0].get_str().c_str())));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Tu = com_Tu + (com_Ti[i] * (mie::Vuint(u[i].get_str().c_str())));
    }
    //test Fu commitment
    //1 com_Fi *u
    Ec1 com_Fu = com_Fi[0] * (mie::Vuint(u[0].get_str().c_str()));
    for (int i = 1; i < num_eval; ++i)
    {
        com_Fu = com_Fu + (com_Fi[i] * (mie::Vuint(u[i].get_str().c_str())));
    }

    com_R = com_Tu - (com_Fu * x);

    Ec1 com_u2 = (com_l_F*x) + (com_r_F*x_inv);//verifier compute u2's commitment

    com_F = com_L + com_R + com_u2;

    return com_F;
}
Ec1 verifier_compute_com_phase3(Ec1* w_prime, mpz_class* L, mpz_class* R, int size){
    Ec1 com_phase3 = w[0] * (mie::Vuint(L[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_phase3 = com_phase3 + w[i] * (mie::Vuint(L[i].get_str().c_str()));
    }
    for (int i = 0; i < size; ++i)
    {
        com_phase3 = com_phase3 + w_prime[i] * (mie::Vuint(R[i].get_str().c_str()));
    }
    com_phase3 = com_phase3 + com_d;
    return com_phase3;
}


void compute_Li_Ri_com_Li_com_Ri(Ec1* w_prime,int size, int num_eval){
    Li = new mpz_class*[2*size];
    Ri = new mpz_class*[2*size];

    com_Li = new Ec1[2*size];
    com_Ri = new Ec1[2*size];

    //compute -R  and  com(-R) for only once
    mpz_class* R_neg = new mpz_class[size-1];
    R_neg =  mpzTimesMpzVector(neg_mpz, R, size-1);

    Ec1 com_R_neg = innerEc1Mpz(w_prime,R_neg,size-1); 
    for (int i = 0; i < size; ++i)
    {
        if(i<num_eval){
            mpz_class* Xi = new mpz_class[size-1];
            mpz_class* Yi = new mpz_class[size-1];
            mpz_class* Zi = new mpz_class[size-1];

            for (int j = 0; j < size-1; ++j)
            {
                Xi[j] = matrix_X[j][i];
                Yi[j] = matrix_Y[j][i];
                Zi[j] = matrix_Z[j][i];

            }
            mpz_class Wi = exponetiate(W,i);
            Li[i] = mpzTimesMpzVector(Wi,Xi,size-1);
            Li[i+size] = mpzTimesMpzVector(Wi,Zi,size-1);
            Ri[i] = mpzhadamardmpz(Yi, R, size-1);
            Ri[i+size] = R_neg;

            //compute each commitment
            com_Li[i] = com_Xi[i] * (mie::Vuint(Wi.get_str().c_str()));
            com_Li[i+size] = com_Zi[i] * (mie::Vuint(Wi.get_str().c_str()));
            com_Ri[i] = com_Yi[i];
            com_Ri[i+size] = com_R_neg;
        }
        else{

            mpz_class* zeros = new mpz_class[size-1];
            Ec1 g_zero = w[0] * ling;

            for (int j = 0; j < size-1; ++j)
            {
                zeros[j] = zero;
            }
            Li[i] = zeros;
            Li[i+size] = zeros;
            Ri[i] = zeros;
            Ri[i+size] = zeros;

            //compute each commitment
            com_Li[i] =g_zero;
            com_Li[i+size] = g_zero;
            com_Ri[i] = g_zero;
            com_Ri[i+size] = g_zero;
        }
        
    }
    

}

void printEc1Vector(vector<Ec1> v){
    for (int i = 0; i < v.size(); ++i)
    {
        cout<<v[i]<<", ";
    }
}
void printVector(vector<int> v){
    for (int i = 0; i < v.size(); ++i)
    {
        cout<<v[i]<<", ";
    }
    printf("%s\n","" );
}

void compute_X_Y_Z_F(int num_eval, int size){//prover compute the three matrices
    //compute X = A*T 
    int matrix_size =size - 1;
    matrix_X = new mpz_class*[matrix_size];
    matrix_Y = new mpz_class*[matrix_size];
    matrix_Z = new mpz_class*[matrix_size];
    matrix_F = new mpz_class*[matrix_size];

   
    for (int i = 0; i < matrix_size; ++i)
    {
        matrix_X[i] = new mpz_class[num_eval];
        matrix_Y[i] = new mpz_class[num_eval];
        matrix_Z[i] = new mpz_class[num_eval];
        matrix_F[i] = new mpz_class[num_eval];

    }
    for (int i = 0; i < num_eval; ++i)
    {
        mpz_class* T_each = T[i];
      
        for (int j = 0; j < matrix_size; ++j)
        {
           
            matrix_X[j][i] = T_each[j]%mod+mod;  
            matrix_Y[j][i] = T_each[1]%mod+mod;              
            matrix_Z[j][i] = T_each[j+1]%mod+mod;              
            matrix_F[j][i] = T_each[0]%mod+mod;

        }
        

    }
    return;
    
}
void compute_look_up_table2(Ec1* g, int size){
    mie::Vuint yi(one.get_str().c_str());
    int num_table = size/G;
    table = new Ec1*[num_table];
    for (int i = 0; i < num_table; ++i)
    {
        int num_per_group_up = pow(2,G);
        table[i] = new Ec1[num_per_group_up];    
        int idx = G*i;

        int possibilities = pow(2,G);
        for (int j = 0; j < possibilities; ++j)
        {
            string s = bitset<G>(j).to_string();
            table[i][j] = g[0]*ling;
            for (int k = 0; k <G; k++)
            {
                string bit;
                bit.push_back(s[k]);
                table[i][j] = table[i][j] + g[idx+k] *mie::Vuint(bit);
            }
           
        }
    }
}

void multiPointsEval(Ec1* A, mpz_class* points, int A_len,int deg_A,int num_points, int *indices, Ec1* output, Ec1 g1){

    // if(test_round!=1 and num_points==num_eval_pt/2){
    //     printf("%s\n","====================HALF WAY THRU!!!!!==================" );
    // }

    test_round++;
    if(num_points==1){
        Ec1 eval = A[0]*(mie::Vuint(exponetiate(points[0],zero).get_str().c_str()));
        mpz_class idx;
        for (int i = 1; i < A_len; i++)
        {
            idx = mpz_class(i);
            eval += A[i]*(mie::Vuint(exponetiate(points[0],idx).get_str().c_str()));
        }
        
        output[indices[0]] = eval;

        return;
    }
    else{
        mpz_class* roots_l = new mpz_class[num_points/2];
        mpz_class* roots_r = new mpz_class[num_points/2];
        int indices_len = num_points/2;
        int* indices_l = new int [indices_len];
        int* indices_r = new int [indices_len];


        int idx = 0;
        for (int i = 0; i < num_points; i++)
        {
            if(i<num_points/2){
                roots_l[i] = points[i];
                indices_l[i] = indices[i];

            }
            else{
                roots_r[idx] = points[i];
                indices_r[idx] = indices[i];
                idx =  idx+1;
            }
        }


        //compute f0/f1
        num_points = num_points/2;


        Ec1* f0 = dividePoly(A,roots_l,A_len,deg_A,num_points,"left",g1);
        // exit(0);
        multiPointsEval(f0,roots_l,A_len/2+1,deg_A/2,num_points,indices_l,output,g1);
        delete [] f0;
        delete [] roots_l;
        delete [] indices_l;
        //exit(1);

        Ec1* f1 = dividePoly(A,roots_r,A_len,deg_A,num_points,"right",g1);
        multiPointsEval(f1,roots_r,A_len/2+1,deg_A/2,num_points,indices_r,output,g1);
        delete [] f1;
        delete [] roots_r;
        delete [] indices_r;

       
        return;

    }
}


std::vector<Ec1*> results_thread(2);
Ec1* multiPointsEval_naive_thread(Ec1* A, vector<mpz_class*> eval_points, int A_size, int deg_A, int num_eval_pt, Ec1 g1, int thread_idx, Ec1* results, int round, int side){
    Ec1 g_fR_cy_;
    mpz_class idx;
    for(int k = 0; k<num_eval_pt;k++){
        //cout<<"thread "<<thread_idx<<" eval_points: "<<eval_points[thread_idx][0]<<endl;
        Ec1 result = A[0]*(mie::Vuint(exponetiate(eval_points[thread_idx][k],zero).get_str().c_str()));
        int right_degree = deg_A;
        //Ec1 g_fR_cy_;
        if(round==0 && side ==1){
            g_fR_cy_ = result*(mie::Vuint(exponetiate(eval_points[thread_idx][k],right_degree).get_str().c_str()));
        }
        for (int l = 1; l < A_size; l++)
        {
            idx = mpz_class(l);
            result += A[l]*(mie::Vuint(exponetiate(eval_points[thread_idx][k],idx).get_str().c_str()));
            if(round==0){
                g_fR_cy_ = result*(mie::Vuint(exponetiate(eval_points[thread_idx][k],right_degree).get_str().c_str()));
            }
        }
        results[k] =  result;
    }
    //cout<<"thread_idx: "<<thread_idx<<endl;
    if(side == 0){
        for (int i = 0; i < num_eval_pt; ++i)
        {
            g_fL[thread_idx*num_eval_pt+i] =  results[i];
        }
    }
    else{
        for (int i = 0; i < num_eval_pt; ++i)
        {
            g_fR[thread_idx*num_eval_pt+i] =  results[i];
        }
    }

    return results;

}
Ec1* multiPointsEval_naive(Ec1* A, mpz_class* eval_points, int A_size, int deg_A, int num_eval_pt, Ec1 g1){
    Ec1* results = new Ec1[num_eval_pt];
    mpz_class idx;
    for(int k = 0;k<num_eval_pt;k++){
        Ec1 result = A[0]*(mie::Vuint(exponetiate(eval_points[k],zero).get_str().c_str()));
        for (int l = 1; l < A_size; l++)
        {
            idx = mpz_class(l);
            result += A[l]*(mie::Vuint(exponetiate(eval_points[k],idx).get_str().c_str()));
        }
        results[k] =  result;
    }
    return results;

}
Ec1* compute_gf(Ec1* A, int deg_A,int A_size ,mpz_class* eval_points, int num_eval_pt,Ec1 g1){
    int threshold_degree = 11;
    if(num_eval_pt<pow(2,threshold_degree)){//2^10
        Ec1* results = multiPointsEval_naive(A, eval_points, A_size, deg_A, num_eval_pt, g1);
        return results;
    }

    cout<<"poly size: "<<A_size<<endl;
    cout<<"deg: "<<deg_A<<endl;
    
    int A_size2 = deg_A*2+1;// I don't know why this is. But ok....
    Ec1* A2 = new Ec1[A_size2];
    for (int i = 0; i < A_size2; ++i)
    {
        if(i<A_size){
            A2[i] = A[i];
        }
        else{
            A2[i] = g1* ling;
        }
    }

    //eval FFT
    //cout<<"=============FFT result==============="<<endl;
    Ec1* eval_FFT = new Ec1[num_eval_pt];
    int num_eval_batch =1;
    if(deg_A!=0){
        int log_newsize = floor(log2(deg_A));
        num_eval_batch = pow(2,log_newsize);
    }
    
    int iteration = (num_eval_pt)/(num_eval_batch);
    cout<<"num_eval_pt: "<<num_eval_pt<<endl;
    cout<<"num iteration: "<<iteration<<endl;
    int *indices_test = new int [num_eval_batch];
    cout<<"num_eval_batch: "<<num_eval_batch<<endl;

    for (int j = 0; j < num_eval_batch; ++j)
    {
        indices_test[j] = j;
    }
    Ec1* eval_FFT_batch = new Ec1[num_eval_batch];
    int index = 0;
    for (int i = 0; i < iteration; ++i)
    {  
        mpz_class* eval_points_batch = new mpz_class[num_eval_batch];
        for (int j = 0; j < num_eval_batch; ++j)
        {
            eval_points_batch[j] = eval_points[i*num_eval_batch + j];
        }

        if(num_eval_batch<=pow(2,threshold_degree)){//use single point eval
           eval_FFT_batch = multiPointsEval_naive(A, eval_points_batch, A_size, deg_A, num_eval_batch, g1);
        }
        else{//use double point eval;
            multiPointsEval(A2, eval_points_batch, A_size2, deg_A, num_eval_batch, indices_test,eval_FFT_batch, g1);
        }

        // clock_t t = clock();
        // multiPointsEval(A2, eval_points_batch, A_size2, deg_A, num_eval_batch, indices_test,eval_FFT_batch, g1);
        // multieval_time += clock()-t;
        for (int j = 0; j < num_eval_batch; ++j)
        {
            eval_FFT[index] = eval_FFT_batch[j];
            index++;
        }
        delete [] eval_points_batch;
    }
    delete [] eval_FFT_batch;
    delete [] indices_test;
    return eval_FFT;
}
Ec1* prover_set_broadcast_values_and_return_gy(int newsize,Ec1 g1,Ec1* a,  int deg_A, mpz_class* eval_points, int num_eval){
    //BBM.coommit initial

    //step 1: prover compute the eval vectors


    T = new mpz_class*[num_eval]; //all eval vectors
    for (int i = 0; i < num_eval; ++i)
    {
        T[i] = new mpz_class[newsize];
        for (int j = 0; j < newsize; ++j)
        {
            if(j > deg_A){
                T[i][j] = zero;
            }else{
                T[i][j] = exponetiate(eval_points[i],j);
            }
        }
    }


    matrix_T = new mpz_class*[newsize];
    for (int i = 0; i < newsize; ++i)
    {
        matrix_T[i] = new mpz_class[num_eval];
        for (int j = 0; j < num_eval; ++j)
        {
            matrix_T[i][j] = T[j][i];
        }
    }

    //batch all the T
    r = new mpz_class[num_eval];//randomness for batching
    b = new mpz_class[newsize];
    for (int i = 0; i < num_eval; ++i)
    {
        int val = rand()%100+1;
        r[i] = mpz_class(val);
    }

    for (int i = 0; i < num_eval; ++i)
    {
        for (int j = 0; j < newsize; ++j)
        {
            //T[i][j] = T[i][j] * r[i];
            b[j] = (b[j] + T[i][j] * r[i]);
        }
    }

    for (int i = 0; i < newsize; ++i)
    {
        b[i] = b[i]%mod;
    }
    

    Ec1 d = innerEc1Mpz(a,b,newsize);
    com = bbm_enc(a,b,d,v,w,newsize);
    clock_t t_1 = clock();
    Ec1* evaluations = new Ec1[num_eval];
    auto t_eval = make_timer("Phase 1 evaluation");
    evaluations = compute_gf(a,deg_A,newsize,eval_points,num_eval,g1);  
    t_eval.stop();

    float eval1_time = clock() - t_1;

    // clock_t t_2 = clock();
    // Ec1* evaluations2 = new Ec1[num_eval];
    // for (int i = 0; i < num_eval; ++i)
    // {
    //     evaluations2[i] = innerEc1Mpz(a,T[i],newsize);
    // }
    // float eval2_time = clock() - t_2;
    // printCipherPoly(evaluations,5);
    // cout<<"================================"<<endl;
    // printCipherPoly(evaluations2,5);
    // evaluations = evaluations2;
    // exit(1);

    cout<<"eval 1 time: "<<(double)(eval1_time)/CLOCKS_PER_SEC<<"s\n";
    // cout<<"eval 2 time: "<<(double)(eval2_time)/CLOCKS_PER_SEC<<"s\n";
    // ===================debug========================
    // cout<<com<<endl;
    clock_t t_com_Ti = clock();
    compute_com_Ti(num_eval_pt, newsize);    

    //compute_com_Ti(num_eval_pt, newsize);    
    comTi_time = clock() - t_com_Ti;

    Ec1 d_2 = evaluations[0] * (mie::Vuint(r[0].get_str().c_str()));
    for (int i = 1; i < num_eval; ++i)
    {
        d_2 = d_2 + (evaluations[i] * (mie::Vuint(r[i].get_str().c_str())));
    }

    Fp12 com_2 = compute_com_phase1(a,b,d_2,v,w,r,newsize,num_eval);
    com = com_2;
   
    return evaluations;

}

void Prover(Ec2 g2, Ec1 g1,int d,int deg_A, Ec1* A ,Ec2* v, Ec1* w, mpz_class* eval_points, mpz_class* b, int num_eval,Fp12 com,int round){
    //printCipher2Poly(v,d);
    //cout<<com<<endl;

    if(d>=2){
    printf("=============Prover iteration : %d==================\n",round);

        //cout<<"degree_A: "<<deg_A<<endl<<endl;
        //cout<<"round: "<<round<<endl;
        int d_prime = d/2;
        //split the A, b into 2 halves
        clock_t time = clock();

        Ec1* A_L = new Ec1[d_prime];
        Ec1* A_R = new Ec1[d_prime];

        mpz_class* b_L = new mpz_class[d_prime];
        mpz_class* b_R = new mpz_class[d_prime];


        for (int i = 0; i < d_prime; i++)
        {
            A_L[i] = A[i];
            A_R[i] = A[i+d_prime];
            b_L[i] = b[i]%mod;
            b_R[i] = b[i+d_prime]%mod;
        }

        //split the pk v,w into 2 halves
        Ec2* v_L = new Ec2[d_prime];
        Ec2* v_R = new Ec2[d_prime];
        
        Ec1* w_L = new Ec1[d_prime];
        Ec1* w_R = new Ec1[d_prime];
        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];
            w_L[i] = w[i];
            w_R[i] = w[i+d_prime];
        }

        //compute cl, cr
        Ec1 cl,cr;

      
        cl = innerEc1Mpz(A_R,b_L,d_prime);
        cr = innerEc1Mpz(A_L,b_R,d_prime);


        // cl = innerEc1Mpz(A_L,b_R,d_prime);
        // cr = innerEc1Mpz(A_R,b_L,d_prime);

        Fp12 L, R;
        //compute L = bbm.encrytion

        L = bbm_enc(A_R, b_L,cl,v_L,w_R,d_prime);
        R = bbm_enc(A_L, b_R,cr,v_R,w_L,d_prime);


        // L = bbm_enc_A(A_R,v_L,d_prime);
        // R = bbm_enc_A(A_L,v_R,d_prime);
        // L = bbm_enc(A_L, b_R,cl,v_R,w_L,d_prime);
        // R = bbm_enc(A_R, b_L,cr,v_L,w_R,d_prime);

        //store the proofs
        proof_L[round] = L;
        proof_R[round] = R;


        //grab challenge from chals;
        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mpz_class x_inv_mpz = invert_mpz_class(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());


        //proofsize
        int proof_size_each = 2* sizeof(Fp12);
        proof_size += proof_size_each;


        //"fold" the secret vector A,b and v,w
        time = clock();
        Ec1* A_prime = new Ec1[d_prime];
        Ec2* v_prime = new Ec2[d_prime];
        mpz_class* b_prime = new mpz_class[d_prime];
        Ec1* w_prime = new Ec1[d_prime];

        mpz_class x_temp2 = x_temp*x_temp;
        mie::Vuint x_minu2 = invert(x_temp2);
        mie::Vuint x_2(x_temp2.get_str().c_str());

        for (int i = 0; i < d_prime; i++)
        {   
            A_prime[i] = A_L[i] + A_R[i]*x;
            v_prime[i] = v_L[i] + v_R[i]*x_inv;

            b_prime[i] = b_L[i] + b_R[i]*x_inv_mpz;
            w_prime[i] = w_L[i] + w_R[i]*x;

        }
        //store debug proof

        proof_A[round] = A_prime;
        proof_b[round] = b_prime;
        Fp12 com_prime = power(L , x) * com * power(R ,x_inv);
        folding_time += (double)(clock()-time)/CLOCKS_PER_SEC;




       
        //clean memory
        delete[] v_R;
        delete[] v_L;
        delete[] A_R;
        delete[] A_L;


        //recurse 
        round = round + 1;
        deg_A = d_prime-1;
        Prover(g2,g1,d_prime, deg_A,A_prime,v_prime,w_prime,eval_points,b_prime,num_eval,com_prime,round);
        return;
        
    }
    else{
    printf("=============Prover Final iteration==================\n");

        //create final proof
        final_proof_1 = A[0];
        final_proof_2 = b[0]%mod;
        
        proof_size += sizeof(final_proof_1)/3;
        proof_size += sizeof(final_proof_2);

        return;
    }
}

void Verifier(Ec2 g2, Ec1 g1, int d, Ec2* v, Ec1* w, Fp12 com, int num_eval, int round){
    //cout<<com<<endl;
    if(d>=2){
    printf("=============Verifier iteration : %d==================\n",round);

        int d_prime = d/2;
        //receive msgs from the prover
        Fp12 L = proof_L[round];
        Fp12 R = proof_R[round];

        clock_t time = clock();
        //compute v_prime, w_prime
        Ec2* v_prime = new Ec2[d_prime];
        Ec2* v_L = new Ec2[d_prime];
        Ec2* v_R = new Ec2[d_prime];


        Ec1* w_prime = new Ec1[d_prime];
        Ec1* w_L = new Ec1[d_prime];
        Ec1* w_R = new Ec1[d_prime];
        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];

            w_L[i] = w[i];
            w_R[i] = w[i+d_prime];
        }
        

        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mpz_class x_inv_mpz = invert_mpz_class(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());

        mpz_class x_temp2 = x_temp*x_temp;
        mie::Vuint x_minu2 = invert(x_temp2);
        mie::Vuint x_2(x_temp2.get_str().c_str());

        for (int i = 0; i < d_prime; i++)
        {
            v_prime[i] = v_L[i] + v_R[i]*x_inv;

            w_prime[i] = w_L[i] + w_R[i]*x;
        }


        folding_time_veri += (double)(clock()-time)/CLOCKS_PER_SEC;

        Fp12 com_prime = power(L , x) * com * power(R ,x_inv);

        //com_prime = com;
        

       
          
        
        delete[] v_L;
        delete[] v_R;

        //recursion
        round = round +1;
        Verifier(g2, g1,d_prime, v_prime,w_prime, com_prime,num_eval, round);
        return;
    }
    else{//check com = ped.commit(C,b,<C,b>)
        printf("=============Verifier Final iteration ==================\n");
       Ec1 * vec1 = new Ec1[1];
       mpz_class * vec2 = new mpz_class[1];
       vec1[0] = final_proof_1;
       vec2[0] = final_proof_2%mod;


       Ec1 final_innerpro = vec1[0]*(mie::Vuint(vec2[0].get_str().c_str()));
    
       Fp12 com2 = bbm_enc(vec1,vec2,final_innerpro,v,w,1);

       //cout<<com2<<endl;
       if (com == com2)
       {
         cout<<"Verification Result: PASS"<<endl;
       }
       else{
         cout<<"Verifier: FAIL"<<endl;

       }
      
    }
}

//Prover_phase2(g2,g1,newsize,deg_A,L_phase3,w,w_prime,R_phase3,num_eval_pt,0);

void Prover_phase2(Ec2 g2, Ec1 g1,int d,int deg_A, mpz_class* A ,Ec1* v, Ec1* w, mpz_class* b, int num_eval,int round){
   
    printf("=============Prover iteration : %d==================\n",round);
    if(d>=2){
        //cout<<"degree_A: "<<deg_A<<endl<<endl;
        //cout<<"round: "<<round<<endl;
        int d_prime = d/2;
        //split the A, b into 2 halves
        clock_t time = clock();

        mpz_class* A_L = new mpz_class[d_prime];
        mpz_class* A_R = new mpz_class[d_prime];

        mpz_class* b_L = new mpz_class[d_prime];
        mpz_class* b_R = new mpz_class[d_prime];

        for (int i = 0; i < d_prime; i++)
        {
            A_L[i] = A[i]%mod;
            A_R[i] = A[i+d_prime]%mod;
            b_L[i] = b[i]%mod;
            b_R[i] = b[i+d_prime]%mod;
        }
        //split the pk v,w into 2 halves
        Ec1* v_L = new Ec1[d_prime];
        Ec1* v_R = new Ec1[d_prime];
        
        Ec1* w_L = new Ec1[d_prime];
        Ec1* w_R = new Ec1[d_prime];
        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];
            w_L[i] = w[i];
            w_R[i] = w[i+d_prime];
        }


        //compute cl, cr
        mpz_class cl,cr;
        // cl = innerEc1Mpz(A_R,b_L,d_prime);
        // cr = innerEc1Mpz(A_L,b_R,d_prime);


        cl = innerMpzMpz(A_L,b_R,d_prime);
        cr = innerMpzMpz(A_R,b_L,d_prime);

        Ec1 L, R;
        //compute L = bbm.encrytion

        // L = bbm_enc(A_R, b_L,cl,v_L,w_R,d_prime);
        // R = bbm_enc(A_L, b_R,cr,v_R,w_L,d_prime);

        L = ped_enc(A_L, b_R,cl,v_R,w_L,d_prime);
        R = ped_enc(A_R, b_L,cr,v_L,w_R,d_prime);

        //store the proofs
        proof_L_phase2[round] = L;
        proof_R_phase2[round] = R;

        //proofsize
        int proof_size_each = 2* sizeof(Ec1)/3;
        proof_size += proof_size_each;

        //grab challenge from chals;
        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mpz_class x_inv_mpz = invert_mpz_class(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());





        //"fold" the secret vector A,b and v,w
        time = clock();
        mpz_class* A_prime = new mpz_class[d_prime];
        Ec1* v_prime = new Ec1[d_prime];
        mpz_class* b_prime = new mpz_class[d_prime];
        Ec1* w_prime = new Ec1[d_prime];

        mpz_class x_temp2 = x_temp*x_temp;
        mie::Vuint x_minu2 = invert(x_temp2);
        mie::Vuint x_2(x_temp2.get_str().c_str());

        for (int i = 0; i < d_prime; i++)
        {   
            A_prime[i] = A_L[i]*x_temp + A_R[i]*x_inv_mpz;
            v_prime[i] = v_L[i]*x_inv + v_R[i]*x;

            b_prime[i] = b_L[i] * x_inv_mpz + b_R[i]*x_temp;
            w_prime[i] = w_L[i] * x + w_R[i]*x_inv;

        }


        folding_time += (double)(clock()-time)/CLOCKS_PER_SEC;
        //printCipher2Poly(v,d);
        //printCipherPoly(w_prime,d_prime);
        //clean memory
        delete[] v_R;
        delete[] v_L;
        delete[] A_R;
        delete[] A_L;


        //recurse 
        round = round + 1;
        deg_A = d_prime-1;
        Prover_phase2(g2,g1,d_prime, deg_A,A_prime,v_prime,w_prime,b_prime,num_eval,round);
        return;
        
    }
    else{
        //create final proof
        final_proof_1_phase2 = A[0];
        final_proof_2_phase2 = b[0];
        proof_size+= sizeof(final_proof_1_phase2);
        proof_size+= sizeof(final_proof_2_phase2);

        return;
    }
}

void Verifier_phase2(Ec2 g2, Ec1 g1, int d, Ec1* v, Ec1* w, Ec1 com, int num_eval, int round){
    if(d>=2){
    printf("=============Verifier iteration : %d==================\n",round);

        int d_prime = d/2;
        //receive msgs from the prover
        Ec1 L = proof_L_phase2[round];
        Ec1 R = proof_R_phase2[round];

        clock_t time = clock();
        //compute v_prime, w_prime
        Ec1* v_prime = new Ec1[d_prime];
        Ec1* v_L = new Ec1[d_prime];
        Ec1* v_R = new Ec1[d_prime];


        Ec1* w_prime = new Ec1[d_prime];
        Ec1* w_L = new Ec1[d_prime];
        Ec1* w_R = new Ec1[d_prime];
        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];

            w_L[i] = w[i];
            w_R[i] = w[i+d_prime];
        }
        

        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mpz_class x_inv_mpz = invert_mpz_class(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());

        mpz_class x_temp2 = x_temp*x_temp;
        mie::Vuint x_minu2 = invert(x_temp2);
        mie::Vuint x_2(x_temp2.get_str().c_str());

        for (int i = 0; i < d_prime; i++)
        {
            v_prime[i] = v_L[i]*x_inv + v_R[i]*x;

            w_prime[i] = w_L[i]*x + w_R[i]*x_inv;
        }

        //printCipher2Poly(v,d);
        //printCipherPoly(w_prime,d_prime);
        folding_time_veri += (double)(clock()-time)/CLOCKS_PER_SEC;

        Ec1 com_prime;
        com_prime = L * x_2 + com + R * x_minu2;
        //com_prime = L * x_minu2 + com + R * x_2;

               
        delete[] v_L;
        delete[] v_R;

        //recursion
        round = round +1;
        Verifier_phase2(g2, g1,d_prime, v_prime,w_prime, com_prime,num_eval, round);
        return;
    }
    else{//check com = ped.commit(C,b,<C,b>)
        printf("=============Verifier Final iteration ==================\n");
       mpz_class * vec1 = new mpz_class[1];
       mpz_class * vec2 = new mpz_class[1];
       vec1[0] = final_proof_1_phase2%mod;
       vec2[0] = final_proof_2_phase2%mod;


       mpz_class final_innerpro = (vec1[0]*vec2[0])%mod;
       Ec1 com2 = ped_enc(vec1,vec2,final_innerpro,v,w,1);

       if (com == com2)
       {
         cout<<"Linear constriant Test: PASS"<<endl;
       }
       else{
         cout<<"Linear constriant Test: FAIL"<<endl;
       }
    }
}
void Prover_Phase3(mpz_class &d,Ec1 com_d,mpz_class x,int m,int Li_size,int round){
    if(m>=1){
        //compute_l_r
        mpz_class l,r;
        l = 0;
        r = 0;
        for (int i = 0; i < m; ++i)
        {
            l = (l + innerMpzMpz(Li[i],Ri[i+m],Li_size))%mod+mod;
            r = (r + innerMpzMpz(Li[i+m],Ri[i],Li_size))%mod+mod;
        }
        Ec1 com_l = cku * (mie::Vuint(l.get_str().c_str()));
        Ec1 com_r = cku * (mie::Vuint(r.get_str().c_str()));

        //store proof
        proof_coml[round] = com_l;
        proof_comr[round] = com_r;

        //proofsize
        int proof_size_each = 2* sizeof(Ec1)/3;
        proof_size += proof_size_each;


        //fold Li Ri
        for (int i = 0; i < m; ++i)
        {

            Li[i] = mpzAddmpz(Li[i], mpzTimesMpzVector(invert_mpz_class(x), Li[i+m],Li_size),Li_size);
            Ri[i] = mpzAddmpz(Ri[i], mpzTimesMpzVector(x, Ri[i+m],Li_size),Li_size);
        }
           



        //update d, m
        d = (d + x*l + invert_mpz_class(x)*r)%mod;



        // mie::Vuint x_(x.get_str().c_str());
        // mie::Vuint x_inv = invert(x);
        // Ec1 com_d2 = com_d +  com_l*x_ + com_r*x_inv;
        // cout<<com_d2<<endl;

        mie::Vuint d_(d.get_str().c_str());
        //cout<<cku*d_<<endl;


        if(m==1){
            return;
        }
        m = m/2;
        round = round+1;
        Prover_Phase3(d,com_d,x,m,Li_size,round);
    }
    else{

        return;
    }
    
}

void Verifier_Phase3(mpz_class d,Ec1 &com_d,mpz_class x,int m,int Li_size,int round){
    cout<<m<<endl;
    if(m>=1){

        //take proof out
        Ec1 com_l = proof_coml[round];
        Ec1 com_r = proof_comr[round];

        mie::Vuint x_(x.get_str().c_str());
        mie::Vuint x_inv = invert(x);




        //fold com_Li, com_Ri;
        for (int i = 0; i < m; ++i)
        {

            com_Li[i] = com_Li[i] + com_Li[i+m]*x_inv;
            com_Ri[i] = com_Ri[i] + com_Ri[i+m]*x_;
        }
        //w = ec1AddEc1(w, mpzTimesEc1Vector(invert_mpz_class(x),w,Li_size),Li_size);


        //fold cross term 
        com_d = com_d +  com_l*x_ + com_r*x_inv;
        //cout<<com_d<<endl;

        if(m==1){
            return;
        }
        m = m/2;
        round = round+1;
        Verifier_Phase3(d,com_d,x,m,Li_size,round);
    }
    else{
        return;
    }


}
void run_server(std::chrono::seconds accept_time, const std::vector<mpz_class>& input, const std::string bind_addr, uint16_t bind_port) {

    /*  Initialize a executor */
    asio::io_context executor;

    /*  Construct a TCP server */
    tcp_server sock(executor);

    /*  Optional: set worker threads */
    sock.threads(2);

    DEBUG << "Binding to 127.0.0.1:5555";

    /*  Bind the server to localhost:5555 */
    auto accepted = sock.bind(bind_addr, bind_port, accept_time);
    //auto accepted = sock.bind("127.0.0.1", 5555, accept_time);

    std::cout << "Accepted " << accepted << " connections" << std::endl;
    INFO << "Accepted " << accepted << " connections";

    /*  Notify clients to start */
    sock.broadcast("start"s);

    auto t = make_timer("P1 protocol");
    auto t1 = make_timer("Key Generation");
    /*  Prepare random seed */
    mpz_class mpr = mpz_class((Param::r).toString());
    DEBUG <<"mpr: "<<mpr<<endl;
    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));




    /*  Deserialize received messages and send back */
    {
        auto [msgs, kickedout] = sock.collect(100s);

        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }

    // Key Generation
    std::cout << "Key Generation" << std::endl;
    Ec1 Pusum,Pvsum,Phsum;
    bool flag = true;
        for (auto& msg : msgs) {
            Ec1 u, v, h;
            sock.archive().unpack(msg, u, v, h);

            std::stringstream os;
            os << "Get g^u = " << u << std::endl
               << "Get g^v = " << v << std::endl
               << "Get g^h = " << h << std::endl;
            DEBUG << os.str();
        //std::cout << os.str();

        if(flag) {
          Pusum = u;
          Pvsum = v;
          Phsum = h;
              flag = false;
        }
        else {
              Pusum += u;
              Pvsum += v;
              Phsum += h;
        }
        }

        std::stringstream os;
        os << "Sending Public key u = " << Pusum << std::endl
           << "Sending Public key v = " << Pvsum << std::endl
           << "Sending Public key h = " << Phsum << std::endl;
        DEBUG << os.str();
    //std::cout << os.str();
        sock.broadcast(Pusum, Pvsum, Phsum);
        total_communication += sizeof(Ec1);
        pubkey key{ Pusum, Pvsum, Phsum };
    t1.stop();
        //auto t2 = make_timer("Aggregate all polynomials");

        // Combined polynomial
        std::optional<std::vector<ciphertext>> enc_poly_sum = std::nullopt;
        std::tie(msgs, kickedout) = sock.collect(100s);
        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
        auto t_agg = make_timer("Aggregate all polynomials");
        for (int j = 0; j < 100; ++j)
        {
            for (auto& msg : msgs) {
                std::vector<ciphertext> enc_poly;
                sock.archive().unpack(msg, enc_poly);
                DEBUG << "Receiving ..." << std::endl;

                if (enc_poly_sum) {
                    
                         for (auto i = 0; i < enc_poly.size(); i++) {
                            enc_poly_sum.value()[i] = enc_poly_sum.value()[i] + enc_poly[i];
                        }
               
                }
                
                else {
                    enc_poly_sum = std::move(enc_poly);
                }
            }
        }
        
        t_agg.stop();

        //t2.stop();
        auto t3 = make_timer("Interactive Evaluate at lambda");

    // Evaluate on a random point lambda
    std::cout << "Lambda test" << std::endl;

    mpz_class lambda = rr.get_z_bits(400) % mpr;
    auto lambda_eval = homorphic_polyeval(key,g1,lambda,enc_poly_sum.value());
    //cout<<"lamda eval: "<<lambda_eval.get(2)<<endl;

//===========================================================================
//================lamda eval using phase1 code =============================
//===========================================================================


    int poly_size = enc_poly_sum->size();
    keygen_public(poly_size,g1, g2);
    mpz_class* t_lambda = new mpz_class[poly_size]; //all eval vectors
    mpz_class* lamdas = new mpz_class[1];
    for (int i = 0; i < poly_size; ++i)
    {
        t_lambda[i] = exponetiate(lambda,i);
    }
    Ec1* poly = new Ec1[poly_size];
    for (int i = 0; i < poly_size; ++i)
    {
        poly[i] = enc_poly_sum.value()[i][2];
    }
    Ec1 d = innerEc1Mpz(poly,t_lambda,poly_size);
    Fp12 com_lamda = bbm_enc(poly,t_lambda,d,v,w,poly_size);

    int rounds = log2(poly_size);//rounds = how many chals or proofs
    //generate challenges before hand
    chals = generate_phase1_random_challenge(rounds);  
    init_proof(rounds);//assign the correct size to each proof vector

    Prover(g2,g1,poly_size,0,poly,v,w,lamdas,t_lambda,1,com_lamda,0);
    std::vector<Ec2> v_send(poly_size);
    for (int i = 0; i < poly_size; ++i)
    {
        v_send[i] = v[i];
    }
    std::vector<Ec1> w_send(poly_size);
    for (int i = 0; i < poly_size; ++i)
    {
        w_send[i] = w[i];
    }

    sock.broadcast(poly_size, v_send, w_send,h1,h2,cku, rounds, com_lamda);
    sock.broadcast(proof_L, proof_R);
    sock.broadcast(final_proof_1, final_proof_2.get_str());

    vector<string>chals_lambda(rounds);
    for (int i = 0; i < rounds; ++i)
    {
        chals_lambda[i] = chals[i].get_str();

    }
    sock.broadcast(chals_lambda);
//===========================================================================
//===========================================================================

    // sock.broadcast(lambda.get_str());

    // // Combined evaluation at lambda 
    //     std::optional<ciphertext> lambda_eval_sum = std::nullopt;
    //     std::tie(msgs, kickedout) = sock.collect(100s);
    //     if (kickedout) {
    //         WARNING << kickedout << " clients have been kicked out";
    //     }
        
    //     for (auto& msg : msgs) {
    //         ciphertext enc_lambda_eval;
    //         sock.archive().unpack(msg, enc_lambda_eval);
    //         DEBUG << "Receiving ..." << std::endl;

    //         if (lambda_eval_sum) {
    //                 lambda_eval_sum.value() = lambda_eval_sum.value() + enc_lambda_eval;
    //         }
    //         else {
    //             lambda_eval_sum = std::move(enc_lambda_eval);
    //         }
    // }

    /*
    std::cout<<lambda_eval_sum.value()[0] << std::endl;
    std::cout<<            lambda_eval[0] << std::endl;
    std::cout<<lambda_eval_sum.value()[1] << std::endl;
    std::cout<<            lambda_eval[1] << std::endl;
    std::cout<<lambda_eval_sum.value()[2] << std::endl;
    std::cout<<            lambda_eval[2] << std::endl;*/
    //sock.broadcast(lambda_eval_sum.value());
    t3.stop();

    //=================WENXUAN malicious=================================//

    std::cout << "Proof Generation" << std::endl;
    
    


    int deg_A = enc_poly_sum->size()-1;
    DEBUG<<"deg_A"<<deg_A<<endl;
    int coeff_len = deg_A+1;
    int newsize = coeff_len;
    int oldsize = newsize;
    std::vector<Ec1> A_(coeff_len);
    for (int i = 0; i < enc_poly_sum->size(); ++i)
    {
        A_[i] = enc_poly_sum.value()[i][2];
    }
    if(isPowerOfTwo(coeff_len)==0){
        int log_newsize = ceil(log2(coeff_len));
        newsize = pow(2,log_newsize);
        A_.resize(newsize);
    }

    Ec1* A = new Ec1[newsize];// this is the true f(x)
    for (int i = 0; i < newsize; ++i)
    {
        if(i<oldsize){
            A[i] = A_[i];
        }else{
            A[i] = g1*mie::Vuint(0);
        }
    }
    int num_eval = input.size();
    DEBUG <<"input size: "<<num_eval<<endl;
    eval_points = new mpz_class[num_eval];
    for (int i = 0; i < num_eval; ++i)
     {
         eval_points[i] = input[i];
         //cout<<"a[i]"<<a[i]<<endl;
     } 
    keygen_public(newsize,g1, g2);
    rounds = log2(newsize);//rounds = how many chals or proofs
    //generate challenges before hand
    chals = generate_phase1_random_challenge(rounds);  
    init_proof(rounds);//assign the correct size to each proof vector
    compute_look_up_table2(w, newsize);

    cout<<"newsize: "<<newsize<<endl;
    cout<<"num_eval: "<<num_eval<<endl;
    cout<<"deg_A: "<<deg_A<<endl;

    num_eval_pt = num_eval;
    cout<<"==============================Begin PHASE 1================================"<<endl;

    auto t4 = make_timer("Interactive Proof Generation");
    auto t4a = make_timer("Proof Generation Phase 1");
    clock_t t_prover = clock();
    Ec1* gy = prover_set_broadcast_values_and_return_gy(newsize,g1,A,deg_A,eval_points,num_eval_pt);//set com return gy

    commitment_size += num_eval*sizeof(Ec1)/3;
    commitment_size += num_eval*sizeof(com);

    Prover(g2, g1,newsize,deg_A, A,v,w,eval_points,b,num_eval_pt,com,0);
    cout<<"prover time: "<<(double)(clock()-t_prover)/CLOCKS_PER_SEC<<"s\n";

    final_proof_2 = final_proof_2%mod;
    std::vector<Ec1> A_send(newsize);
    for (int i = 0; i < newsize; ++i)
    {
        A_send[i] = A[i];
    }
    t4a.stop();

    sock.broadcast(newsize, num_eval,
                   deg_A, A_send, rounds, com);
    // vector<Ec1> com_Ti_send;
    // com_Ti_send.resize(num_eval_pt);
    // for (int i = 0; i < num_eval_pt; ++i)
    // {
    //     com_Ti_send[i] = com_Ti[i];
    // }

    // sock.broadcast(com_Ti_send);   
    sock.broadcast(proof_L, proof_R);

    final_proof_2 = final_proof_2%mod;

    sock.broadcast(final_proof_1, final_proof_2.get_str());
    
    //send chals
    // vector<mpz_class> chals;
    vector<string>chals_str(rounds);
    for (int i = 0; i < rounds; ++i)
    {
        chals_str[i] = chals[i].get_str();
        //cout<<chals[i]<<endl;
    total_communication+= rounds * sizeof(chals_str[i]);

    }
    auto t4d = make_timer("Proof Generation Phase 4");
    sock.broadcast(chals_str);
    t4d.stop();

    cout<<"==============================END PHASE 1================================"<<endl;
    cout<<"==============================Begin PHASE 2================================"<<endl;
    clock_t t_phase2 = clock();

    mpz_class chal_x;
    mpz_urandomm(chal_x.get_mpz_t(),r_state,p.get_mpz_t()); 
    init_selector_matrices(newsize);
    init_phase2_randomness(newsize,num_eval_pt);
 
    clock_t t_XYZ = clock();    
    compute_X_Y_Z_F(num_eval_pt, newsize);
    float computematrixtime = clock() - t_XYZ;

    //prover commits all Xi;
    clock_t t_com_Xi = clock();
    compute_com_Xi_com_Yi_com_Zi_com_Fi(newsize,num_eval_pt);  
    commitment_size+= 4* num_eval_pt * sizeof(Ec1)/3;
    //send com_Xi com_Yi com_Zi com_Fi
    // vector<Ec1> com_Xi_send;
    // vector<Ec1> com_Yi_send;
    // vector<Ec1> com_Zi_send;
    // vector<Ec1> com_Fi_send;

    // com_Xi_send.resize(num_eval_pt);
    // com_Yi_send.resize(num_eval_pt);
    // com_Zi_send.resize(num_eval_pt);
    // com_Fi_send.resize(num_eval_pt);

    // for (int i = 0; i < num_eval_pt; ++i)
    // {
    //     com_Xi_send[i] = com_Xi[i];
    //     com_Yi_send[i] = com_Yi[i];
    //     com_Zi_send[i] = com_Zi[i];
    //     com_Fi_send[i] = com_Fi[i];

    // }
    // vector<string> s_send;
    // vector<string> u_send;

    // for (int i = 0; i < newsize-1; ++i)
    // {
    //     s_send[i] = s[i].get_str();
    // }
    // for (int i = 0; i < num_eval; ++i)
    // {
    //     u_send[i] = u[i].get_str();
    // }
    //sock.broadcast(com_Xi_send,com_Yi_send,com_Zi_send,com_Fi_send,chal_x.get_str());
    float comXi_time = clock() - t_com_Xi;

    compute_Sa_Xu_Tu_Sb_Yu_Sc_Zu_Sd_Fu( newsize,num_eval_pt);  
    //send sa,sb,sc,sd. this is public.
    // vector<string> Sa_send;
    // vector<string> Sb_send;
    // vector<string> Sc_send;
    // vector<string> Sd_send;
    // Sa_send.resize(newsize);
    // Sb_send.resize(newsize);
    // Sc_send.resize(newsize);
    // Sd_send.resize(newsize);
    // for (int i = 0; i < newsize; ++i)
    // {
    //     Sa_send[i] = Sa[i].get_str();
    //     Sb_send[i] = Sb[i].get_str();
    //     Sc_send[i] = Sc[i].get_str();
    //     Sd_send[i] = Sd[i].get_str();

    // }
    // sock.broadcast(Sa_send,Sb_send,Sc_send,Sd_send);
    //send commitments
    compute_l_r_commitments(cku,newsize,num_eval_pt,chal_x);
    //sock.broadcast(com_l_X,com_r_X,com_l_Y,com_r_Y,com_l_Z,com_r_Z,com_l_F,com_r_F);

    mpz_class* L_X = new mpz_class[newsize];
    mpz_class* R_X = new mpz_class[newsize];
    mpz_class* L_Y = new mpz_class[newsize];
    mpz_class* R_Y = new mpz_class[newsize];
    mpz_class* L_Z = new mpz_class[newsize];
    mpz_class* R_Z = new mpz_class[newsize];
    mpz_class* L_F = new mpz_class[newsize];
    mpz_class* R_F = new mpz_class[newsize];
    compute_L_and_R(L_X,R_X,chal_x,newsize,0);
    compute_L_and_R(L_Y,R_Y,chal_x,newsize,1);
    compute_L_and_R(L_Z,R_Z,chal_x,newsize,2);
    compute_L_and_R(L_F,R_F,chal_x,newsize,3);

    for (int i = 0; i < newsize; ++i)
    {
        L_X[i] = L_X[i]+mod;
        R_X[i] = R_X[i]+mod;
        L_Y[i] = L_Y[i]+mod;
        R_Y[i] = R_Y[i]+mod;
        L_Z[i] = L_Z[i]+mod;
        R_Z[i] = R_Z[i]+mod;
        L_F[i] = L_F[i]+mod;
        R_F[i] = R_F[i]+mod;
    }



    auto t_verfier_phase2_additional = make_timer("verfier pahse 2 additional");
    //Ec1 com_u = cku*(mie::Vuint(LR.get_str().c_str()));
    mie::Vuint x(chal_x.get_str().c_str());
    mie::Vuint x_inv = invert(chal_x);
    mpz_class x_inv_mpz = invert_mpz_class(chal_x);
    total_communication += sizeof(mpz_class);
    com_phase2_X = verifier_compute_com_phase2_X(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    com_phase2_Y = verifier_compute_com_phase2_Y(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    com_phase2_Z = verifier_compute_com_phase2_Z(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    com_phase2_F = verifier_compute_com_phase2_F(newsize,num_eval_pt,x,x_inv,x_inv_mpz);

    t_verfier_phase2_additional.stop();

    float phase2_verifier_time;
    clock_t t_innerprod = clock();
    std::vector<Ec1> v2_send(newsize);
    for (int i = 0; i < newsize; ++i)
    {
        v2_send[i] = v_2[i];
    }
    
    sock.broadcast(v2_send,com_phase2_X,com_phase2_Y,com_phase2_Z,com_phase2_F);


    Prover_phase2(g2,g1,newsize,deg_A,L_X,v_2,w,R_X,num_eval_pt,0);
    sock.broadcast(proof_L_phase2,proof_R_phase2,final_proof_1_phase2.get_str(),final_proof_2_phase2.get_str());
cout<<"=====================X sent ============================="<<endl;
    Prover_phase2(g2,g1,newsize,deg_A,L_Y,v_2,w,R_Y,num_eval_pt,0);
    sock.broadcast(proof_L_phase2,proof_R_phase2,final_proof_1_phase2.get_str(),final_proof_2_phase2.get_str());
cout<<"=====================Y sent ============================="<<endl;
    Prover_phase2(g2,g1,newsize,deg_A,L_Z,v_2,w,R_Z,num_eval_pt,0);
    sock.broadcast(proof_L_phase2,proof_R_phase2,final_proof_1_phase2.get_str(),final_proof_2_phase2.get_str());
cout<<"=====================Z sent ============================="<<endl;
    Prover_phase2(g2,g1,newsize,deg_A,L_F,v_2,w,R_F,num_eval_pt,0);
    sock.broadcast(proof_L_phase2,proof_R_phase2,final_proof_1_phase2.get_str(),final_proof_2_phase2.get_str());
cout<<"=====================F sent ============================="<<endl;

    float innerProducttime = clock() - t_innerprod;

    float phase2_time = clock() - t_phase2;
    cout<<"compute X,Y,Z time: "<<(double)(computematrixtime)/CLOCKS_PER_SEC<<"s\n";
    cout<<"Com Xi time: "<<(double)(comXi_time)/CLOCKS_PER_SEC<<"s\n";
    cout<<"Com innerProduct time: "<<(double)(innerProducttime)/CLOCKS_PER_SEC<<"s\n";
    cout<<"Verifier time: "<<(double)(phase2_verifier_time)/CLOCKS_PER_SEC<<"s\n";
    cout<<"Phase2 time: "<<(double)(phase2_time)/CLOCKS_PER_SEC<<"s\n";

    cout<<"==============================End PHASE 2================================"<<endl;
    cout<<"==============================Begin PHASE 3 (quadratic constraint)================================"<<endl;

    clock_t t_phase3 = clock();

    //step 1 verifier sends a random vector R, and randomness w
    init_phase3_randomness(newsize);

    Ec1* w_prime = new Ec1[newsize];
    for (int i = 0; i < newsize; ++i)
    {
        w_prime[i] = w[i] * invert(R[i]);
    }
    cout<<"check: "<<num_eval_pt<<endl;
    compute_Li_Ri_com_Li_com_Ri(w_prime,newsize, num_eval_pt);
    //send com_Li, com_Ri to verifier
    vector<Ec1> com_Li_vec;
    vector<Ec1> com_Ri_vec;
    com_Li_vec.resize(2*newsize);
    com_Ri_vec.resize(2*newsize);
    for (int i = 0; i < 2*newsize; ++i)
    {
        com_Li_vec[i] = com_Li[i];
        com_Ri_vec[i] = com_Ri[i];
    }
    vector<Ec1> w_prime_send;
    w_prime_send.resize(newsize);
    for (int i = 0; i < newsize; ++i)
    {
        w_prime_send[i] = w_prime[i];
    }
    sock.broadcast(com_Li_vec,com_Ri_vec,w_prime_send,cku, x_phase3.get_str());



    int Li_size =  newsize-1;
    com_d = cku* ling;
    Prover_Phase3(zero,com_d,x_phase3,newsize,Li_size,0);
    cout<<"prover phase 3 quadratic DONE"<<endl;
  
    //send phase3 proof
    sock.broadcast(proof_coml,proof_comr);


    mpz_class* L_phase3 = Li[0];
    mpz_class* R_phase3 = Ri[0];

   
    //compute final innerproduct
    Ec1 com_phase3 = verifier_compute_com_phase3(w_prime, L_phase3,R_phase3,num_eval_pt);
    cout<<com_phase3<<endl;

    Prover_phase2(g2,g1,newsize,deg_A,L_phase3,w,w_prime,R_phase3,num_eval_pt,0);
    sock.broadcast(proof_L_phase2,proof_R_phase2,final_proof_1_phase2.get_str(),final_proof_2_phase2.get_str());

    float phase3_time = clock() - t_phase3;
    cout<<"Phase3 time: "<<(double)(phase3_time)/CLOCKS_PER_SEC<<"s\n";
    cout <<"proof_size: "<<proof_size<<endl;

    t4.stop();

    cout<<"==============================END PHASE 3================================"<<endl;
    //=======================================================================//
    //=================END WENXUAN malicious=================================//
    //=======================================================================//

    auto t5 = make_timer("Multipoint evaluation");
    //wenxuan debug

    //clock_t test_muthu = clock();
    auto enc_eval = homorphic_polyeval(key,g1,input,enc_poly_sum.value());
    

    t5.stop();
    std::cout << "Decryption" << std::endl;
    auto t6 = make_timer("Interactive Decryption");

    auto reenc_eval = rerandomize(key,g1,enc_eval); 
    total_communication +=  reenc_eval.size()* sizeof(Ec1);   

    sock.broadcast(reenc_eval);
    // Rerandomization 
        std::optional<std::vector<ciphertext>> reenc_sum = std::nullopt;
        std::tie(msgs, kickedout) = sock.collect(100s);
        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
        
        for (auto& msg : msgs) {
        std::vector<ciphertext> reenc;
            sock.archive().unpack(msg, reenc);
        DEBUG << "Receiving ..." << std::endl;

            if (reenc_sum) {
                for (auto i = 0; i < reenc.size(); i++) {
                    reenc_sum.value()[i] = reenc_sum.value()[i] + reenc[i];
                }
            }
            else {
                reenc_sum = std::move(reenc);
            }
    }
    sock.broadcast(reenc_sum.value());
    total_communication +=  reenc_sum->size()* sizeof(Ec1);   

    std::optional<std::vector<Ec1>> dprime_sum = std::nullopt;
        std::tie(msgs, kickedout) = sock.collect(100s);
        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
        // Distributed decryption receiving dprime 
        for (auto& msg : msgs) {
        std::vector<Ec1> dprime;
            sock.archive().unpack(msg, dprime);
        DEBUG << "Receiving dprime" << std::endl;

            if (dprime_sum) {
                for (auto i = 0; i < dprime.size(); i++) {
                    dprime_sum.value()[i] = dprime_sum.value()[i] + dprime[i];
                }
            }
            else {
                dprime_sum = std::move(dprime);
            }
    }
    sock.broadcast(dprime_sum.value());
    total_communication +=  reenc_sum->size()* sizeof(Ec1)/3;   


        std::optional<std::vector<Ec1>> eprime_sum = std::nullopt;
        std::tie(msgs, kickedout) = sock.collect(100s);
        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
        // Distributed decryption receiving eprime 
        for (auto& msg : msgs) {
        std::vector<Ec1> eprime;
            sock.archive().unpack(msg, eprime);
        DEBUG << "Receiving eprime" << std::endl;

            if (eprime_sum) {
                for (auto i = 0; i < eprime.size(); i++) {
                    eprime_sum.value()[i] = eprime_sum.value()[i] + eprime[i];
                }
            }
            else {
                eprime_sum = std::move(eprime);
            }
    }
        for (auto i = 0;i < eprime_sum.value().size(); i++)
        {
          if(eprime_sum.value()[i].isZero()) std::cout << input[i] << ", ";
        }
        std::cout << std::endl;

        sock.broadcast(true);

        t6.stop();
        t.stop();
        show_timer();
        }

    INFO << "Done!";
    total_communication = total_communication + proof_size + commitment_size;
        cout<<"Communication Size: "<<total_communication<<endl;

}

// ------------------------------------------------------------

void run_client(const std::vector<mpz_class>& input, const std::string bind_addr, uint16_t bind_port) {

    mpz_class mpr = mpz_class((Param::r).toString());
    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));

    /*  Additional thread hint, optimized for single thread */
    asio::io_context executor(1);

    tcp_client sock(executor);

    /*  Connect to remote server  */
    std::string port(std::to_string(bind_port));
    std::cout << "Connecting..." << std::endl;
    sock.connect(bind_addr, port);
    //sock.connect("127.0.0.1", "5555");
    std::cout << "Connected" << std::endl;

    /*  Wait until server is ready... */
    auto [str] = sock.receive<std::string>(100s);

    auto t = make_timer("Client protocol");
    auto t1 = make_timer("Key Generation");
    mpz_class mpzx = rr.get_z_bits(400);
    mpzx = mpzx % mpr;
    mpz_class mpzy = rr.get_z_bits(400);
    mpzy = mpzy % mpr;
    mpz_class mpzz = rr.get_z_bits(400);
    mpzz = mpzz % mpr;
    //mpz_class mpzx(1),mpzy(2),mpzz(3);

    Ec1 Pu = g1 * mpzx;
    Ec1 Pv = g1 * mpzy;
    Ec1 Ph = g1 * mpzz;

    DEBUG << "Server says: " << str;

    // ------------------------------------------------------------
    // Key Generation
    // ------------------------------------------------------------

    std::stringstream os;
    os << "Sending g^u = " << Pu << std::endl
       << "Sending g^v = " << Pv << std::endl
       << "Sending g^h = " << Ph << std::endl;
    DEBUG << os.str();

    os.clear();
    os.str("");

    // ============================================================
    /*  Sending a vector like this  */
    // sock.send(std::vector<Ec1>{ Pu, Pv, Ph });
    // auto [vec] = sock.receive<std::vector<Ec1>>(100s);
    // ============================================================

    sock.send(Pu, Pv, Ph);
    total_communication += sizeof(Ec1);
    auto [sum_Pu, sum_Pv, sum_Ph] = sock.receive<Ec1, Ec1, Ec1>(100s);

    os << "Public key u = " << sum_Pu << std::endl
       << "Public key v = " << sum_Pv << std::endl
       << "Public key h = " << sum_Ph << std::endl;
    DEBUG << os.str();
    //std::cout << os.str();

    pubkey key{ sum_Pu, sum_Pv, sum_Ph };
    t1.stop();

    // ------------------------------------------------------------
    // Sending encrypted polynomials
    // ------------------------------------------------------------
    
    std::cout << "Sending encrypted polynomials." << std::endl;

    auto t2 = make_timer("Encryption");
    std::vector<ciphertext> enc_poly;

    // mpz_class* input_ptr = new mpz_class[input.size()];
    for(auto i=0; i < input.size(); i++) {
       // input_ptr[i] = input[i];
       DEBUG << "root [" << i << "] = " << input[i];

    }
    mpz_class* coeffs;
    coeffs = getCoeff(input.data(),input.size());



    for (auto i=0; i<input.size()+1; i++) {
        //enc_poly.emplace_back(encrypt(key, g1, coeffs[i]));
        DEBUG <<sock.identity() << ": c_" <<i << " = " << coeffs[i] % mpr << std::endl;
        enc_poly.emplace_back(encrypt(key, g1, coeffs[i] % mpr));
    }

    DEBUG << "Sending " << enc_poly.size() << " elements" << std::endl;

    //Scale the encrypted polynomial by lambda_i
    mpz_class scalar = rr.get_z_bits(400) % mpr;
    for(auto& i: enc_poly) i*=scalar;
    sock.send(enc_poly);
    total_communication += 3*enc_poly.size()*sizeof(Ec1)/3;
    t2.stop();
    std::cout << "Lambda test" << std::endl;
   
//====================================================================================================
//====================verfier check if lambda is computed correctly =============================
//====================================================================================================
    auto t3 = make_timer("Evaluate at lambda");
    //Evaluate at a random point lambda
    //sock.broadcast(poly_size, v_send, w_send, com_lamda);
    auto [poly_size_lamda,v_pk, w_pk, h1_pk,h2_pk,cku_pk, rounds_lambda, com_lamda] 
        = sock.receive<int,vector<Ec2>, vector<Ec1>, Ec2, Ec2, Ec1, int, Fp12>(6000s);

    auto [proof_L_lambda, proof_R_lambda] = sock.receive<std::vector<Fp12>, std::vector<Fp12>>(100s);

    auto [final_proof_1_lambda, final_proof_2_lambda] = sock.receive<Ec1, string>(100s);


    final_proof_1 = final_proof_1_lambda;
    final_proof_2 = mpz_class(final_proof_2_lambda);

    proof_L = proof_L_lambda;
    proof_R = proof_R_lambda;
    
    
    //receive chals
    auto [chals_lamda] = sock.receive<vector<string>>(100s);

    chals.resize(rounds_lambda);
    for (int i = 0; i < rounds_lambda; ++i)
    {
        chals[i] = mpz_class(chals_lamda[i]);
        cout<<chals[i] <<endl;

    }
    Verifier(g2,g1,poly_size_lamda,v_pk.data(),w_pk.data(),com_lamda,1,0);

//====================================================================================================
//====================================================================================================
  

    t3.stop();


   

   
    //====================wenxuan's malicious part ====================/
    cout<<"==============================Begin PHASE 1================================"<<endl;


    auto t4 = make_timer("Verifier time");


    auto [newsize, num_eval, deg_A, A, rounds, com_v] =
        sock.receive<int,int,
                     int, std::vector<Ec1>, int, Fp12>(6000s);


    //auto[com_Ti_recv] = sock.receive< vector<Ec1> >(5000s);
    init_proof(rounds);//assign the correct size to each proof vector

    h1 = h1_pk;
    h2 = h2_pk;
    cku = cku_pk;
    //receive proofs
    auto [proof_L_recv, proof_R_recv] = sock.receive<std::vector<Fp12>, std::vector<Fp12>>(100s);
    auto [final_proof_1_v,final_proof_2_v] = sock.receive<Ec1, string>(100s);


    final_proof_1 = final_proof_1_v;
    final_proof_2 = mpz_class(final_proof_2_v);

    proof_L = proof_L_recv;
    proof_R = proof_R_recv;
    
    
    //receive chals
    auto [chals_str] = sock.receive<vector<string>>(100s);
    chals.resize(rounds);
    for (int i = 0; i < rounds; ++i)
    {
        chals[i] = mpz_class(chals_str[i]);
        DEBUG<<chals[i]<<endl;
    }

    //gy = prover_set_broadcast_values_and_return_gy(newsize,g1,A,deg_A,eval_points,num_eval);//set com return gy
    DEBUG<<"begin verifier"<<endl;
    DEBUG<<"g1: "<<g1<<endl;
    DEBUG<<"g2: "<<g2<<endl;
    DEBUG<<"newsize: "<<newsize<<endl;
    DEBUG<<"num_eval: "<<num_eval<<endl;
    DEBUG<<"gy:"<<endl;
    
    std::cout << "Verifying proof" << std::endl;
    auto t_verfier_phase1= make_timer("Verifier phase1");
    // Verifier(g2, g1, newsize, v_pk.data(),
    //          eval_points.data(), com_v, gy_verfier.data(), num_eval, 0);
    Verifier(g2,g1,newsize,v_pk.data(),w_pk.data(),com_v,num_eval_pt,0);

    t_verfier_phase1.stop();
    cout<<"==============================END PHASE 1================================"<<endl;
    cout<<"==============================Begin PHASE 2================================"<<endl;

    //revceive the commitments
    // auto [com_Xi_recv,  com_Yi_recv, com_Zi_recv, com_Fi_recv, x_chal_str] =
    //     sock.receive<vector<Ec1>, vector<Ec1>,vector<Ec1>,vector<Ec1>,string>(6000s);
    //mpz_class chal_x = mpz_class(x_chal_str);
    // com_Xi = new Ec1[num_eval_pt];
    // com_Yi = new Ec1[num_eval_pt];
    // com_Zi = new Ec1[num_eval_pt];
    // com_Fi = new Ec1[num_eval_pt];

    // for (int i = 0; i < num_eval_pt; ++i)
    // {
    //     com_Xi[i] = com_Xi_recv[i];
    //     com_Yi[i] = com_Yi_recv[i];
    //     com_Zi[i] = com_Zi_recv[i];
    //     com_Fi[i] = com_Fi_recv[i];
    // }
//  sock.broadcast(Sa_send,Sb_send,Sc_send,Sd_send);
    // auto [Sa_recv,  Sb_recv, Sc_recv, Sd_recv] =
    //     sock.receive<vector<string>, vector<string>,vector<string>,vector<string>>(6000s);

    // Sa = new mpz_class[newsize];
    // Sb = new mpz_class[newsize];
    // Sc = new mpz_class[newsize];
    // Sd = new mpz_class[newsize];
    // for (int i = 0; i < newsize; ++i)
    // {
    //     Sa[i] = mpz_class(Sa_recv[i]);
    //     Sb[i] = mpz_class(Sb_recv[i]);
    //     Sc[i] = mpz_class(Sc_recv[i]);
    //     Sd[i] = mpz_class(Sd_recv[i]);

    // }
    //sock.broadcast(com_l_X,com_r_X,com_l_Y,com_r_Y,com_l_Z,com_r_Z,com_l_F,com_r_F);
    // auto [com_l_X_recv, com_r_X_recv, com_l_Y_recv, com_r_Y_recv, com_l_Z_recv, com_r_Z_recv, com_l_F_recv, com_r_F_recv] =
    //     sock.receive<Ec1,Ec1,Ec1,Ec1,Ec1,Ec1,Ec1,Ec1>(6000s);

    // com_l_X = com_l_X_recv;
    // com_r_X = com_r_X_recv;
    // com_l_Y = com_l_Y_recv;
    // com_r_Y = com_r_Y_recv;
    // com_l_Z = com_l_Z_recv;
    // com_r_Z = com_r_Z_recv;
    // com_l_F = com_l_F_recv;
    // com_r_F = com_r_F_recv;

    //Ec1 com_u = cku*(mie::Vuint(LR.get_str().c_str()));
    // mie::Vuint x(chal_x.get_str().c_str());
    // mie::Vuint x_inv = invert(chal_x);
    // mpz_class x_inv_mpz = invert_mpz_class(chal_x);

    // com_phase2_X = verifier_compute_com_phase2_X(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    // com_phase2_Y = verifier_compute_com_phase2_Y(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    // com_phase2_Z = verifier_compute_com_phase2_Z(newsize,num_eval_pt,x,x_inv,x_inv_mpz);
    // com_phase2_F = verifier_compute_com_phase2_F(newsize,num_eval_pt,x,x_inv,x_inv_mpz);



    auto [v2_recv,  com_phase2_X_recv, com_phase2_Y_recv, com_phase2_Z_recv, com_phase2_F_recv] =
        sock.receive<std::vector<Ec1>, Ec1,Ec1,Ec1,Ec1>(6000s);
    com_phase2_X = com_phase2_X_recv;
    com_phase2_Y = com_phase2_Y_recv;
    com_phase2_Z = com_phase2_Z_recv;
    com_phase2_F = com_phase2_F_recv;
    v_2 = new Ec1[newsize];
    for (int i = 0; i < newsize; ++i)
    {
        v_2[i] = v2_recv[i];
    }
    //receive proofs for X
    auto [proof_L_phase2_recv_X,  proof_R_phase2_recv_X, final_proof_1_phase2_str_X, final_proof_2_phase2_str_X] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>,string,string>(6000s);

    proof_L_phase2 = proof_L_phase2_recv_X;
    proof_R_phase2 = proof_R_phase2_recv_X;
    final_proof_1_phase2 = mpz_class(final_proof_1_phase2_str_X);
    final_proof_2_phase2 = mpz_class(final_proof_2_phase2_str_X);
    auto t6 = make_timer("Verifier phase2 x");
    Verifier_phase2(g2,g1,newsize,v_2,w_pk.data(),com_phase2_X,num_eval_pt,0);
    t6.stop();
    auto [proof_L_phase2_recv_Y,  proof_R_phase2_recv_Y, final_proof_1_phase2_str_Y, final_proof_2_phase2_str_Y] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>,string,string>(6000s);

    proof_L_phase2 = proof_L_phase2_recv_Y;
    proof_R_phase2 = proof_R_phase2_recv_Y;
    final_proof_1_phase2 = mpz_class(final_proof_1_phase2_str_Y);
    final_proof_2_phase2 = mpz_class(final_proof_2_phase2_str_Y);
    auto t7 = make_timer("Verifier phase2 y");
    Verifier_phase2(g2,g1,newsize,v_2,w_pk.data(),com_phase2_Y,num_eval_pt,0);
     t7.stop();
     auto [proof_L_phase2_recv_Z,  proof_R_phase2_recv_Z, final_proof_1_phase2_str_Z, final_proof_2_phase2_str_Z] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>,string,string>(6000s);

    proof_L_phase2 = proof_L_phase2_recv_Z;
    proof_R_phase2 = proof_R_phase2_recv_Z;
    final_proof_1_phase2 = mpz_class(final_proof_1_phase2_str_Z);
    final_proof_2_phase2 = mpz_class(final_proof_2_phase2_str_Z);
    auto t8 = make_timer("Verifier phase2 z");
    Verifier_phase2(g2,g1,newsize,v_2,w_pk.data(),com_phase2_Z,num_eval_pt,0);
    t8.stop();
    auto [proof_L_phase2_recv_F,  proof_R_phase2_recv_F, final_proof_1_phase2_str_F, final_proof_2_phase2_str_F] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>,string,string>(6000s);

    proof_L_phase2 = proof_L_phase2_recv_F;
    proof_R_phase2 = proof_R_phase2_recv_F;
    final_proof_1_phase2 = mpz_class(final_proof_1_phase2_str_F);
    final_proof_2_phase2 = mpz_class(final_proof_2_phase2_str_F);
    auto t9 = make_timer("Verifier phase2 F");
    Verifier_phase2(g2,g1,newsize,v_2,w_pk.data(),com_phase2_F,num_eval_pt,0);
    t9.stop();
    cout<<"==============================End PHASE 2================================"<<endl;
    cout<<"==============================Begin PHASE 3 (quadratic constraint)================================"<<endl;

    auto [com_Li_rec,  com_Ri_rec, w_prime_rec,cku_rec, x_phase3_str] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>, std::vector<Ec1>,Ec1,string>(6000s);
    x_phase3 = mpz_class(x_phase3_str);
    auto t11 = make_timer("Verifier phase 3");    
    com_Li = new Ec1[2*newsize];
    com_Ri = new Ec1[2*newsize];
    for (int i = 0; i < 2*newsize; ++i)
    {
        com_Li[i] = com_Li_rec[i];
        com_Ri[i] = com_Ri_rec[i];
    }

    auto [proof_coml_rec,  proof_comr_rec] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>>(6000s);
    proof_coml = proof_coml_rec;
    proof_comr = proof_comr_rec;

    com_d = cku_rec * ling;

    Verifier_Phase3(zero,com_d,x_phase3,newsize,newsize-1,0);
    cout<<"Verifier_Phase3 out"<<endl;

    Ec1 com_test = com_Li[0] + com_Ri[0] + com_d;
    cout<<com_test<<endl;
    auto [proof_L_phase2_recv_,  proof_R_phase2_recv_, final_proof_1_phase2_str_, final_proof_2_phase2_str_] =
        sock.receive<std::vector<Ec1>, std::vector<Ec1>,string,string>(6000s);

    proof_L_phase2 = proof_L_phase2_recv_;
    proof_R_phase2 = proof_R_phase2_recv_;
    final_proof_1_phase2 = mpz_class(final_proof_1_phase2_str_);
    final_proof_2_phase2 = mpz_class(final_proof_2_phase2_str_);

    Verifier_phase2(g2,g1,newsize,w_pk.data(),w_prime_rec.data(),com_test,num_eval,0);
    t11.stop();
    t4.stop();

    cout<<"==============================END PHASE 3================================"<<endl;
    //=======================================================================//
    //====================end wenxuan's malicious part ====================/
    //=======================================================================//

    std::cout << "Waiting to receive encrypted evaluations" << std::endl;
    auto tc1 = make_timer("Waiting for encrypted evaluations");

     //Decrypt and check if zero
    //
    auto [enc_eval] = sock.receive<std::vector<ciphertext>>(6000s);
    tc1.stop();
    std::cout << "Begin decryption" << std::endl;
    auto t10 = make_timer("Interactive Decryption");
    std::vector<ciphertext> reenc = rerandomize(key,g1,enc_eval);
    sock.send(reenc);
    total_communication += reenc.size() * sizeof(Ec1);
    auto [reenc_eval] = sock.receive<std::vector<ciphertext>>(100s);

    

    std::vector<Ec1> dprime, eprime;
    for (auto i = 0;i < reenc_eval.size(); i++)
    {
      dprime.emplace_back(reenc_eval[i][2] * mpzx);
    }
    sock.send(dprime);
    total_communication += dprime.size() * sizeof(Ec1)/3;
    auto [d] = sock.receive<std::vector<Ec1>>(100s);
    for (auto i = 0;i < d.size(); i++)
    {
      eprime.emplace_back(d[i] * mpzy - reenc_eval[i][0] * mpzy - reenc_eval[i][1] * mpzx);
    }
    sock.send(eprime);
    total_communication += eprime.size() * sizeof(Ec1)/3;

    auto [hresult] = sock.receive<bool>(100s);

    t10.stop();
    t.stop();
    show_timer();

    total_communication = total_communication + proof_size + commitment_size;
    cout<<"Communication Size: "<<total_communication<<endl;
    auto [up, down] = sock.communication_cost(storage_unit::Byte);

    std::cout << "Uploaded: " << up << " Bytes, " << std::endl
         << "Downloaded: " << down << " Bytes" << std::endl;
    INFO << "Done!";
}



int main(int argc, char *argv[])
{

    disable_logging();
    seed = rand();  
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
    p.set_str("21888242871839275222246405745257275088548364400416034343698204186575808495617",10);
    root_of_unity.set_str("12237580043724246686314448501040597564946944530769647865992922846663565257669", 10);
    bn::CurveParam cp = bn::CurveSNARK1;
    Param::init(cp);
    const Point& pt = selectPoint(cp);
    const Ec1 gt1(pt.g1.a, pt.g1.b);
    g1.set(pt.g1.a,pt.g1.b);
    g2.set(
        Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
        Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
    );



    namespace po = boost::program_options;

    po::options_description desc("Options");

    desc.add_options()
        ("client,c",   "Act like client")
        ("server,s",   "Act like server")
        ("duration,d", po::value<int>(), "Registration time")
    ("bind",        po::value<std::string>(), "Ip address and port the server will listen to. Default: 127.0.0.1:5555")
        ("file,f",    po::value<std::string>(), "Input file name");

    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        int duration_sec = 5;
        std::string path;
    
    std::string bind_addr("127.0.0.1");
    uint16_t bind_port = 5555;
        if (vm.count("client")) {
            if (vm.count("file")) { path = vm["file"].as<std::string>(); }

            std::fstream fs = psi::util::open_file(path);
            auto numbers = psi::util::read_line<mpz_class>(fs);

        if (vm.count("bind")) {
               std::vector<std::string> split_result;
               bind_addr = vm["bind"].as<std::string>();
               boost::split(split_result, bind_addr, [](char c){ return c == ':'; });

               if (split_result.size() > 1) {
                  bind_addr = std::move(split_result[0]);
                  bind_port = static_cast<uint16_t>(std::stoi(split_result[1]));
               }

               DEBUG << "Bind to " << bind_addr << ":" << bind_port;
            }


            run_client(std::move(numbers),bind_addr,bind_port);
        }
        else if (vm.count("server")) {
            if (vm.count("file")) { path = vm["file"].as<std::string>(); }

            std::fstream fs = psi::util::open_file(path);
            auto numbers = psi::util::read_line<mpz_class>(fs);
            if (vm.count("duration")) { duration_sec = vm["duration"].as<int>(); }


        if (vm.count("bind")) {
        std::vector<std::string> split_result;
        bind_addr = vm["bind"].as<std::string>();
        boost::split(split_result, bind_addr, [](char c){ return c == ':'; });

            if (split_result.size() > 1) {
              bind_addr = std::move(split_result[0]);
                      bind_port = static_cast<uint16_t>(std::stoi(split_result[1]));
        }

        DEBUG << "Bind to " << bind_addr << ":" << bind_port;
        
        }     
        run_server(std::chrono::seconds(duration_sec),std::move(numbers),bind_addr,bind_port);
        }
        else {
            std::cout << desc << std::endl;
        }
    }
    catch (...) {
        ERROR << boost::current_exception_diagnostic_information();
        return 1;
    }

   
}

/*
    a tiny sample of optimal ate pairing
*/

#include <thread>
#include <future>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>
#include <signal.h>

#include "bn.h"
#include "test_point.hpp"
// /#include "fft.h"
//#include "enviroment.h"

using namespace std;
using namespace bn;
static int errNum = 0;
unsigned long int seed;
mpz_class mod("21888242871839275222246405745257275088548364400416034343698204186575808495617");
mpz_class root_of_unity;

gmp_randstate_t r_state; //included in environment.h
mpz_class p;             //included in environment.h
Ec2* v;          //ck1
Ec1* w;         //ck2
Ec2 h1;         //ck3
Ec2 h2;         //ck3
mpz_class* b;    //batched eval vectors
int m,n,f_len,g_len;




//initialize all the proofs
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

//challenges x
vector<mpz_class> chals;



//Prover's broadcast msgs
Fp12 com;

//tools to calculate runtime
float eval_time = 0;
float pairing_time = 0;
float folding_time = 0;
float multieval_time = 0;
float singleeval_time = 0;

float folding_time_veri = 0;
float eval_time_veri = 0;
int test_round;


float Ec1_fft_time = 0;
//calculate f0/f1
Ec1* f0_;
Ec1* f1_;
int num_eval_pt; // = deg(f)


//some repetitive values
//mpz_class zero = mod;
mpz_class one = mpz_class(1);
mie::Vuint neg("21888242871839275222246405745257275088548364400416034343698204186575808495616");
mie::Vuint ling("21888242871839275222246405745257275088548364400416034343698204186575808495617");
Ec1 g1_p;
mpz_class zero(0);


//evaluation results
vector<Ec1> eval_results;


//randomness for batching
mpz_class* r;


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

vector<mpz_class> generate_random_challenge(int rounds){
    vector<mpz_class> x(rounds);

    mpz_class x_each;
    for(int i=0;i<rounds;i++){
        mpz_urandomm(x_each.get_mpz_t(),r_state,p.get_mpz_t());
        x[i] = x_each;
    }
    return x;
}
void init_proof(int rounds){
    proof_L.resize(rounds);
    proof_R.resize(rounds);
    proof_g_fL.resize(rounds);
    proof_g_fR.resize(rounds);

    proof_A.resize(rounds);
    proof_b.resize(rounds);
}
mpz_class exponetiate(mpz_class base, mpz_class pow){
        mpz_class result;
        //mpz_class mod("16798108731015832284940804142231733909759579603404752749028378864165570215949");
        mpz_powm(result.get_mpz_t(),base.get_mpz_t(),pow.get_mpz_t(),mod.get_mpz_t());
        return result;  
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
        cout<<poly[i]<<endl;
    }
}

Ec1 innerEc1Mpz(Ec1* vec1, mpz_class* vec2, int size){

    Ec1 innerProduct = vec1[0] * (mie::Vuint(vec2[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        innerProduct = innerProduct + vec1[i] * (mie::Vuint(vec2[i].get_str().c_str()));;
    }
    return innerProduct;
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
    for(int i=0;i<size*2;i++){
        eval[i] = eval1[i] * (mie::Vuint(eval2[i].get_str().c_str()));
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
mpz_class** polys;
mpz_class* getCoeff(mpz_class* roots, int m){
    polys = new mpz_class*[m];
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


void keygen_public(int size, Ec1 g1, Ec2 g2){
    //construct v,w,h1,h2

    mpz_class rand;
    //mpz_urandomm(rand.get_mpz_t(), r_state,rand.get_mpz_t());
    mpz_class alpha_mpz;
    mpz_class beta_mpz;
    mpz_class c1_mpz;
    mpz_class c2_mpz;


    //generating alpha and beta for public keys
    mpz_urandomm(alpha_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(beta_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(c1_mpz.get_mpz_t(),r_state,p.get_mpz_t());
    mpz_urandomm(c2_mpz.get_mpz_t(),r_state,p.get_mpz_t());

    //mie::Vuint alpha(alpha_mpz.get_str().c_str());

    v = new Ec2[size];
    w = new Ec1[size];
    for(int i=0;i<size;i++){
        mpz_class beta_each = exponetiate(beta_mpz,2*i);
        mpz_class alpha_each = exponetiate(beta_mpz,2*i);

        mie::Vuint beta(beta_each.get_str().c_str());
        mie::Vuint alpha(alpha_each.get_str().c_str());

        v[i] = g2 *beta;
        w[i] = g1 *alpha;    
    
    }

    mie::Vuint c1(c1_mpz.get_str().c_str());
    mie::Vuint c2(c2_mpz.get_str().c_str()); 
    h1 = g2*c1;
    h2 = g2*c2;
}

Ec1* dividePoly(Ec1* A, mpz_class *roots,int A_len, int deg_A,int num_roots, string side, Ec1 g1){
    Ec1* A_flip = flip_Poly_Ec1_Coeff(A,deg_A+1,A_len);
    clock_t f0 = clock();
    clock_t t = clock();


    g_len = num_roots+1;
    mpz_class*g = new mpz_class[g_len];
    g = getCoeff(roots, num_roots);


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
    //cout<<"reached here"<<endl;
    int l = ceil(log2(n-m+1));
    int h_size = pow(2,l);
    //int h_size = n-m+1;
    mpz_class * h = new mpz_class[h_size];
    h = compute_h(g_flip,g_len,n,m);
    //cout<<"after h"<<endl;
    


    //=========================================//
    //==========calcualte r (f0 or f1) ========//
    //=========================================//
    
    //compute q_flip = A_flip * h;   A_flip is cipher, but h is plain text
    int q_size = n-m+1;

    t = clock();
    Ec1* q_flip =  multiplyEc1_mpz(A_flip,A_len,h,h_size,g1);
    //cout<<"after h"<<endl;

    Ec1_fft_time += clock()-t;


    float q_flip_time = clock() - t;


    //compute q = flip(q_flip);
    Ec1 *q =flip_Poly_Ec1_Coeff(q_flip,q_size,q_size);
    //cout<<"here"<<endl;

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
void multiPointsEval(Ec1* A, mpz_class* points, int A_len,int deg_A,int num_points, int *indices, Ec1* output, Ec1 g1){

    // if(test_round!=1 and num_points==num_eval_pt/2){
    //     printf("%s\n","====================HALF WAY THRU!!!!!==================" );
    // }
    // cout<<"here"<<endl;
    // cout<<"deg_A: "<<deg_A<<endl;
    // cout<<"A_len: "<<A_len<<endl;
    // cout<<"num_points: "<<num_points<<endl;

    // exit(0);
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
       // Ec1* results = multiPointsEval_naive(A, eval_points, A_size, deg_A, num_eval_pt, g1);
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
        clock_t single = clock();
        Ec1* results = multiPointsEval_naive(A, eval_points, A_size, deg_A, num_eval_pt, g1);
        singleeval_time += clock() - single;
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
           clock_t single = clock();
           eval_FFT_batch = multiPointsEval_naive(A, eval_points_batch, A_size, deg_A, num_eval_batch, g1);
           singleeval_time += clock() - single;
        }
        else{//use double point eval;
            clock_t t = clock();
            multiPointsEval(A2, eval_points_batch, A_size2, deg_A, num_eval_batch, indices_test,eval_FFT_batch, g1);
            multieval_time += clock()-t;
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

Fp12 bbm_enc(Ec1* C, mpz_class* b, Ec1 d, Ec2* ck1, Ec1* ck2, int size){
    Fp12 com_bbm;
    opt_atePairing(com_bbm,  ck1[0],C[0]);
    for (int i = 1; i < size; i++)
    {
        Fp12 e;
        opt_atePairing(e, ck1[i],C[i]);
        com_bbm*=e;
    }
    Ec1 com_b = ck2[0]*(mie::Vuint(b[0].get_str().c_str()));
    for (int i = 1; i < size; ++i)
    {
        com_b = com_b + ck2[i]*(mie::Vuint(b[i].get_str().c_str()));
    }


    Fp12 e;
    opt_atePairing(e, h1,com_b);
    com_bbm *=e;
    opt_atePairing(e, h2, d);
    com_bbm *=e;
    return com_bbm;

}

Ec1* prover_set_broadcast_values_and_return_gy(int newsize,Ec1 g1,Ec1* A,  int deg_A, mpz_class* eval_points, int num_eval){
    

    //BBM.coommit initial

    //step 1: prover compute the eval vectors

    mpz_class** T = new mpz_class*[num_eval]; //all eval vectors
    for (int i = 0; i < num_eval; ++i)
    {
        T[i] = new mpz_class[num_eval];
        for (int j = 0; j < newsize; ++j)
        {
            T[i][j] = exponetiate(eval_points[i],j);
        }
    }

    //batch all the T
    r = new mpz_class[num_eval];//randomness for batching
    b = new mpz_class[num_eval];
    for (int i = 0; i < num_eval; ++i)
    {
        int val = rand()%100+1;
        r[i] = mpz_class(val);
    }

    for (int i = 0; i < num_eval; ++i)
    {
        for (int j = 0; j < newsize; ++j)
        {
            T[i][j] = T[i][j] * r[i];
            b[j] = (b[j] + T[i][j]);
        }
    }
    for (int i = 0; i < newsize; ++i)
    {
        b[i] = b[i]%mod;
    }
    opt_atePairing(com,v[0],A[0]);
    for (int i = 1; i < newsize; ++i)
    {
        Fp12 e;
        opt_atePairing(e, v[i],A[i]);
        com *=e;
    }

    Ec1 com_b = w[0]*(mie::Vuint(b[0].get_str().c_str()));
    for (int i = 1; i < newsize; ++i)
    {
        com_b = com_b + w[i]*(mie::Vuint(b[i].get_str().c_str()));
    }
    // Ec2* v;          //ck1
    // Ec1* w;         //ck2
    // Ec2 h1;         //ck3
    // Ec2 h2;         //ck3
    
    //eval polynomials on all eval points
    Ec1* enc_f = compute_gf(A,deg_A,newsize,eval_points,num_eval,g1);  
    Ec1 d; //encrypted evaluation on batched vector b.
    d = enc_f[0] * (mie::Vuint(r[0].get_str().c_str()));
    for (int i = 1; i < newsize; ++i)
    {
        d = d + enc_f[i] * (mie::Vuint(r[i].get_str().c_str()));
    }


    //complete committing com

    Fp12 e;
    opt_atePairing(e, h1,com_b);
    com *=e;
    opt_atePairing(e, h2,d);
    com *=e;    
    // cout<<com<<endl;

    // cout<<"=========================="<<endl;

    // Fp12 com2 = bbm_enc(A,b,d,v,w,newsize);
    // cout<<com2<<endl;



    //compute y and gy
  
    int *indices = new int [num_eval];
    for (int i = 0; i < num_eval; ++i)
    {
        indices[i] = i;
    }



    Ec1* gy_broadcast = compute_gf(A,deg_A,newsize,eval_points,num_eval,g1);  
    return gy_broadcast;

    //return gy_broadcast;
}


void Prover(Ec2 g2, Ec1 g1,int d,int deg_A, Ec1* A ,Ec2* v, Ec1* w, mpz_class* eval_points, mpz_class* b, int num_eval,int round){
    printf("=============Prover iteration : %d==================\n",round);

    if(d>=2){
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

      
        // cl = innerEc1Mpz(A_R,b_L,d_prime);
        // cr = innerEc1Mpz(A_L,b_R,d_prime);


        cl = innerEc1Mpz(A_L,b_R,d_prime);
        cr = innerEc1Mpz(A_R,b_L,d_prime);

        Fp12 L, R;
        //compute L = bbm.encrytion

        // L = bbm_enc(A_R, b_L,cl,v_L,w_R,d_prime);
        // R = bbm_enc(A_L, b_R,cr,v_R,w_L,d_prime);


        L = bbm_enc(A_L, b_R,cl,v_R,w_L,d_prime);
        R = bbm_enc(A_R, b_L,cr,v_L,w_R,d_prime);

        //store the proofs
        proof_L[round] = L;
        proof_R[round] = R;


        //grab challenge from chals;
        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mpz_class x_inv_mpz = invert_mpz_class(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());





        //"fold" the secret vector A,b and v,w
        time = clock();
        Ec1* A_prime = new Ec1[d_prime];
        Ec2* v_prime = new Ec2[d_prime];
        mpz_class* b_prime = new mpz_class[d_prime];
        Ec1* w_prime = new Ec1[d_prime];



        for (int i = 0; i < d_prime; i++)
        {   
            A_prime[i] = A_L[i]*x + A_R[i]*x_inv;
            v_prime[i] = v_L[i]*x_inv + v_R[i]*x;

            b_prime[i] = b_L[i] * x_inv_mpz + b_R[i]*x_temp;
            w_prime[i] = w_L[i] * x + w_R[i]*x_inv;

        }
        //store debug proof

        proof_A[round] = A_prime;
        proof_b[round] = b_prime;
       
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
        Prover(g2,g1,d_prime, deg_A,A_prime,v_prime,w_prime,eval_points,b_prime,num_eval,round);
        return;
        
    }
    else{
        //create final proof
        cout<<"here"<<endl;
        final_proof_1 = A[0];
        final_proof_2 = b[0];
        cout<<"done"<<endl;

        return;
    }
}
//Verifier(g2,g1,newsize,v,w,eval_points,com,gy,num_eval_pt,0);

void Verifier(Ec2 g2, Ec1 g1, int d, Ec2* v, Ec1* w, mpz_class* eval_points, Fp12 com,Ec1* g_y, int num_eval, int round){
    
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
            v_prime[i] = v_L[i]*x_inv + v_R[i]*x;

            w_prime[i] = w_L[i]*x + w_R[i]*x_inv;
        }

        //printCipher2Poly(v,d);
        //printCipherPoly(w_prime,d_prime);
        folding_time_veri += (double)(clock()-time)/CLOCKS_PER_SEC;

        Fp12 com_prime;
        com_prime = power(L, x_2) * com * power(R, x_minu2); 

        //=-===============begin degbug====================
        Ec1* A_prime = proof_A[round];
        mpz_class* b_prime = proof_b[round];
        
  
        Ec1 innerpro = innerEc1Mpz(A_prime,b_prime,d_prime);
      
        Ec1* A_test = new Ec1[d];
        Ec1* A_left = new Ec1[d_prime];
        Ec1* A_right = new Ec1[d_prime];
        for (int i = 0; i < d_prime; ++i)
        {
            A_left[i] = A_prime[i]*x_inv;
            A_right[i] = A_prime[i]*x;
            A_test[i] = A_left[i];
            A_test[i+d_prime] = A_right[i];
        }


        mpz_class* b_test = new mpz_class[d];
        mpz_class* b_left = new mpz_class[d_prime];
        mpz_class* b_right = new mpz_class[d_prime];
        for (int i = 0; i < d_prime; ++i)
        {
            b_left[i] = (b_prime[i]*x_temp)%mod;
            b_right[i] = (b_prime[i]*x_inv_mpz)%mod;
            b_test[i] = b_left[i];
            b_test[i+d_prime] = b_right[i];
        }

  

        Fp12 com_check = bbm_enc(A_prime,b_prime,innerpro,v_prime,w_prime,d_prime);
        Fp12 com_test = bbm_enc(A_test,b_test,innerpro,v,w,d);



        printf("%s\n","L com R ============================" );
        cout<<com_prime<<endl;
        printf("%s\n","com computed============================" );
        cout<<com_test<<endl;
        //=-===============end debug====================

        exit(1); 

      
        
        delete[] v_L;
        delete[] v_R;

        //recursion
        round = round +1;
        Verifier(g2, g1,d_prime, v_prime,w_prime,eval_points, com_prime, g_y,num_eval, round);
        return;
    }
    else{//check com = bbm.commit(C,b,<C,b>)
        cout<<"============================ Verifier Final iteration ============================"<<endl;
        cout<<com<<endl;
        
        Ec1 final_inner_prod = final_proof_1*(mie::Vuint(final_proof_2.get_str().c_str()));
        Ec1* final_proof_1_vec = new Ec1[1];
        mpz_class* final_proof_2_vec = new mpz_class[1];
        final_proof_1_vec[0] = final_proof_1;
        final_proof_2_vec[0] = final_proof_2;


        Fp12 com_2 = bbm_enc(final_proof_1_vec,final_proof_2_vec,final_inner_prod,v,w,1);

        printf("%s\n","================================================" );
        cout<<com_2<<endl;

    }
}
int main(int argc, char *argv[])
{


    seed = rand();  
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
    p.set_str("21888242871839275222246405745257275088548364400416034343698204186575808495617",10);
    root_of_unity.set_str("12237580043724246686314448501040597564946944530769647865992922846663565257669", 10);
    bn::CurveParam cp = bn::CurveSNARK1;
    Param::init(cp);
    const Point& pt = selectPoint(cp);
    const Ec2 g2(
        Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
        Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
    );
    const Ec1 g1(pt.g1.a, pt.g1.b);
    g1_p = g1*ling;
    mpz_class zero(0);

    //=========================================//
    //=====Init A, a and some sizes============//
    //=========================================//
    int log_f_pow;
    sscanf(argv[1], "%d", &log_f_pow);

    
    cout<<"log_deg(f): "<<log_f_pow<<endl;
    int log_g_pow = log_f_pow-1;
    
    n = pow(2,log_f_pow); //deg(f)
    m = pow(2,log_g_pow); //deg(g)    n >= m


    f_len = n + 1;


    //n = log_f_pow;
    n--;
    int deg_A = n;
    // g_len = m + 1;

     //initialize f(x)  or A in the malicious protocol
    int coeff_len = n+1;
    vector<Ec1> A_(coeff_len);
    vector<mpz_class> a_(coeff_len);

    for (int i = 0; i < coeff_len; ++i)//for an easy example, our polynomial's coeff are within 100
    {
        int val = rand()%100+1;
        A_[i] = g1*(val);
        a_[i] = mpz_class(val); 
    }


    //if the vector is not a power of 2, append 0s
    int newsize = coeff_len;
    if(isPowerOfTwo(coeff_len)==0){
        int log_newsize = ceil(log2(coeff_len));
        newsize = pow(2,log_newsize);
        A_.resize(newsize);
        a_.resize(newsize);
    }
    cout<<"newsize: "<<newsize<<endl;
    int oldsize = newsize;
    Ec1* A = new Ec1[newsize];// this is the true f(x)
    mpz_class* a = new mpz_class[newsize];
    for (int i = 0; i < newsize; ++i)
    {
        if(i<oldsize){
            A[i] = A_[i];
            a[i] = a_[i];
        }else{
            A[i] = g1*mie::Vuint(0);
            a[i] = zero;
        }
    }


    //generate random eval points
    //num_eval_pt = 8; // = deg(f)
    num_eval_pt = n+1; // = deg(f)

    mpz_class* eval_points = new mpz_class[num_eval_pt];
    for (int i = 0; i < num_eval_pt; ++i)
    {
        eval_points[i] = p - rand();
    }

    //=========================================//
    //==================MPSI Begin=============//
    //=========================================//
    mie::Vuint test("16798108731015832284940804142231733909759579603404752749028378864165570215949");


    cout<<"polysize: "<<newsize<<endl;
    cout<<"num points: "<<num_eval_pt<<endl;

    keygen_public(newsize,g1, g2);
    //x: all the chals needed for this protocol
    int rounds = log2(newsize);//rounds = how many chals or proofs

    //generate challenges before hand
    chals = generate_random_challenge(rounds);  
    init_proof(rounds);//assign the correct size to each proof vector
    


    //eval all the points, compute batch eval points
    Ec1* gy = prover_set_broadcast_values_and_return_gy(newsize,g1,A,deg_A,eval_points,num_eval_pt);//set com return gy
    cout<<"gy g_y_computed"<<endl;
    //Run Prover

    clock_t t_prover = clock();
    Prover(g2,g1,newsize,deg_A,A,v,w,eval_points,b,num_eval_pt,0);
    float prover_time = clock() - t_prover;
    cout<<"prover time: "<<(double)(prover_time)/CLOCKS_PER_SEC<<"s\n";
    cout<<"prover pairing time: "<<pairing_time<<"s\n";
    cout<<"prover evaluation time: "<<eval_time<<"s\n";
    cout<<"prover Eval_FFT time: "<<(double)multieval_time/CLOCKS_PER_SEC<<"s\n";
    cout<<"prover Eval_naive time: "<<(double)singleeval_time/CLOCKS_PER_SEC<<"s\n";
    cout<<"prover fold time: "<<folding_time<<"s\n";

    //float other_time =prover_time - pairing_time - eval_time - folding_time; 
    //cout<<"prover other time:"<<other_time<<"s\n";
    printf("%s\n\n\n\n","");


    //Run Verifier (assume verifier has the evaluation point t)
    
    clock_t t_verifier = clock();
    Verifier(g2,g1,newsize,v,w,eval_points,com,gy,num_eval_pt,0);
    cout<<"verfier time: "<<(double)(clock()-t_verifier)/CLOCKS_PER_SEC<<"s\n";
    cout<<"verfier evaluation time: "<<eval_time_veri<<"s\n";
    cout<<"verfier fold time: "<<folding_time_veri<<"s\n";
    
}

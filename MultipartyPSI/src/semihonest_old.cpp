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
#include <chrono>

#include <boost/program_options.hpp>

#include <ligero/tcp/server.hpp>
#include <ligero/tcp/client.hpp>
#include <ligero/network/serialize.hpp>

#define MIE_ATE_USE_GMP 1


using namespace ligero::net;
using namespace std::chrono_literals;
using namespace std::string_literals;

using serializer_t = serializer<boost_portable_binary>;
using tcp_server = basic_tcp_server<serializer_t>;
using tcp_client = basic_tcp_client<serializer_t>;

using namespace std;
using namespace bn;
static int errNum = 0;
unsigned long int seed;
mpz_class mod("21888242871839275222246405745257275088548364400416034343698204186575808495617");
mpz_class root_of_unity;

gmp_randstate_t r_state; //included in environment.h
mpz_class p;             //included in environment.h
Ec2* v;          //pk
int m,n,f_len,g_len;
Ec1 g1;
 



//initialize all the proofs
vector<Fp12> proof_L;
vector<Fp12> proof_R;
vector<Ec1*> proof_g_fL;
vector<Ec1*> proof_g_fR;

vector<Ec1> p_;
vector<Ec1> r_;

Ec1 final_proof;
//challenges x
vector<mpz_class> chals;



//Prover's broadcast msgs
Fp12 com;

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


//some repetitive values
//mpz_class zero = mod;
mpz_class one = mpz_class(1);
mie::Vuint neg("21888242871839275222246405745257275088548364400416034343698204186575808495616");
mie::Vuint ling("21888242871839275222246405745257275088548364400416034343698204186575808495617");
Ec1 g1_p;
mpz_class zero(0);

template<class T, class S>
void verify(const char *msg, const T& a, const S& b)
{
	if (a == b) {
		printf("%s : ok\n", msg);
	} else {
		printf("%s : ng\n", msg);
		PUT(a);
		PUT(b);
		errNum++;
	}
}


// Client / Server (Party P1) code
//
//

void run_server(std::chrono::seconds accept_time) {

    /*  Initialize a executor */
    asio::io_context executor;

    /*  Construct a TCP server */
    tcp_server sock(executor);

    /*  Optional: set worker threads */
    sock.threads(2);

    DEBUG << "Binding to 127.0.0.1:5555";

    /*  Bind the server to localhost:5555 */
    auto accepted = sock.bind("127.0.0.1", 5555, accept_time);

    INFO << "Accepted " << accepted << " connections";

    /*  Notify clients to start */
    sock.broadcast("start"s);

    /*  Deserialize received messages and send back */
    {
        auto [msgs, kickedout] = sock.collect(1s);

        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
       
	Ec1 Pusum,Pvsum,Phsum,Ptemp;
	{	
		for (auto& msg : msgs) {
		    /*  Deserialize using server's archive. Alternatively, one could write:
		     *  
		     *  int a, b, c;
		     *  sock.archive().unpack(msg, a, b, c);
		     */
		    auto [pktemp] = sock.archive().unpack<string>(msg);
		    std::stringstream temp;
		    temp.str(pktemp);
		    temp >> Ptemp;

		    Pusum = Pusum + Ptemp;

		    std::cout << "received g^u = " << pktemp <<  std::endl;
		    
		}
		DEBUG << "Seding Public Key u: " << Pusum << std::endl ;

		stringstream temp;
		temp << Pusum;
		string stemp;
		temp >> stemp;
		sock.broadcast(stemp);

	}

	{
		for (auto& msg : msgs) {
		    /*  Deserialize using server's archive. Alternatively, one could write:
		     *  
		     *  int a, b, c;
		     *  sock.archive().unpack(msg, a, b, c);
		     */
		    auto [pktemp] = sock.archive().unpack<string>(msg);
		    std::stringstream temp;
		    temp.str(pktemp);
		    temp >> Ptemp;

		    Pvsum = Pvsum + Ptemp;

		    std::cout << "received g^v = " << pktemp <<  std::endl;
		    
		}
		DEBUG << "Sending Public Key v: " << Pvsum << std::endl ;

		stringstream temp;
		temp << Pvsum;
		string stemp;
		temp >> stemp;
		sock.broadcast(stemp);
	}

	{
		for (auto& msg : msgs) {
		    /*  Deserialize using server's archive. Alternatively, one could write:
		     *  
		     *  int a, b, c;
		     *  sock.archive().unpack(msg, a, b, c);
		     */
		    auto [pktemp] = sock.archive().unpack<string>(msg);
		    std::stringstream temp;
		    temp.str(pktemp);
		    temp >> Ptemp;

		    Phsum = Phsum + Ptemp;

		    DEBUG << "received g^h = " << pktemp <<  std::endl;
		    
		}
		DEBUG << "Sending Public Key h: " << Phsum << std::endl ;

		stringstream temp;
		temp << Phsum;
		string stemp;
		temp >> stemp;
		sock.broadcast(stemp);
	}
    }

    INFO << "Done!";
}

// ------------------------------------------------------------

void run_client(int input) {

    mpz_class mpr = mpz_class((Param::r).toString());
    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));
/*    mpz_class ran;
    ran =rr.get_z_bits(400);
    std::cout<<"ran before = " << ran.get_str() << std::endl;
    ran = ran % mpr;
    std::cout<<"ran after  = " << ran.get_str() << std::endl;
    Ec1 Ptemp = g1 * ran;
    PUT(Ptemp);
*/

    /*  Additional thread hint, optimized for single thread */
    asio::io_context executor(1);

    tcp_client sock(executor);

    /*  Connect to remote server  */
    sock.connect("127.0.0.1", "5555");

    /*  Wait until server is ready... */
    auto [str] = sock.receive<std::string>(10s);


    const char *xStr = "1456789012345";
    const char *yStr = "942273498237498237498237489322342342424242421";
    const char *wStr = "94238749237492374982734982734982739487298472398472938472389472389742837428937429837429834729837428937429837429837429837492837492837498237489237489237498273498237498237498237489322342342424242421";
    const char *zStr = "28409802834209384202";

    mpz_class mpzx = rr.get_z_bits(400);
    mpzx = mpzx % mpr;
    mpz_class mpzy = rr.get_z_bits(400);
    mpzy = mpzy % mpr;
    mpz_class mpzz = rr.get_z_bits(400);
    mpzz = mpzz % mpr;


    const mie::Vuint x(mpzx.get_str());
    const mie::Vuint y(mpzy.get_str());
    const mie::Vuint z(mpzz.get_str());
	
    Ec1 Pu = g1 * x;
    Ec1 Pv = g1 * y;
    Ec1 Ph = g1 * z;

    /*
    const char *rStr = "1456789012345";
    const char *sStr = "92342342424242421";
    const mie::Vuint r(rStr);
    const mie::Vuint s(sStr);

    const char *mStr = "0";
    const mie::Vuint m(mStr);
    Ec1 Pa = Pu * r;
    Ec1 Pb = Pv * s;
    Ec1 Pc = (g1 * (r+s)) + (Ph * m);

    const char *tStr = "123456789012345";
    const mie::Vuint t(tStr);
    Ec1 Paa = Pa * t;
    Ec1 Pbb = Pb * t;
    Ec1 Pcc = Pc * t;

    Ec1 Pdd = ((Pc * x) * y); 
    Ec1 Pee = Pdd - Pb * x - Pa * y; 

    Ec1 temp = (Pc - g1) * (r+s);
    std::cout << "Value = " << temp << std::endl;
    std::cout << "Value = " << ((Ph * m) * (r+s)) << std::endl;




    std::cout << "Encryption is " << ((((Ph * m) * x) * y)) << std::endl;
    std::cout << "Decryption = " << Pee << std::endl;

    PUT(Ph * Param::r);
    PUT(Ph * Param::p);
    */


    Ec1 Ppk_u, Ppk_v, Ppk_h;
    DEBUG << "Server says: " << str;

    {
        //sock.send(input);
	string pktemp;
	stringstream temp;
	temp << Pu;
	temp >> pktemp;
	DEBUG << "Sending g^u = " << pktemp << std::endl;
	sock.send(pktemp);
        auto [sum] = sock.receive<string>(10s);
	stringstream temp2;
	temp2.str(sum);
	temp2 >> Ppk_u;
	DEBUG << "Public key u = " << sum << std::endl;
    }

    {
	string pktemp;
	stringstream temp;
	temp << Pv;
	temp >> pktemp;
	DEBUG << "Sending g^v = " << pktemp << std::endl;
	sock.send(pktemp);
        auto [sum] = sock.receive<string>(10s);
	stringstream temp2;
	temp2.str(sum);
	temp2 >> Ppk_v;
	DEBUG << "Public key v = " << sum << std::endl;
    }


    {
	string pktemp;
	stringstream temp;
	temp << Ph;
	temp >> pktemp;
	DEBUG << "Sending g^h = " << pktemp << std::endl;
	sock.send(pktemp);
        auto [sum] = sock.receive<string>(10s);
	stringstream temp2;
	temp2.str(sum);
	temp2 >> Ppk_h;
	DEBUG << "Public key h = " << sum << std::endl;
    }


    INFO << "Done!";
}




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
        //cout<<poly[i]%mod<<endl;
    }
}
void printCipherPoly(Ec1* poly, int size){
    for (int i = 0; i < size; i++)
    {
        //cout<<poly[i]<<endl;
    }
}
void printCipherPoly_vec(vector<Ec1> poly, int size){
    for (int i = 0; i < size; i++)
    {
        //cout<<poly[i]<<endl;
    }
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

Ec1* prover_set_broadcast_values_and_return_gy(int newsize,Ec1 g1,Ec1* A,  int deg_A, mpz_class* eval_points, int num_eval){
    //compute com = <A,v> 


    opt_atePairing(com,  v[0],A[0]);
    for (int i = 1; i < newsize; i++)
    {
        Fp12 e;
        opt_atePairing(e, v[i],A[i]);
        com*=e;
    }
    //compute y and gy
  
    int *indices = new int [num_eval];
    for (int i = 0; i < num_eval; ++i)
    {
        indices[i] = i;
    }
    Ec1* gy_broadcast = new Ec1[num_eval];

    
    multiPointsEval(A, eval_points,newsize,deg_A,num_eval,indices, gy_broadcast,g1);
    return gy_broadcast;
}

Ec1* compute_gf(Ec1* A, int deg_A,int A_size ,mpz_class* eval_points, int num_eval_pt,Ec1 g1){

    Ec1* results = new Ec1[num_eval_pt];
    if(num_eval_pt>deg_A){
        int iteration = num_eval_pt/A_size;
        int num_eval_each = num_eval_pt/iteration;
      
        int *indices_test = new int [num_eval_each];
        for (int i = 0; i < num_eval_each; ++i)
        {
            indices_test[i] = i;
        }

        for (int i = 0; i < iteration; ++i)
        {
            Ec1 *results_each = new Ec1[num_eval_each];
            mpz_class *eval_point_each = new mpz_class[num_eval_each];

            for (int j = 0; j < num_eval_each; ++j)
            {
                eval_point_each[j] = eval_points[j+i*num_eval_each];
                
            }
            if(deg_A>0){
                multiPointsEval(A, eval_point_each, A_size,deg_A,num_eval_each,indices_test,results_each, g1);
            }
            else{
                mpz_class idx;
                for(int k = 0;k<num_eval_each;k++){

                    Ec1 result = A[0]*(mie::Vuint(exponetiate(eval_point_each[k],zero).get_str().c_str()));
                    for (int l = 1; l < A_size; l++)
                    {
                        idx = mpz_class(l);
                        result += A[l]*(mie::Vuint(exponetiate(eval_point_each[k],idx).get_str().c_str()));
                    }
                    results_each[k] =  result;

                }
            }
            for (int j = 0; j < num_eval_each; ++j)
            {
                results[j+i*num_eval_each] = results_each[j];
                
            }
        }
    }
    return results;
}
void Prover(Ec2 g2, Ec1 g1,int d,int deg_A, Ec1* A ,Ec2* v, mpz_class* eval_points, int num_eval,int round){
    printf("=============Prover iteration : %d==================\n",round);
    if(d>=2){
        //cout<<"round: "<<round<<endl;
        int d_prime = d/2;
        //split the secret vector into 2
        clock_t time = clock();

        Ec1* A_L = new Ec1[d_prime];
        Ec1* A_R = new Ec1[d_prime];

        for (int i = 0; i < d_prime; i++)
        {
            A_L[i] = A[i];
            A_R[i] = A[i+d_prime];
        }

        //split the pk v in to 2
        Ec2* v_L = new Ec2[d_prime];
        Ec2* v_R = new Ec2[d_prime];
    
        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];
        }
        Fp12 L, R;
        //compute L = <A_L,v_R>

        opt_atePairing(L, v_R[0],A_L[0]);
        for (int i = 1; i < d_prime; i++)
        {
            Fp12 e;
            opt_atePairing(e, v_R[i],A_L[i]);
            L*=e;
        }

        //compute R
        opt_atePairing(R, v_L[0],A_R[0]);
        for (int i = 0; i < d_prime; i++)
        {
            Fp12 e;
            opt_atePairing(e,v_L[i], A_R[i]);
            R*=e;
        }
        pairing_time += (double)(clock()-time)/CLOCKS_PER_SEC;
    

        //compute g_fL, g_fR
        time = clock();
        int deg_A_l = d/2-1;
        int deg_A_r = deg_A-deg_A_l-1;
        Ec1* g_fL = compute_gf(A_L, deg_A_l, d_prime ,eval_points,  num_eval, g1);
        Ec1* g_fR = compute_gf(A_R, deg_A_r, d_prime ,eval_points,  num_eval, g1);
        eval_time += (double)(clock()-time)/CLOCKS_PER_SEC;

        //create the proofs
        proof_L[round] = L;
        proof_R[round] = R;
        proof_g_fL[round] = g_fL;
        proof_g_fR[round] = g_fR;



        //grab challenge from chals;
        mpz_class x_temp = chals[round];
        
        mie::Vuint x_inv = invert(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());




        //"fold" the secret vector A and pk v
        time = clock();
        Ec1* A_prime = new Ec1[d_prime];
        Ec2* v_prime = new Ec2[d_prime];
        Ec1 AR_x; //AR_x = A_R ^ (x^-1)
        Ec2 vR_x;
    

        for (int i = 0; i < d_prime; i++)
        {   
            AR_x = A_R[i]*x_inv;
            A_prime[i] = A_L[i] + AR_x;

            vR_x = v_R[i]*x;
            v_prime[i] = v_L[i] + vR_x;
        }
        folding_time += (double)(clock()-time)/CLOCKS_PER_SEC;

        //recurse 
        round = round + 1;
        deg_A = d_prime-1;    
        Prover(g2,g1,d_prime, deg_A,A_prime,v_prime,eval_points,num_eval,round);
        //printf("Return: %s\n", "round");
        return;
        
    }
    else{
        //create final proof
        final_proof = A[0];
        return;
    }
}

void Verifier(Ec2 g2, Ec1 g1, int d, Ec2* v,  mpz_class* eval_points, Fp12 com,Ec1* g_y, int num_eval, int round){
    //printf("=============Verifier iteration : %d==================\n",round);

    if(d>=2){
        int d_prime = d/2;
        //receive msgs from the prover
        Fp12 L = proof_L[round];
        Fp12 R = proof_R[round];
        Ec1 *g_fL = proof_g_fL[round];
        Ec1 *g_fR = proof_g_fR[round];


        clock_t time = clock();
        //compute v_prime
        Ec2* v_prime = new Ec2[d_prime];
        Ec2* v_L = new Ec2[d_prime];
        Ec2* v_R = new Ec2[d_prime];

        for (int i = 0; i < d_prime; i++)
        {
            v_L[i] = v[i];
            v_R[i] = v[i+d_prime];
        }

        Ec2 vR_x;
        mpz_class x_temp = chals[round];
        mie::Vuint x_inv = invert(x_temp);
        mie::Vuint x(x_temp.get_str().c_str());

        for (int i = 0; i < d_prime; i++)
        {
            vR_x = v_R[i]*x;
            v_prime[i] = v_L[i] + vR_x;
        }

        folding_time_veri += (double)(clock()-time)/CLOCKS_PER_SEC;

        Fp12 com_prime;

        com_prime = power(L, x) * com * power(R, x_inv); 

        time = clock();
        //verfier can compute r = (g_fR)^(t^d')
        mpz_class d_prime_mpz(to_string(d_prime).c_str());



        Ec1* g_y_prime = new Ec1[num_eval];
        //verify every single eval points
        int correct_count = 0;
        for (int i = 0; i < num_eval; ++i)
        {
            mpz_class t_mpz = eval_points[i];
            mie::Vuint t_d_prime = mie::Vuint(exponetiate(t_mpz,d_prime_mpz).get_str().c_str());
            Ec1 r = g_fR[i] * t_d_prime;
            Ec1 g_y_computed = g_fL[i]+r;
            eval_time_veri += (double)(clock()-time)/CLOCKS_PER_SEC;

            //perform round verification
            if(g_y[i]==g_y_computed){
                correct_count++;
            }
            else{
                printf("Verifier iteration %d FAIL\n", round);
                exit(1);
            }
            
            g_y_prime[i] = g_fL[i] + g_fR[i]*x_inv;
        }
        if(correct_count==num_eval){
            printf("Verifier iteration %d PASS\n", round);
        }
        else{
            printf("Verifier iteration %d FAIL\n", round);
            
        }
        
        

        //recursion
        round = round +1;
        Verifier(g2, g1,d_prime, v_prime,eval_points, com_prime, g_y_prime,num_eval, round);
        return;
    }
    else{
        //Final round: verify ga = gy
        int correct_count = 0;

        for (int i = 0; i < num_eval; ++i)
        {
            if(g_y[i]==final_proof){
                correct_count++;
            }else{
                printf("Verifier Final check FAIL\n");
                exit(1);
            }            
        }
        if(correct_count==num_eval){
           printf("Verifier Final check PASS\n");
        }else{
            printf("Verifier Final check FAIL\n");
        }      
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
    const Ec1 gt1(pt.g1.a, pt.g1.b);
    g1.set(pt.g1.a,pt.g1.b);
    g1_p = g1*ling;
    mpz_class zero(0);

    namespace po = boost::program_options;

    po::options_description desc("Options");
    desc.add_options()
        ("client,c",   "Act like client")
        ("server,s",   "Act like server")
        ("duration,d", po::value<int>(), "Registration time")
        ("input,i",    po::value<int>(), "Input value");

    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        int duration_sec = 5, input = 1;

        if (vm.count("client")) {
            if (vm.count("input")) { input = vm["input"].as<int>(); }
            run_client(input);
        }
        else if (vm.count("server")) {
            if (vm.count("duration")) { duration_sec = vm["duration"].as<int>(); }
            run_server(std::chrono::seconds(duration_sec));
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

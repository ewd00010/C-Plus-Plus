/**
 * @file
 * @brief finds the Greatest common divisor using [extended Euclidean algorithm]
 * (https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm)
 *
 * @details
 * Finding coefficients of a and b ie x and y in  Bézout's identity
 * \f[\text{gcd}(a, b) = a \times x + b \times y \f]
 * This is also used in finding the Modular
 * multiplicative inverse (MMI) of a number. \f[(A * B)%M == 1\f] Here B is the MMI of A for
 * given M, so extendedEuclid (A, M) gives B.
 */

#include <algorithm>  // for swap function
#include <assert.h> // for assert
#include <iostream> // for IO operations

/**
 * @brief Mathematical algorithms
 * @namespace
 */
namespace math {

    /**
    * @brief Function to update the coefficients per iteration
    * \f[r_0,\,r = r,\, r_0 - \text{quotient}\times r\f]
    * @param r
    * @param r0
    * @param quotient
    */
    inline void update_step(uint64_t *r, uint64_t *r0, const uint64_t quotient) {
        uint64_t temp = *r;
        *r = *r0 - (quotient * temp);
        *r0 = temp;
    }

    /**
    * @brief Implementation using iterative algorithm from
    * [Wikipedia](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode)
    * @param A
    * @param B
    * @param GCD
    * @param x
    * @param y
    */
    void extendedEuclid_1(uint64_t A, uint64_t B, uint64_t *GCD, int64_t*x, int64_t *y) {
        if (B > A)
            std::swap(A, B);  // Ensure that A >= B

        uint64_t s = 0, s0 = 1;
        uint64_t t = 1, t0 = 0;
        uint64_t r = B, r0 = A;

        while (r != 0) {
            uint64_t quotient = r0 / r;
            update_step(&r, &r0, quotient);
            update_step(&s, &s0, quotient);
            update_step(&t, &t0, quotient);
        }
        *GCD = r0;
        *x = s0;
        *y = t0;
    }

    /**
    * @brief Implementation using recursive algorithm
    * @param A
    * @param B
    * @param GCD
    * @param x
    * @param y
    */
    void extendedEuclid(uint64_t A, uint64_t B, uint64_t *GCD, int64_t *x, int64_t *y) {
        if (B > A)
            std::swap(A, B);  // Ensure that A >= B

        if (B == 0) {
            *GCD = A;
            *x = 1;
            *y = 0;
        } else {
            extendedEuclid(B, A % B, GCD, x, y);
            int64_t temp = *x;
            *x = *y;
            *y = temp - (A / B) * (*y);
        }
    }
}  // namespace math

/**
 * @brief Self-test implementations
 * @returns void
 */
static void tests(){
    assert(true);
}

/**
 * @brief Main function
 * @return 0 on exit
 */
int main() {
    tests();
    uint64_t a, b, gcd;
    int64_t x, y;
    std::cin >> a >> b;
    math::extendedEuclid(a, b, &gcd, &x, &y);
    std::cout << gcd << " " << x << " " << y << std::endl;
    math::extendedEuclid_1(a, b, &gcd, &x, &y);
    std::cout << gcd << " " << x << " " << y << std::endl;
    return 0;
}

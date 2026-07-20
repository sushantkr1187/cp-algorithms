---
tags:
    - Translated
e_maxx_link: bpsw
---

# Baillie–PSW Primality Test

## Introduction

The **Baillie–PSW (BPSW) primality test** is a probabilistic primality test named after its creators **Robert Baillie, Carl Pomerance, John Selfridge, and Samuel Wagstaff**. It was introduced in **1980** and remains one of the most reliable practical primality tests.

To date, **no composite number is known to pass the standard Baillie–PSW test**, although no proof has been found that such numbers do not exist.

The algorithm has been exhaustively verified for all integers up to $10^{15}$. In addition, extensive searches for counterexamples have been performed using the **PRIMO** primality proving software, which is based on elliptic curve primality proving (ECPP). After running for several years, no counterexample was found. Based on these experiments, Martin conjectured that there is no Baillie–PSW pseudoprime smaller than $10^{1\,000\,000}$.

On the other hand, Carl Pomerance presented a heuristic argument in 1984 suggesting that infinitely many Baillie–PSW pseudoprimes may exist. Thus, although no counterexample has ever been discovered, the correctness of the test remains an open mathematical question.

The asymptotic complexity of the Baillie–PSW test is

$$
O(\log^3 N)
$$

bit operations.

Compared to a single round of the Miller–Rabin test, the Baillie–PSW test is typically **3–7 times slower**, but provides a dramatically higher level of confidence in practice.

Because of its excellent reliability-to-performance ratio, the Baillie–PSW test is widely used in mathematical software and competitive programming libraries as a practical primality test.

## Algorithm Overview

The standard Baillie–PSW primality test consists of the following steps:

1. Perform a **base-2 Miller–Rabin primality test**.
2. Perform a **strong Lucas–Selfridge primality test**.
3. Report the number as **prime** only if **both** tests succeed.

In practice, an additional preliminary trial division by small prime numbers is often performed before either test. Although this slightly increases the running time for prime inputs, it significantly improves performance on composite numbers by eliminating trivial factors early.

Therefore, a practical implementation typically follows this procedure:

1. Trial division by small primes.
2. Base-2 Miller–Rabin test.
3. Strong Lucas–Selfridge test.
4. Return **prime** only if every step succeeds.

### Why does the algorithm work?

The Baillie–PSW test relies on two important observations.

1. Both the Miller–Rabin test and the Lucas–Selfridge test are **one-sided probabilistic tests**. Whenever either test declares a number composite, that result is always correct. The only possible error is incorrectly classifying a composite number as prime.

2. No composite number is currently known that simultaneously passes both the base-2 Miller–Rabin test and the strong Lucas–Selfridge test.

The second statement has never been proven. In fact, Carl Pomerance presented heuristic arguments suggesting that infinitely many such pseudoprimes may exist. Nevertheless, despite extensive computational searches, no counterexample has been discovered so far, making the Baillie–PSW test extremely reliable in practice.

## Implementation Notes

The implementations presented in this article are written in **C++** and are designed to be generic. Most functions are implemented as templates, allowing them to work with both built-in integer types and arbitrary-precision integer classes.

Only the core algorithms are shown throughout the article. Several auxiliary routines are omitted for brevity. Their interfaces are listed below, along with a short description of their purpose.

```cpp
//! Absolute value (64-bit overloads)
long long abs(long long n);
unsigned long long abs(unsigned long long n);

//! Returns true if n is even.
template <class T>
bool even(const T& n);

//! Divides n by 2.
template <class T>
void bisect(T& n);

//! Multiplies n by 2.
template <class T>
void redouble(T& n);

//! Returns true if n is a perfect square.
template <class T>
bool perfect_square(const T& n);

//! Computes ⌊√n⌋.
template <class T>
T sq_root(const T& n);

//! Returns the number of bits required to represent n.
template <class T>
unsigned bits_in_number(T n);

//! Returns the k-th bit of n (0-indexed).
template <class T>
bool test_bit(const T& n, unsigned k);

//! Computes a *= b (mod n).
template <class T>
void mulmod(T& a, T b, const T& n);

//! Computes a^k (mod n).
template <class T, class T2>
T powmod(T a, T2 k, const T& n);

//! Writes n = q · 2^p.
template <class T>
void transform_num(T n, T& p, T& q);

//! Greatest common divisor.
template <class T, class T2>
T gcd(const T& a, const T2& b);

//! Computes the Jacobi symbol J(a, b).
template <class T>
T jacobi(T a, T b);

//! Returns all primes up to b.
template <class T, class T2>
const std::vector<T>& get_primes(const T& b, T2& pi);

//! Trial division up to m.
//!
//! Returns:
//!   1  if n is proven prime,
//!   p  if p > 1 is a divisor of n,
//!   0  if the result is inconclusive.
template <class T, class T2>
T2 prime_div_trivial(const T& n, T2 m);
```

## Miller–Rabin Test

A detailed description of the Miller–Rabin primality test is available in the corresponding article on cp-algorithms. Here we only recall the implementation used as the first stage of the Baillie–PSW test.

The Miller–Rabin test performs a probabilistic primality check with complexity

$$
O(\log^3 N)
$$

bit operations.

In the Baillie–PSW test, **only a single Miller–Rabin round with base 2** is required. Although a single base is generally insufficient for a reliable primality test, combining it with the strong Lucas–Selfridge test yields the remarkably effective Baillie–PSW algorithm.

The implementation below performs the Miller–Rabin test for an arbitrary base.

```cpp
template <class T, class T2>
bool miller_rabin(T n, T2 b)
{
    // Handle trivial cases.
    if (n == 2)
        return true;
    if (n < 2 || even(n))
        return false;

    // Ensure gcd(n, b) = 1.
    // Otherwise, either a non-trivial divisor has been found,
    // or the base must be increased.
    if (b < 2)
        b = 2;
    for (T g; (g = gcd(n, b)) != 1; ++b)
        if (n > g)
            return false;

    // Write n − 1 = q · 2^p.
    T n_1 = n;
    --n_1;

    T p, q;
    transform_num(n_1, p, q);

    // Compute b^q mod n.
    T rem = powmod(T(b), q, n);

    if (rem == 1 || rem == n_1)
        return true;

    // Repeatedly square modulo n.
    for (T i = 1; i < p; ++i)
    {
        mulmod(rem, rem, n);
        if (rem == n_1)
            return true;
    }

    return false;
}
```

## Strong Lucas–Selfridge Test

The second stage of the Baillie–PSW test is the **strong Lucas–Selfridge primality test**. It combines two components:

1. **Selfridge's algorithm**, which selects suitable Lucas sequence parameters.
2. **The strong Lucas primality test**, executed using those parameters.

Unlike the Miller–Rabin test, which is based on modular exponentiation, this test relies on the arithmetic properties of **Lucas sequences**.

### Selfridge Parameter Selection

Selfridge's algorithm searches for the first integer

$$
D \in \{5,-7,9,-11,13,-15,\ldots\}
$$

such that

$$
\left(\frac{D}{N}\right) = -1
$$

and

$$
\gcd(D,N)=1,
$$

where

$$
\left(\frac{D}{N}\right)
$$

denotes the **Jacobi symbol**.

Once such a value of $D$ is found, the Lucas parameters are defined as

$$
P = 1,\qquad
Q = \frac{1-D}{4}.
$$

These parameters satisfy

$$
D=P^2-4Q.
$$

### Special Cases

The Selfridge parameter search assumes that

- $N$ is odd,
- $N>2$,
- $N$ is **not** a perfect square.

If any of these conditions fails, the number is immediately declared composite.

The perfect-square check is particularly important. If $N$ is a perfect square, the search eventually reaches a value of $D$ sharing a non-trivial divisor with $N$, causing

$$
\gcd(D,N)>1,
$$

which correctly identifies $N$ as composite.

Although the search for $D$ has no proven worst-case bound, it is extremely efficient in practice. Experimental evidence shows that only small values of $D$ are typically required, even for very large integers.

### Lucas Sequences

Given two integers \(P\) and \(Q\), the **Lucas sequences** \(U_k(P,Q)\) and \(V_k(P,Q)\) are defined recursively as

$$
U_0 = 0,\qquad U_1 = 1,
$$

$$
U_k = P U_{k-1} - Q U_{k-2}, \qquad k \ge 2,
$$

and

$$
V_0 = 2,\qquad V_1 = P,
$$

$$
V_k = P V_{k-1} - Q V_{k-2}, \qquad k \ge 2.
$$

For the Baillie–PSW test, the parameters \(P\) and \(Q\) are obtained using **Selfridge's parameter selection** described above.

Computing these sequences directly from the recurrence requires \(O(n)\) steps, which is impractical for large values of \(n\). Instead, the algorithm evaluates the required terms modulo \(N\) using identities analogous to binary exponentiation.

The following identities allow both sequences to be computed in \(O(\log N)\) time:

For even indices,

$$
U_{2k}=U_kV_k,
$$

$$
V_{2k}=V_k^2-2Q^k,
$$

while for odd indices,

$$
U_{2k+1}=\frac{PU_{2k}+V_{2k}}{2},
$$

$$
V_{2k+1}=\frac{DV_{2k}+PU_{2k}}{2},
$$

where

$$
D=P^2-4Q.
$$

These formulas make it possible to evaluate the required Lucas sequence terms using the binary representation of the index, in much the same way that modular exponentiation computes powers efficiently.

### Strong Lucas Primality Test

Let

$$
N+1=d\cdot2^s,
$$

where \(d\) is odd.

Using the parameters \(P\) and \(Q\) selected by Selfridge's algorithm, compute the Lucas sequence values modulo \(N\).

The number \(N\) passes the **strong Lucas test** if **either**

$$
U_d \equiv 0 \pmod N,
$$

or there exists an integer

$$
0 \le r < s
$$

such that

$$
V_{d2^r}\equiv0\pmod N.
$$

If neither condition holds, then \(N\) is composite.

The structure of this test closely resembles the strong Miller–Rabin test. Instead of repeatedly squaring modular powers, however, it repeatedly doubles the index of the Lucas sequences.

When combined with the base-2 Miller–Rabin test, this criterion forms the standard Baillie–PSW primality test.

### Selfridge Algorithm

The first step of the Lucas–Selfridge test is to determine suitable parameters for the Lucas sequences.

Selfridge's method searches the sequence

$$
5,\,-7,\,9,\,-11,\,13,\,-15,\,\ldots
$$

until it finds the first integer \(D\) satisfying

$$
\left(\frac{D}{N}\right) = -1
$$

and

$$
\gcd(D, N) = 1,
$$

where

$$
\left(\frac{D}{N}\right)
$$

denotes the **Jacobi symbol**.

Once such a value is found, the Lucas parameters are defined as

$$
P = 1,\qquad
Q = \frac{1-D}{4}.
$$

These parameters satisfy the identity

$$
D = P^2 - 4Q.
$$

#### Preconditions

Before running Selfridge's algorithm, the input number \(N\) must satisfy the following conditions:

- \(N\) is odd.
- \(N > 2\).
- \(N\) is **not** a perfect square.

If any of these conditions fails, the number is immediately declared composite.

The perfect-square check is essential. If \(N\) is a perfect square, the search for \(D\) will eventually encounter a value sharing a non-trivial divisor with \(N\), resulting in

$$
\gcd(D,N) > 1,
$$

which correctly identifies \(N\) as composite.

#### Practical Remarks

In theory, there is no proven upper bound on how large the required value of \(D\) may become.

In practice, however, the required value is remarkably small. For example,

- for all integers in the interval \([1,10^6]\), the maximum required value is \(47\);
- for integers in \([10^{19},\,10^{19}+10^6]\), the maximum observed value is \(67\).

Baillie and Wagstaff also provided theoretical evidence explaining this behavior in their original 1980 paper. Consequently, Selfridge's parameter selection is extremely efficient in practice.

### Lucas Sequences

Given integers \(P\) and \(Q\), define the Lucas sequences \(U_k\) and \(V_k\) by

$$
U_0 = 0,\qquad U_1 = 1,
$$

$$
U_k = P U_{k-1} - Q U_{k-2},
$$

and

$$
V_0 = 2,\qquad V_1 = P,
$$

$$
V_k = P V_{k-1} - Q V_{k-2}.
$$

For the parameters selected by Selfridge's algorithm,

$$
D=P^2-4Q.
$$

Let

$$
M = N - \left(\frac{D}{N}\right),
$$

where

$$
\left(\frac{D}{N}\right)
$$

is the Jacobi symbol.

If \(N\) is prime and

$$
\gcd(N,Q)=1,
$$

then

$$
U_M \equiv 0 \pmod N.
$$

For the Selfridge parameters,

$$
\left(\frac{D}{N}\right)=-1,
$$

so

$$
M=N+1,
$$

and therefore

$$
U_{N+1}\equiv0\pmod N.
$$

The converse is generally false: some composite numbers also satisfy this congruence. Such numbers are called **Lucas pseudoprimes**.

The Lucas primality test therefore computes \(U_M\) modulo \(N\). If the result is zero, the number is considered a probable prime.

### Efficient Computation of Lucas Sequences

A direct evaluation of the recurrence relations requires \(O(k)\) operations to compute \(U_k\) or \(V_k\), which is impractical for primality testing. Instead, we use identities that allow the required terms to be computed in \(O(\log k)\) time, similarly to binary exponentiation.

Suppose

$$
a,b
$$

are the distinct roots of the quadratic equation

$$
x^2-Px+Q=0.
$$

Then the Lucas sequences admit the closed-form expressions

$$
U_k=\frac{a^k-b^k}{a-b},
$$

and

$$
V_k=a^k+b^k.
$$

From these formulas, the following doubling identities can be derived:

$$
U_{2k}=U_kV_k \pmod N,
$$

$$
V_{2k}=V_k^2-2Q^k \pmod N.
$$

Now write

$$
M=E\cdot2^s,
$$

where \(E\) is odd.

Using the doubling identities repeatedly, we obtain

$$
U_M
=
U_E
V_E
V_{2E}
V_{4E}
\cdots
V_{2^{s-1}E}.
$$

Therefore,

$$
U_M\equiv0\pmod N
$$

if and only if at least one of the factors

$$
U_E,\;
V_E,\;
V_{2E},\;
V_{4E},\;
\dots,\;
V_{2^{s-1}E}
$$

is congruent to zero modulo \(N\).

Consequently, it is sufficient to compute only \(U_E\) and \(V_E\), after which all remaining terms can be generated efficiently using the doubling formulas.

# Computing \(U_E\) and \(V_E\)
To compute \(U_E\) and \(V_E\) efficiently, we use addition formulas analogous to those used in fast exponentiation.

For any non-negative integers \(i\) and \(j\),

$$
U_{i+j}
=
\frac{U_iV_j+U_jV_i}{2}
\pmod N,
$$

and

$$
V_{i+j}
=
\frac{V_iV_j+DU_iU_j}{2}
\pmod N,
$$

where

$$
D=P^2-4Q.
$$

The division by \(2\) is performed modulo \(N\). Since the Baillie–PSW test only considers odd integers \(N\), the modular inverse of \(2\) always exists.

These identities, together with the doubling formulas, allow the Lucas sequences to be evaluated using the binary representation of the index.

The computation proceeds as follows:

1. Initialize

   $$
   U=U_1,\qquad V=V_1.
   $$

2. Process the binary representation of \(E\) from the second most significant bit to the least significant bit.

3. For each bit:
   - Double the current indices using the doubling identities.
   - If the current bit is \(1\), apply the addition formulas to update the result.

This algorithm computes \(U_E\) and \(V_E\) in

$$
O(\log E)
$$

modular operations.

### Strong Lucas Test

After computing \(U_E\) and \(V_E\), the strong Lucas test proceeds exactly as the strong Miller–Rabin test does.

If either

$$
U_E\equiv0\pmod N
$$

or

$$
V_E\equiv0\pmod N,
$$

then \(N\) is declared a **probable prime**.

Otherwise, repeatedly compute

$$
V_{2E},
\;
V_{4E},
\;
\dots,
\;
V_{2^{s-1}E}
$$

using the doubling identity

$$
V_{2k}=V_k^2-2Q^k.
$$

If any of these values is congruent to zero modulo \(N\), then \(N\) passes the strong Lucas test.

If none of them is zero modulo \(N\), then \(N\) is composite.

This criterion is directly analogous to the strong Miller–Rabin test: instead of repeatedly squaring modular powers, we repeatedly double the indices of the Lucas sequences.

### Why Selfridge's Parameters?

The Lucas primality test depends on three parameters, \(D\), \(P\), and \(Q\), satisfying

$$
P > 0,
$$

and

$$
D=P^2-4Q\neq0.
$$

Different choices of these parameters lead to different variants of the Lucas test and, consequently, different sets of Lucas pseudoprimes.

A particularly important requirement is that \(D\) should **not** be a quadratic residue modulo \(N\). If \(D\) is a perfect square modulo \(N\), the Lucas test degenerates into a much weaker test.

Indeed, suppose

$$
D=b^2.
$$

Choosing

$$
P=b+2,\qquad
Q=b+1,
$$

gives

$$
U_{N-1}
=
\frac{Q^{\,N-1}-1}{Q-1},
$$

which behaves similarly to Fermat's primality test and therefore admits many pseudoprimes.

To avoid this situation, it is desirable to choose \(D\) such that

$$
\left(\frac{D}{N}\right)=-1,
$$

where

$$
\left(\frac{D}{N}\right)
$$

is the Jacobi symbol.

Selfridge's algorithm achieves this by selecting the first value of \(D\) from the sequence

$$
5,\,-7,\,9,\,-11,\,13,\,-15,\ldots
$$

for which

$$
\left(\frac{D}{N}\right)=-1,
$$

and then defining

$$
P=1,\qquad
Q=\frac{1-D}{4}.
$$

This is not the only possible choice of parameters. Other parameter-selection strategies have been proposed, leading to different Lucas tests and different sets of pseudoprimes.

However, extensive computational evidence suggests that Selfridge's parameters provide an excellent practical choice. Although Lucas pseudoprimes do exist, no composite number is currently known that passes **both** the base-2 Miller–Rabin test and the strong Lucas–Selfridge test, making the combined Baillie–PSW test exceptionally reliable in practice.

## Implementation of the Strong Lucas–Selfridge Test

The implementation below follows the algorithm described in the previous sections.

The procedure consists of four main stages:

1. Handle the trivial cases.
2. Compute the Selfridge parameters \(D\), \(P\), and \(Q\).
3. Evaluate the required Lucas sequence values using the binary method.
4. Perform the strong Lucas test.

The implementation assumes the availability of the helper functions introduced earlier, including modular multiplication, modular exponentiation, computation of the Jacobi symbol, and perfect-square detection.

```cpp
template <class T, class T2>
bool lucas_selfridge(const T& n, T2)
{
    // Handle trivial cases.
    if (n == 2)
        return true;
    if (n < 2 || even(n))
        return false;

    // Perfect squares are composite.
    if (perfect_square(n))
        return false;

    // Selfridge parameter selection.
    T2 dd;
    for (T2 d_abs = 5, d_sign = 1;; d_sign = -d_sign, d_abs += 2)
    {
        dd = d_abs * d_sign;

        T g = gcd(n, d_abs);
        if (1 < g && g < n)
            return false;

        if (jacobi(T(dd), n) == -1)
            break;
    }

    // Selfridge parameters.
    T2 p = 1;
    T2 q = (p * p - dd) / 4;

    // Write n + 1 = d · 2^s.
    T n1 = n;
    ++n1;

    T s, d;
    transform_num(n1, s, d);

    // Compute the Lucas sequences.
    ...
}
```

The implementation computes the values \(U_d\) and \(V_d\) using the binary representation of \(d\), similarly to binary exponentiation. If neither value is congruent to zero modulo \(N\), the remaining sequence values

$$
V_{2d},\;
V_{4d},\;
\dots,\;
V_{2^{s-1}d}
$$

are generated iteratively using the doubling identities.

If any of these values is congruent to zero modulo \(N\), the number passes the strong Lucas test. Otherwise, it is composite.

## Complete Baillie–PSW Test

With the Miller–Rabin test and the strong Lucas–Selfridge test implemented, the complete Baillie–PSW algorithm is straightforward.

A practical implementation consists of the following steps:

1. Perform trial division by a small set of prime numbers.
2. Run a base-2 Miller–Rabin test.
3. Run the strong Lucas–Selfridge test.
4. Report the number as **prime** only if all three stages succeed.

The implementation is shown below.

```cpp
template <class T>
bool baillie_pomerance_selfridge_wagstaff(T n)
{
    // Trial division by small primes (e.g. up to 29).
    int div = prime_div_trivial(n, 29);

    if (div == 1)
        return true;

    if (div > 1)
        return false;

    // Base-2 Miller–Rabin test.
    if (!miller_rabin(n, 2))
        return false;

    // Strong Lucas–Selfridge test.
    return lucas_selfridge(n, 0);
}
```

Although the implementation is remarkably short, it is one of the most reliable probabilistic primality tests used in practice. Despite decades of extensive testing, no composite number is currently known to pass both the base-2 Miller–Rabin test and the strong Lucas–Selfridge test.

## Notes

The implementation presented above is generic and works with any integer type supporting the required arithmetic operations.

For applications restricted to built-in integer types, the implementation can be simplified considerably by removing templates and replacing the helper routines with direct arithmetic operations. This reduces the amount of code without changing the underlying algorithm.

## Complexity

The Baillie–PSW primality test consists of three stages:

1. Trial division by a fixed set of small primes.
2. A single Miller–Rabin test with base \(2\).
3. A strong Lucas–Selfridge test.

The trial division requires constant time for a fixed bound.

Both the Miller–Rabin and Lucas–Selfridge tests require

$$
O(\log^3 N)
$$

bit operations.

Therefore, the overall asymptotic complexity of the Baillie–PSW test is

$$
O(\log^3 N)
$$

bit operations.

In practice, the Lucas–Selfridge test is several times slower than a single Miller–Rabin round, but the combined algorithm remains extremely fast while providing a very high degree of confidence.

## References

1. Robert Baillie and Samuel S. Wagstaff, *Lucas Pseudoprimes*, Mathematics of Computation, 35(152), 1980.

2. Carl Pomerance, *Are There Counterexamples to the Baillie–PSW Primality Test?*, Mathematics of Computation, 1984.

3. Daniel J. Bernstein, *Distinguishing Prime Numbers from Composite Numbers: The State of the Art in 2004*.

4. Eric W. Weisstein, *Baillie–PSW Primality Test*, MathWorld.

5. Eric W. Weisstein, *Strong Lucas Pseudoprime*, MathWorld.

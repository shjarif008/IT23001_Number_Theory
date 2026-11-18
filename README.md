## Table of Contents

1. [How to run](#1-how-to-run)
2. [Roadmap of the workshop](#2-roadmap-of-the-workshop)
3. [Part I, Number Theory](#3-part-i--number-theory)
   - [Day 0: Divisibility and modular arithmetic](#day-0--divisibility-remainders-and-modular-arithmetic)
   - [Day 1: GCD, Euclid, Diophantine equations, modular division, CRT, fast exponentiation, Fermat and Euler](#day-1--gcd-euclid-diophantine-equations-crt-and-modular-exponentiation)
   - [Day 2: Primes and primality testing](#day-2--primes-and-primality-testing)
4. [Part II, Cryptography](#4-part-ii--cryptography)
   - [Vocabulary and threat model](#41-vocabulary-and-the-alice-bob-eve-model)
   - [Classical ciphers: XOR, Caesar, Vigenère](#42-classical-ciphers)
   - [Perfect secrecy, one-time pad, many-time pad](#43-perfect-secrecy-the-one-time-pad-and-why-reuse-is-fatal)
   - [Computational security and brute force](#44-computational-security)
   - [Diffie-Hellman](#45-diffie-hellman-key-exchange)
   - [RSA](#46-rsa)
   - [Breaking badly-used RSA (Assignment Q1-Q7)](#47-breaking-badly-used-rsa-rsa-assignment-q1q7)
   - [Shor's algorithm](#48-shors-algorithm-the-classical-skeleton)
5. [Assignment answers at a glance](#5-assignment-answers-at-a-glance)
6. [Cheat sheet](#6-cheat-sheet)
7. [References](#7-references)

---

## Part I: Number Theory

### Day 0: Divisibility, Remainders and Modular Arithmetic

#### Definitions

**Divisibility.** $b \mid a$ iff there is an integer $k$ with $a = b\cdot k$.

**Division with remainder.** For $b>0$ there are unique integers $q, r$ with

$$a = q\,b + r,\qquad 0 \le r < b .$$

Python's `divmod(a, b)` and `a % b` follow exactly this convention for $b>0$ (even for negative $a$):
$-13 = (-5)\cdot 3 + 2$.

**Congruence.** $a \equiv b \pmod m$ means $a,b$ have the same remainder mod $m$, equivalently $m \mid (a-b)$.

#### Key lemmas (with one-line proofs)

| Lemma | Proof idea |
|---|---|
| $c\mid a,\ c\mid b \Rightarrow c \mid a\pm b$ | $a=ck_1,\ b=ck_2 \Rightarrow a\pm b = c(k_1\pm k_2)$ |
| $b \mid a \Rightarrow b \mid ac$ | $ac = b(kc)$ |
| $a\equiv b,\ c\equiv d \Rightarrow a+c \equiv b+d$ and $ac\equiv bd \pmod m$ | $ac\equiv ad \equiv bd$ (use the constant rule twice) |

Because of the last row we may **replace every number by its remainder** (or by $0,\pm1,\dots$) at any time.

#### Worked problems from the slides

```python
# 17*(12*19+5) - 23  (mod 3)  -> 2
print((17 * (12*19 + 5) - 23) % 3)        # 2
# 14+41+20+13+29 (mod 4) -> 1
print(sum([14, 41, 20, 13, 29]) % 4)      # 1
# last two digits of 99^99 -> 99   because 99 ≡ -1 (mod 100) and (-1)^99 = -1
print(pow(99, 99, 100))                   # 99
```

**Divisibility by 3.** Since $10\equiv 1 \pmod 3$, we have $10^k\equiv 1$, hence

$$3475 = 3\cdot10^3+4\cdot10^2+7\cdot 10+5 \equiv 3+4+7+5 = 19 \equiv 1 \pmod 3 .$$

So a number is congruent mod 3 (and mod 9) to its digit sum.

#### Operation tables and the "division problem"

Addition mod $m$ always has an inverse (subtraction), since each row of the addition table contains every
residue. **Multiplication does not**: look at the tables for $m=7$ and $m=6$.

![Multiplication tables mod 7 and mod 6](images/mult_tables.png)

* mod 7 (prime): every non-zero row is a permutation of $0..6$, so division by any $a\ne0$ works.
* mod 6: rows $2,3,4$ miss residues, so e.g. there is **no** $x$ with $3x\equiv 1 \pmod 6$.

This is the puzzle that Day 1 resolves: *division mod $n$ works exactly when $\gcd(a,n)=1$.*

---

### Day 1: GCD, Euclid, Diophantine Equations, CRT and Modular Exponentiation

#### 1.1 Greatest common divisor and Euclid's algorithm

$\gcd(a,b)$ is the largest integer dividing both. The naive method (try every $d$) needs about $\min(a,b)$
divisions: *more than a minute* for the 12-digit example in the slides, and thousands of years for
hundred-digit numbers.

**Euclid's lemma:** $d \mid a$ and $d\mid b \iff d \mid (a-b)$ and $d \mid b$. Repeating subtraction collapses to
taking remainders:

$$\gcd(a,b) = \gcd(b,\ a \bmod b),\qquad \gcd(a,0)=a .$$

```python
def gcd_euclid(a, b):
    while b:
        a, b = b, a % b
    return a

print(gcd_euclid(1849639579327, 790933790547))   # 3416723  (11 loop steps, as in the slides)
```

**Why it is fast.** For $a \ge b>0$: $(a \bmod b) < a/2$.

* If $b \le a/2$, then $a \bmod b < b \le a/2$.
* If $b > a/2$, then $a \bmod b = a-b < a/2$.

So the numbers at least halve every step and the loop runs at most $\log_2 a + \log_2 b$ times (about 660 steps
for two 100-digit numbers; fewer than 16 610 for 5 000-digit numbers).

![Euclid's algorithm iterations vs. the logarithmic bound](images/euclid_steps.png)

#### 1.2 Extended Euclid: Bézout's identity (Assignment Q1)

**Certificate lemma.** If $d\mid a$, $d \mid b$ and $d = ax+by$ for integers $x,y$, then $d=\gcd(a,b)$.
(Any common divisor divides $ax+by=d$, so $\gcd\le d$; and $d$ is itself a common divisor, so $d\le\gcd$.)

**Recursion.** If $d = bp + (a \bmod b)\,q$, then because $a \bmod b = a - \lfloor a/b\rfloor b$,

$$d = a\,q + b\,\Bigl(p - \Bigl\lfloor \tfrac{a}{b}\Bigr\rfloor q\Bigr).$$

```python
def extended_gcd(a, b):
    """Return (d, x, y) with d = gcd(a, b) = a*x + b*y.   Requires a >= b >= 0, a+b > 0."""
    assert a >= b >= 0 and a + b > 0
    if b == 0:
        d, x, y = a, 1, 0
    else:
        d, x1, y1 = extended_gcd(b, a % b)
        x, y = y1, x1 - (a // b) * y1
    assert a % d == 0 and b % d == 0 and d == a * x + b * y   # the template's assertions
    return d, x, y

print(extended_gcd(8, 4))      # (4, 0, 1)
print(extended_gcd(27, 15))    # (3, -1, 2)      3 = 27*(-1) + 15*2
print(extended_gcd(391, 299))  # (23, -3, 4)     matches the slide: 23 = 391*(-3) + 299*4
```

*Edge cases handled:* $b=0$ (base case), $a=b$, and the template's precondition $a \ge b$.

#### 1.3 LCM (Assignment Q2)

```python
def lcm(a, b):
    return a // extended_gcd(max(a, b), min(a, b))[0] * b    # divide first to keep numbers small

assert lcm(9, 3) == 9 and lcm(17, 37) == 629 and lcm(24, 16) == 48
```

#### 1.4 Linear Diophantine equations (Assignment Q3)

> **Theorem.** $ax+by=c$ has an integer solution $\iff \gcd(a,b)\mid c$.

*($\Rightarrow$)* $c=ax+by=d(px+qy)$. *($\Leftarrow$)* If $c=td$, multiply Bézout's $ax_0+by_0=d$ by $t$.

**All solutions.** With $d=\gcd(a,b)$, $a=dp$, $b=dq$ and one solution $(x_0,y_0)$:

$$x = x_0 + t\,q,\qquad y = y_0 - t\,p,\qquad t\in\mathbb Z .$$

(Proof: two solutions differ by $(\Delta x,\Delta y)$ with $p\Delta x = -q\Delta y$; since $\gcd(p,q)=1$,
Euclid's lemma forces $q \mid \Delta x$.)

```python
def diophantine(a, b, c):
    d, x, y = extended_gcd(max(a, b), min(a, b))
    if a < b:                      # extended_gcd needs a >= b, so swap back
        x, y = y, x
    assert c % d == 0
    k = c // d
    x, y = x * k, y * k
    assert a * x + b * y == c
    return x, y

print(diophantine(34, 51, 170))   # (-10, 10)
print(diophantine(91, 117, 65))   # (20, -15)
```

*Slide example:* $3x+5y=22$ has $(x_0,y_0)=(9,-1)$, so all solutions are $(9+5t,\,-1-3t)$.
The "apple costs 22 pesos" story asks for $x\ge0$ and $y \le 0$, which gives $t\ge 0$.

#### 1.5 Modular inverses and division (Assignment Q4)

> **Theorem.** $a$ has an inverse mod $n$ $\iff \gcd(a,n)=1$. The inverse is unique.

$ax\equiv1 \pmod n \iff ax+kn=1$ is a Diophantine equation, solvable iff $\gcd(a,n)\mid 1$.
Uniqueness: $x = x(ay) = (xa)y = y$.

Then $b/a \equiv b\cdot a^{-1} \pmod n$.

```python
def _egcd(a, b):                       # works for any a, b (no ordering requirement)
    if b == 0:
        return a, 1, 0
    g, x, y = _egcd(b, a % b)
    return g, y, x - (a // b) * y

def mod_inverse(a, n):
    g, s, _ = _egcd(a, n)
    assert g == 1, "inverse exists only when gcd(a, n) = 1"
    return s % n                       # normalise into 0..n-1

def divide(b, a, n):
    assert n > 1 and a > 0
    a_inv = mod_inverse(a, n)
    assert (a * a_inv) % n == 1
    return (a_inv * b) % n

print(divide(7, 2, 9))        # 8     (slide: 7/2 ≡ 8 mod 9, check 8*2 = 16 ≡ 7)
print(divide(223, 67, 54))    # 13
print(divide(9, 5, 3))        # 0     (since 9 ≡ 0 mod 3)
```

(Python 3.8+ also offers `pow(a, -1, n)`, which is useful for cross-checking.)

#### 1.6 Chinese Remainder Theorem (Assignment Q6)

> If $\gcd(a,b)=1$, then for every $(r_a,r_b)$ there is $n$ with $n\equiv r_a \pmod a$, $n\equiv r_b \pmod b$,
> and $n$ is unique modulo $ab$.

Existence by counting: the map $r\mapsto (r \bmod a,\ r\bmod b)$ from $\mathbb Z_{ab}$ to $\mathbb Z_a\times\mathbb Z_b$
is injective (if both residues agree, then $a\mid\Delta$ and $b\mid\Delta$, hence $ab\mid\Delta$), and both sets
have $ab$ elements, so it is a bijection.

![CRT grid](images/crt_grid.png)

**Constructive proof.** Take $ax+by=1$. Then $ax\equiv 1 \pmod b$ and $by\equiv1 \pmod a$, so

$$\boxed{\,n = r_a\,b\,y + r_b\,a\,x \pmod{ab}\,}$$

```python
def chinese_remainder(n1, r1, n2, r2):
    g, x, y = _egcd(n1, n2)          # n1*x + n2*y = 1
    assert g == 1
    return (r1 * n2 * y + r2 * n1 * x) % (n1 * n2)

print(chinese_remainder(13, 7, 6, 5))      # 59
print(chinese_remainder(28, 20, 19, 3))    # 440
```

*Consequences from the slides:* if $a\mid b$ then $n\bmod b$ determines $n \bmod a$, whereas for coprime moduli the
residues are **independent**. Large integers can be stored as residues modulo several big primes, which makes
addition and multiplication fast (this is used in practice).

#### 1.7 Fast modular exponentiation (Assignment Q5)

Computing $b^e \bmod m$ by $e$ multiplications is hopeless when $e$ has 1000 digits. Use the binary expansion
$e=\sum \epsilon_k 2^k$:

$$b^e = \prod_{k:\ \epsilon_k=1} b^{2^k},\qquad b^{2^{k+1}} = \bigl(b^{2^k}\bigr)^2 .$$

Squarings: $\log_2 e$; multiplications: at most $\log_2 e$ $\Rightarrow$ **at most $2\log_2 e$ in total**.

Worked example from the slides, $7^{13} \bmod 11$ with $13 = 1101_2 = 8+4+1$:

| $k$ | $7^{2^k} \bmod 11$ | used? |
|---|---|---|
| 0 | 7 | yes |
| 1 | 5 | no |
| 2 | 3 | yes |
| 3 | 9 | yes |

$7\cdot3\cdot9 = 189 \equiv 2 \pmod{11}$.

```python
def fast_pow(b, e, m):
    """b**e mod m with about 2*log2(e) modular multiplications (no built-in pow)."""
    result, b = 1 % m, b % m
    while e:
        if e & 1:                 # this bit of e is set -> multiply the current square in
            result = result * b % m
        b = b * b % m             # square
        e >>= 1
    return result

print(fast_pow(4, 5, 3))          # 1
print(fast_pow(7, 13, 11))        # 2
for _ in range(2000):             # randomised cross-check vs. the built-in
    b, e, m = random.randint(0, 10**6), random.randint(0, 10**6), random.randint(1, 10**6)
    assert fast_pow(b, e, m) == pow(b, e, m)
```

#### 1.8 Fermat's Little Theorem, Euler's totient and Euler's theorem

> **Fermat.** If $p$ is prime and $p\nmid a$, then $a^{p-1}\equiv 1 \pmod p$.

*Proof sketch (cycle argument from the slides).* Multiplication by $a$ permutes the $p-1$ non-zero residues,
so the graph $r\to ar$ splits into disjoint cycles, all of the same length $\ell$ (because
$a^\ell r\equiv r$ means $a^\ell\equiv1$ for every $r$). With $c$ cycles, $c\,\ell=p-1$ and
$a^{p-1}=(a^\ell)^c\equiv1$.

**Corollary (speed-up).** $a^n \equiv a^{\,n \bmod (p-1)} \pmod p$.

**Euler's totient** $\varphi(n)$ = number of $k\in[0,n-1]$ with $\gcd(k,n)=1$.

* $\varphi(p) = p-1$
* $\varphi(pq) = (p-1)(q-1)$ for distinct primes. (CRT: $\gcd(r,pq)=1 \iff r_p\neq0$ and $r_q\ne0$, giving $(p-1)(q-1)$ pairs.)

![Totient function](images/totient.png)

> **Euler.** If $\gcd(a,n)=1$, then $a^{\varphi(n)}\equiv 1 \pmod n$. (Same cycle proof, applied to the $\varphi(n)$ units.)

```python
def phi(n):                                   # slow but transparent
    return sum(1 for k in range(n) if gcd(k, n) == 1) if n > 1 else 1

print(phi(10), phi(55))                       # 4 40
assert all(pow(a, phi(n), n) == 1 for n in range(2, 300) for a in range(1, n) if gcd(a, n) == 1)
```

*Why this matters:* $\varphi(n)$ is easy to compute **if you know the factorisation** of $n$, and no fast
method is known otherwise. That asymmetry is exactly the **trapdoor** used by RSA.

---

### Day 2: Primes and Primality Testing

#### 2.1 How many primes? (infinitely many, in many ways)

* **Euclid:** if $p_1,\dots,p_k$ were all primes, any prime divisor of $\prod p_i+1$ is new. (Also try $\prod p_i-1$.)
* **Factorial:** $n!+1$ has a prime factor $>n$.
* **Analytic:** the harmonic series $\sum 1/n$ diverges, but only finitely many primes would make $\prod(1-1/p)^{-1}$ finite.
* Others in the slides: Goldbach's proof, Saidak's proof, Furstenberg's topological proof.
* **Special forms.** $N=(2p_1\cdots p_k)^2+1$ gives infinitely many primes $\equiv1 \pmod 4$;
  $N=4p_1\cdots p_k-1$ gives infinitely many $\equiv 3 \pmod 4$.
* **Dirichlet:** $ak+b$ contains infinitely many primes iff $\gcd(a,b)=1$ (needs analytic number theory).

**Prime Number Theorem.** With $\pi(x)$ the number of primes $\le x$,

$$\pi(x)\sim \frac{x}{\ln x}.$$

```python
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

P = sieve(100_000)
print(len(P), int(1e5 / math.log(1e5)))     # 9592  vs  8685
```

![Prime number theorem](images/pnt.png)

![Ulam spiral](images/ulam.png)

#### 2.2 Generating and finding primes

* No non-constant polynomial in one variable takes only prime values ($n^2+n+41$ is prime for $n=0..39$, then fails at $n=40$).
  Matiyasevich: a polynomial in 10 variables suffices.
* **Sieve of Eratosthenes** (above): cross out multiples of every prime $\le\sqrt n$.
* **Trial division** up to $\sqrt n$, optimised by *wheeling*: the wheel of 30 only tests residues
  $1,7,11,13,17,19,23,29$.

#### 2.3 Pseudoprimes: Fermat test and its failures

*If $n$ is prime then statement $S$ holds* tells us only that failure $\Rightarrow$ composite. The converse can fail.

**Fermat test:** $a^{n-1}\equiv 1 \pmod n$ for $\gcd(a,n)=1$.

```python
Pset = set(P)
print([n for n in range(3, 1000) if n not in Pset and pow(2, n - 1, n) == 1])   # [341, 561, 645]
print(pow(3, 340, 341), pow(2, 90, 91), pow(3, 90, 91))                         # 56 64 1
```

* $341$ is a pseudoprime base 2 but **not** base 3 ($3^{340}\equiv56$).
* $91$ is a pseudoprime base 3 but **not** base 2 ($2^{90}\equiv 64$).

Combining bases helps, **but**:

> **Carmichael numbers** are composite $n$ with $a^n\equiv a \pmod n$ for *every* $a$.
> Korselt's criterion: $n$ is Carmichael $\iff$ $n$ is composite, squarefree, and $(p-1)\mid(n-1)$ for every prime $p\mid n$.

```python
def is_carmichael(n):
    return n not in Pset and n > 2 and all(pow(a, n - 1, n) == 1 for a in range(2, n) if gcd(a, n) == 1)
print([n for n in range(3, 3000) if is_carmichael(n)])    # [561, 1105, 1729, 2465, 2821]
```

![Fermat liars](images/fermat_liars.png)

The three spikes at fraction $1.0$ are the Carmichael numbers $561$, $1105$ and $1729$. For these, **no** base
exposes compositeness, which is why the Fermat test alone is not enough.

#### 2.4 Miller–Rabin: the standard fix

Write $n-1=2^s d$ with $d$ odd. For a base $a$, $n$ passes if $a^d\equiv 1$ or $a^{2^r d}\equiv -1$ for some $0\le r<s$.
This uses the slide's fact: *$p$ is prime iff the only square roots of $1$ mod $p$ are $\pm1$*. Carmichael numbers
no longer fool it.

```python
def miller_rabin(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
    if n < 2: return False
    for p in bases:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

assert all(miller_rabin(n) == (n in Pset) for n in range(2, 5000))
assert not miller_rabin(561) and not miller_rabin(1729)           # Carmichael numbers caught
print(miller_rabin(2**127 - 1))                                    # True
```

#### 2.5 Lucas' theorem, Lucas–Lehmer, Mersenne primes and perfect numbers

**Lucas' theorem.** If for every prime $q\mid n-1$ there is $a$ with $a^{n-1}\equiv1$ and $a^{(n-1)/q}\not\equiv 1 \pmod n$, then $n$ is prime.
(The $n-1$ test, Pépin's test and Proth's theorem are refinements.)

**Mersenne primes** $M_p=2^p-1$. Mersenne's conjecture "all $2^n-1$ are prime" is false. If $2^n-1$ is prime then $n$ is prime.

**Lucas–Lehmer.** Let $s_1=4$, $s_k=s_{k-1}^2-2$. For odd prime $p$: $M_p$ is prime $\iff M_p \mid s_{p-1}$.

```python
def lucas_lehmer(p):
    if p == 2: return True
    s, M = 4, 2**p - 1
    for _ in range(p - 2):
        s = (s * s - 2) % M
    return s == 0

print([p for p in P if p < 130 and lucas_lehmer(p)])
# [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]   exactly the exponents on the slide
print([2**(p-1) * (2**p - 1) for p in P if p < 14 and lucas_lehmer(p)])
# [6, 28, 496, 8128, 33550336]                      perfect numbers
```

**Antiprimes (highly composite numbers):** if $n=\prod p_i^{a_i}$ then the number of divisors is $d(n)=\prod(a_i+1)$.

#### 2.6 Day 2 assignment: proofs and checks

**Q1: repunits.** Let $a(n)=\underbrace{11\cdots1}_{n}=\frac{10^n-1}{9}$ and $A_p=\{n : p \mid a(n)\}$.

$p\mid a(n)\iff 9p \mid 10^n-1$. If $p\in\{2,5\}$ then $10^n-1$ is coprime to $p$ (it is odd and
$\equiv-1 \pmod{5}$), so $A_p=\emptyset$. Otherwise $\gcd(10,9p)=1$, so by Euler $10^{\varphi(9p)}\equiv1$ and
$A_p$ is non-empty and infinite. **Bonus:**

```python
def repunit_mod(p, n): return ((pow(10, n, 9*p) - 1) // 9) % p
for p in [3, 7, 11, 13, 37]:
    print(p, [n for n in range(1, 60) if repunit_mod(p, n) == 0][:6])
# 3  -> 3, 6, 9, ...   7 -> 6, 12, ...   11 -> 2, 4, ...   13 -> 6, 12, ...   37 -> 3, 6, ...
```

**Q2: even perfect numbers** (Euclid–Euler). $k$ is an even perfect number $\iff k=2^{n-1}(2^n-1)$ with $2^n-1$ prime.

* ($\Leftarrow$) With $q=2^n-1$ prime, $\sigma(k)=\sigma(2^{n-1})\sigma(q)=(2^n-1)(q+1)=(2^n-1)2^n=2k$. ✔
* ($\Rightarrow$) Write $k=2^{m}u$, $u$ odd, $m\ge1$. Perfect means $\sigma(k)=2k$, i.e. $(2^{m+1}-1)\sigma(u)=2^{m+1}u$.
  Since $\gcd(2^{m+1}-1,\,2^{m+1})=1$, we get $u=(2^{m+1}-1)t$ and $\sigma(u)=2^{m+1}t$. If $t>1$ then $u$ has at least the divisors
  $1,\ t,\ u$ so $\sigma(u)\ge u+t+1 > u+t=2^{m+1}t$, contradiction. So $t=1$, $u=2^{m+1}-1$ and $\sigma(u)=u+1$ forces $u$ prime. ∎

**Q3: infinitely many Fermat pseudoprimes to base $a\ge 2$.** For an odd prime $p\nmid a(a^2-1)$, put

$$n=\frac{a^{2p}-1}{a^{2}-1}=1+a^2+a^4+\dots+a^{2(p-1)} .$$

$n$ is composite since it equals $\frac{a^p-1}{a-1}\cdot\frac{a^p+1}{a+1}$ with both factors $>1$.
We have $a^{2p}\equiv1 \pmod n$, It remains to show $2p\mid n-1$:
$n-1=a^2+\dots+a^{2p-2}$ has $p-1$ (an even number of) terms, so it is even; and by Fermat
$p \mid a^{2p-2}-1=(a^2-1)\cdot\frac{a^{2p-2}-1}{a^2-1}$ with $p\nmid a^2-1$, so $p\mid \frac{a^{2p-2}-1}{a^2-1}$ and
$n-1=a^2\cdot\frac{a^{2p-2}-1}{a^2-1}$ is divisible by $p$. Hence $a^{n-1}\equiv1\pmod n$, and different $p$ give different $n$. ∎

---

## 4. Part II: Cryptography

### 4.1 Vocabulary and the Alice-Bob-Eve model

| Term | Meaning |
|---|---|
| plaintext $p$ / ciphertext $c$ | original message / what travels over the channel |
| cipher | the encryption algorithm |
| key | the secret that Alice and Bob share (symmetric) or hold in pairs (asymmetric) |
| cryptanalysis | breaking a cipher without the key |
| **Kerckhoffs's principle** | assume Eve knows *everything except the key* |
| security by obscurity | secrecy of the algorithm itself, which is fragile and discouraged |

**Handout brain-teasers.**

* *"Since Everyone Can Read, Encoding Text In Neutral Sentences Is Doubtfully Effective"* is an acrostic: the first letters spell **SECRET INSIDE**.
* Entropy of 5 fair-die throws: $5\log_2 6\approx 12.92$ bits.
* Entropy of the weather if rain has probability 1: $0$ bits (a deterministic outcome carries no information).
* Average brute-force effort for an $n$-bit key: $2^{n-1}$ trials.

### 4.2 Classical ciphers

#### XOR

$A\oplus0=A,\quad A\oplus A=0,\quad (A\oplus B)\oplus C=A\oplus(B\oplus C)$, so

$$c=m\oplus k\ \Longrightarrow\ c\oplus k=m\oplus k\oplus k=m .$$

```python
print(format(0b011011 ^ 0b010101, "06b"))    # 001110   (handout exercise)
```

#### Caesar cipher and frequency analysis

Shift every letter by $k$: $c\equiv p+k \pmod{26}$. Only 26 keys, and worse: English has a stable letter
"fingerprint" ($e,t,a,o,\dots$), which Al-Kindi used around 800 years after Caesar. A shift-3 ciphertext makes
**h** the most frequent letter instead of **e**.

```python
# approximate English letter frequencies (%), the "fingerprint"
ENG = dict(zip("etaoinshrdlcumwfgypbvkjxqz",
               [12.7, 9.1, 8.2, 7.5, 7.0, 6.7, 6.3, 6.1, 6.0, 4.3, 4.0, 2.8, 2.8, 2.4, 2.4,
                2.2, 2.0, 2.0, 1.9, 1.5, 1.0, 0.8, 0.15, 0.15, 0.1, 0.07]))

def caesar(t, k):
    return "".join(chr((ord(c) - 97 + k) % 26 + 97) if c.isalpha() else c for c in t.lower())

def crack_caesar(ct):                       # score every shift against English letter frequencies
    best = max(range(26),
               key=lambda k: sum(ENG.get(c, 0) for c in caesar(ct, -k) if c.isalpha()))
    return best, caesar(ct, -best)
# -> recovers the shift 3 on a 200-letter English sample
```

#### Vigenère cipher: a flatter fingerprint

Repeat a keyword along the message and shift each letter by the matching key letter (A=1, B=2, …, as in the handout):

| | | | | | | | | | | | | | |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Plaintext | T | H | I | S | T | H | E | M | E | S | S | A | G | E |
| Key | A | B | C | D | A | B | C | D | A | B | C | D | A | B |
| Ciphertext | U | J | L | W | U | J | H | Q | F | U | V | E | H | G |

```python
def vigenere(t, key, dec=False, off=0):
    out, i = [], 0
    for c in t.lower():
        if c.isalpha():
            s = ord(key.lower()[i % len(key)]) - 97 + off
            out.append(chr((ord(c) - 97 + (-s if dec else s)) % 26 + 97)); i += 1
        else:
            out.append(c)
    return "".join(out)

print(vigenere("THISTHEMESSAGE", "ABCD", off=1).upper())   # UJLWUJHQFUVEHG   (matches the handout)
```

![Letter-frequency fingerprints](images/freq.png)

The fingerprint is flatter but a *partial* fingerprint survives because the key repeats. Eve first finds the key
length $L$: the average **index of coincidence**

$$\mathrm{IC}=\frac{\sum_i f_i(f_i-1)}{N(N-1)}$$

of every $L$-th letter jumps to about $0.067$ (English) at the true length (and its multiples), versus about $0.038$ for random text.
Then it is just $L$ Caesar ciphers.

![Index of coincidence reveals key length](images/vigenere_ic.png)

### 4.3 Perfect secrecy, the one-time pad, and why reuse is fatal

If the key is **truly random, as long as the message, and never reused**, the ciphertext carries no information about the
plaintext: Shannon's *perfect secrecy*. The price: the key must be shared in advance and is as long as the message.

**Many-Time Pad.** Reuse the key once and Eve gets

$$c_1\oplus c_2=(m_1\oplus k)\oplus(m_2\oplus k)=m_1\oplus m_2 ,$$

which no longer depends on the key. Now she can **guess a crib** (e.g. `" the "`) in one message and check
whether the XOR reveals letters in the other.

```python
def xor_bytes(a, b): return bytes(x ^ y for x, y in zip(a, b))

m1, m2, key = b"steal the secret", b"the cat the bird", b"supersecretverys"
c1, c2 = xor_bytes(m1, key), xor_bytes(m2, key)
assert xor_bytes(c1, c2) == xor_bytes(m1, m2)            # the key cancels out

x12 = xor_bytes(c1, c2)
def try_crib(crib):
    for pos in range(len(x12) - len(crib) + 1):
        guess = xor_bytes(crib, x12[pos:pos + len(crib)])
        if all(chr(g).isalpha() or chr(g) == " " for g in guess):
            print(pos, crib, "->", guess)

try_crib(b" the ")      # 5 -> b'at th'   and   7 -> b'he se'
try_crib(b"at the ")    # 5 -> b' the se'
```

`" the "` at position 7 of one message and `"at th"` at position 5 of the other: extending the guess to
`"at the "` reveals `" the se"`, i.e. the start of "*steal* the secret" and "*the cat* the bird". With more
ciphertexts under the same key the whole plaintext falls out.

### 4.4 Computational security

Real systems accept *computational* security: brute force is possible in principle but takes absurdly long,
**assuming** the underlying maths problem is hard. At $10^{18}$ decryptions/s (a deliberate over-estimate):

| Cipher | Avg. keys tried $2^{n-1}$ | Time (s) | Time (years) |
|---|---|---|---|
| 56-bit DES | $3.6\times10^{16}$ | $3.6\times10^{-2}$ | $1.1\times10^{-9}$ |
| 128-bit AES | $1.7\times10^{38}$ | $1.7\times10^{20}$ | $5.4\times10^{12}$ |
| 168-bit 3DES | $1.9\times10^{50}$ | $1.9\times10^{32}$ | $5.9\times10^{24}$ |

**One-way functions** are easy forwards, hard backwards; a **trapdoor** (secret extra info) makes the reverse easy.
Two examples used below:

![Factoring vs. modular exponentiation](images/one_way.png)

* multiplying primes is easy, **factoring** is hard (RSA's trapdoor = knowing $p,q$);
* $g^a \bmod N$ is easy, the **discrete logarithm** is hard (Diffie–Hellman).

*Side-channel attacks* (Paul Kocher) skip the maths: timing and power consumption leak key bits, since inputs with
more 1-bits take more time/power.

### 4.5 Diffie-Hellman key exchange

**Cyclic groups.** $(\mathbb Z/N)^\*$ is cyclic (has a *primitive root* / generator $g$ whose powers hit every unit) exactly for
$N=2,\ 4,\ p^k,\ 2p^k$ with $p$ an odd prime. Mod 5: $2$ and $3$ are generators (period 4); $4$ has period 2.

![Powers modulo 11](images/dh_powers.png)

**Protocol** (public: $N,g$):

| Step | Alice | Eve sees | Bob |
|---|---|---|---|
| 1 | picks secret $a$ | $N=11,\ g=7$ | picks secret $b$ |
| 2 | sends $A=g^a \bmod N$ | $A=5$ | sends $B=g^b \bmod N$ |
| 3 | computes $K=B^a \bmod N$ | $B=3$ | computes $K=A^b \bmod N$ |

$$K=(g^b)^a=(g^a)^b=g^{ab} \pmod N .$$

```python
N, g, a, b = 11, 7, 2, 4
A, B = pow(g, a, N), pow(g, b, N)          # 5, 3 (sent in the clear)
k1, k2 = pow(B, a, N), pow(A, b, N)
assert k1 == k2 == 9
print("shared key:", k1)                   # 9
```

Eve must solve the **discrete logarithm problem** (recover $a$ from $g^a \bmod N$), for which no fast classical algorithm is known.
Note that DH solves the *key-distribution* problem that made symmetric ciphers impractical.

### 4.6 RSA

**Key generation**

1. Pick large primes $p\neq q$, set $N=pq$ and $\varphi(N)=(p-1)(q-1)$.
2. Choose $e$ with $\gcd(e,\varphi(N))=1$.
3. Compute $d\equiv e^{-1} \pmod{\varphi(N)}$ using extended Euclid.
4. **Public key** $(e,N)$; **private key** $d$ (or $p,q$).

**Encrypt / decrypt:**

$$c=M^e \bmod N,\qquad M=c^d \bmod N .$$

**Why it works.** $ed=1+k\varphi(N)$, so for $\gcd(M,N)=1$

$$c^d=M^{ed}=M\cdot\bigl(M^{\varphi(N)}\bigr)^k\equiv M\cdot1^k=M \pmod N$$

by Euler's theorem. (It also holds when $\gcd(M,N)>1$, via CRT with Fermat mod each prime.)

```python
# Toy example from the Demo Session: p=5, q=11, e=7
p, q, e = 5, 11, 7
N, phiN = p * q, (p - 1) * (q - 1)       # 55, 40
d = mod_inverse(e, phiN)                 # 23
c = pow(2, e, N)                         # 18
assert pow(c, d, N) == 2
print(N, phiN, d, c)                     # 55 40 23 18
```

Text messages are turned into integers in base 256 (as in the assignment):

```python
def to_int(s):
    r = 0
    for ch in s: r = r * 256 + ord(ch)
    return r

def to_str(n):
    r = ""
    while n > 0: r += chr(n % 256); n //= 256
    return r[::-1]

def Encrypt(message, modulo, exponent):               # RSA Assignment Q1
    return fast_pow(to_int(message), exponent, modulo)

def Decrypt(ciphertext, p, q, exponent):              # RSA Assignment Q2
    d = mod_inverse(exponent, (p - 1) * (q - 1))
    return to_str(fast_pow(ciphertext, d, p * q))

p, q, e = 1000000007, 1000000009, 23917
ct = Encrypt("attack", p * q, e)
print(ct, Decrypt(ct, p, q, e))                       # 577821702011958447 attack
```

### 4.7 Breaking badly-used RSA (RSA Assignment Q1-Q7)

RSA is secure **only if used correctly**. Each question below is a classic implementation mistake.

#### Q3: tiny message space → just encrypt every candidate

RSA (as used here) is deterministic, so if Eve knows the message is one of a few options she encrypts each and compares.
This is why real systems add random padding (OAEP).

```python
def DecipherSimple(ciphertext, modulo, exponent, potential_messages):
    for m in potential_messages:
        if Encrypt(m, modulo, exponent) == ciphertext:
            return m
    return "don't know"

mod3, exp3 = 1110014195838866450995043, 767549
c = Encrypt("wait", mod3, exp3)
print(DecipherSimple(c, mod3, exp3, ["attack", "don't attack", "wait"]))    # wait
```

#### Q4: one prime is small → trial division

```python
def DecipherSmallPrime(ciphertext, modulo, exponent):
    for p in range(2, 10**6):
        if modulo % p == 0:
            return Decrypt(ciphertext, p, modulo // p, exponent)
    return "don't know"
```

A $10^6$-step loop is instant, even when the other prime has hundreds of digits.

#### Q5: primes close together → Fermat factorisation

If $|p-q|$ is small then $N=pq=a^2-b^2$ with $a=\frac{p+q}{2}$ just above $\sqrt N$ and $b=\frac{q-p}{2}$ small, so
$$N=a^2-b^2=(a-b)(a+b).$$
Start at $a=\lceil\sqrt N\rceil$ and increase until $a^2-N$ is a perfect square.

```python
def fermat_factor(n):
    a = isqrt(n) + (0 if isqrt(n) ** 2 == n else 1)
    while True:
        b2 = a * a - n
        b = isqrt(b2)
        if b * b == b2:
            return a - b, a + b
        a += 1

def DecipherSmallDiff(ciphertext, modulo, exponent):
    p, q = fermat_factor(modulo)
    return Decrypt(ciphertext, p, q, exponent)

n = 1000000007 * 1000000009
print(fermat_factor(n))                                    # (1000000007, 1000000009)
print(DecipherSmallDiff(Encrypt("attack", n, 239), n, 239))  # attack
```

#### Q6: two moduli share a prime → one `gcd` breaks both

Insufficient randomness made two keys share $p$. Then $\gcd(N_1,N_2)=p$ (fast thanks to Euclid), and both factorisations follow.

```python
def DecipherCommonDivisor(c1, n1, e1, c2, n2, e2):
    p = gcd(n1, n2)
    return Decrypt(c1, p, n1 // p, e1), Decrypt(c2, p, n2 // p, e2)

pc, qa, qb = 1000000007, 1000000009, 998244353
n1, n2 = pc * qa, pc * qb
print(DecipherCommonDivisor(Encrypt("attack", n1, 239), n1, 239,
                            Encrypt("wait", n2, 65537), n2, 65537))   # ('attack', 'wait')
```

(My test uses different primes and exponents than the notebook's example. Each public exponent must be coprime to
$\varphi$ of its own modulus, otherwise the inverse $d$ does not exist.)

#### Q7: Håstad's broadcast attack ($e=2$)

Bob sends the same $m$ to two people: $c_1\equiv m^2 \pmod{n_1}$ and $c_2\equiv m^2 \pmod{n_2}$. By **CRT** (coprime moduli),
we reconstruct $m^2 \bmod n_1n_2$. If $m^2<n_1n_2$ (true for a short message) this *is* the integer $m^2$, so an
integer square root recovers $m$ with no factoring at all.

```python
def DecipherHastad(c1, n1, c2, n2):
    # CRT wants n1 > n2 (as in the Day 1 template)
    if n1 < n2: c1, n1, c2, n2 = c2, n2, c1, n1
    m_squared = chinese_remainder(n1, c1, n2, c2)
    return to_str(isqrt(m_squared))

p1, q1 = 790383132652258876190399065097, 662503581792812531719955475509
p2, q2 = 656917682542437675078478868539, 1263581691331332127259083713503
n1, n2 = p1 * q1, p2 * q2
print(DecipherHastad(Encrypt("attack", n1, 2), n1, Encrypt("attack", n2, 2), n2))   # attack
```

**Takeaways:** use large random primes far apart, never share primes, never use tiny exponents without padding,
and never encrypt guessable messages deterministically.

### 4.8 Shor's algorithm: the classical skeleton

Quantum computers threaten RSA because they can find the **period** $r$ of $a^x \bmod N$ quickly (that is the only
quantum step). Everything else is the number theory from above:

1. Pick $a<N$ with $\gcd(a,N)=1$.
2. Find the period $r$: smallest $r$ with $a^r\equiv1 \pmod N$. *(quantum step)*
3. If $r$ is even and $a^{r/2}\not\equiv -1\pmod N$, then
$$(a^{r/2}-1)(a^{r/2}+1)=kN=kpq .$$
4. Neither factor is divisible by $N$ (by minimality of $r$ and the assumption), so $p,q$ split between them:
$$p=\gcd(a^{r/2}-1,\,N),\qquad q=\gcd(a^{r/2}+1,\,N).$$

At least 50% of random choices of $a$ work.

```python
def order(x, n): return next(k for k in range(1, n + 1) if pow(x, k, n) == 1)   # the "quantum" step, brute-forced

def shor_classical(N, a):
    assert gcd(a, N) == 1
    r = order(a, N)
    if r % 2 or pow(a, r // 2, N) == N - 1:
        return None                                     # bad luck: pick another a
    return r, gcd(pow(a, r // 2) - 1, N), gcd(pow(a, r // 2) + 1, N)

print(shor_classical(35, 8))     # (4, 7, 5)   ->  8^4 - 1 = 4095 = 35*117 ;  gcd(63,35)=7 , gcd(65,35)=5
```
---

## 7. References

* Workshop repository: <https://github.com/mnp-club/Number_Theory_and_Cryptography_Workshop>
  (slides Day 0-2, assignments, *Cryptography Handout*, *Many-Time Pad*, *Demo Session*, *RSA Assignment*, *RSA Quest*)
* Maths and Physics Club, IIT Bombay: <http://mnp-club.github.io/>
* Reading list from the repo's `Cryptography/References`: Kerckhoffs's principle, one-time pad, Shannon's *Communication Theory of Secrecy Systems*,
  Kocher's timing attacks / physical key extraction, *Algorithms to Live By* (B. Christian), [interactive cipher playground](https://crypto.interactive-maths.com/)
* Course mentioned in the README: MA 523 (Basic Number Theory), Dept. of Mathematics, IIT Bombay

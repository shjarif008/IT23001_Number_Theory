import math, random, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import gcd, isqrt

IMG = "images/"

# ---------- Day 1 ----------
def extended_gcd(a, b):
    assert a >= b >= 0 and a + b > 0
    if b == 0:
        d, x, y = a, 1, 0
    else:
        d, x1, y1 = extended_gcd(b, a % b)
        x, y = y1, x1 - (a // b) * y1
    assert a % d == 0 and b % d == 0 and d == a*x + b*y
    return d, x, y

def lcm(a, b):
    return a // extended_gcd(max(a,b), min(a,b))[0] * b

def diophantine(a, b, c):
    d, x, y = extended_gcd(max(a,b), min(a,b))
    if a < b: x, y = y, x
    assert c % d == 0
    k = c // d
    x, y = x*k, y*k
    assert a*x + b*y == c
    return x, y

def mod_inverse(a, n):
    g, s, t = _egcd(a, n)
    assert g == 1
    return s % n

def _egcd(a, b):
    if b == 0: return a, 1, 0
    g, x, y = _egcd(b, a % b)
    return g, y, x - (a//b)*y

def divide(b, a, n): return (mod_inverse(a, n) * b) % n

def fast_pow(b, e, m):
    r, b = 1 % m, b % m
    while e:
        if e & 1: r = r*b % m
        b = b*b % m
        e >>= 1
    return r

def crt(n1, r1, n2, r2):
    g, x, y = _egcd(n1, n2); assert g == 1
    return (r1*n2*y + r2*n1*x) % (n1*n2)

print("egcd(391,299) =", extended_gcd(391,299))
print("egcd(27,15)   =", extended_gcd(27,15))
print("lcm(17,37)=", lcm(17,37), " lcm(24,16)=", lcm(24,16))
print("diophantine(34,51,170) =", diophantine(34,51,170))
print("diophantine(91,117,65) =", diophantine(91,117,65))
print("divide(7,2,9) =", divide(7,2,9), " divide(223,67,54) =", divide(223,67,54))
print("fast_pow(7,13,11) =", fast_pow(7,13,11), pow(7,13,11))
print("crt(13,7,6,5) =", crt(13,7,6,5), " crt(28,20,19,3) =", crt(28,20,19,3))
assert crt(13,7,6,5)==59 and crt(28,20,19,3)==440 and lcm(17,37)==629
assert fast_pow(4,5,3)==pow(4,5,3)
for _ in range(2000):
    b,e,m=random.randint(0,10**6),random.randint(0,10**6),random.randint(1,10**6)
    assert fast_pow(b,e,m)==pow(b,e,m)

# gcd step counts
def gcd_steps(a,b):
    s=0
    while b: a,b=b,a%b; s+=1
    return s
print("steps for the slide example:", gcd_steps(1849639579327, 790933790547))

# Fig 1: Euclid worst case (Fibonacci) vs bound
xs = list(range(2, 4000, 7))
fig, ax = plt.subplots(figsize=(6,3.6))
steps = [max(gcd_steps(n,k) for k in range(1,n,max(1,n//60))) for n in xs]
ax.plot(xs, steps, label="max steps observed")
ax.plot(xs, [math.log2(n)*2 for n in xs], "--", label=r"bound $\log_2 a+\log_2 b$ (a=b=n)")
ax.set_xlabel("n"); ax.set_ylabel("loop iterations"); ax.legend(); ax.set_title("Euclid's algorithm: logarithmic running time")
fig.tight_layout(); fig.savefig(IMG+"euclid_steps.png", dpi=160); plt.close(fig)

# Fig 2: modular multiplication tables mod 7 vs mod 6
fig, axs = plt.subplots(1,2, figsize=(9,4))
for ax, m in zip(axs,(7,6)):
    T = np.array([[(i*j)%m for j in range(m)] for i in range(m)])
    ax.imshow(T, cmap="viridis")
    for i in range(m):
        for j in range(m):
            ax.text(j,i,T[i,j],ha="center",va="center",color="w",fontsize=9)
    ax.set_title(f"multiplication mod {m}"); ax.set_xticks(range(m)); ax.set_yticks(range(m))
fig.tight_layout(); fig.savefig(IMG+"mult_tables.png", dpi=160); plt.close(fig)

# Fig 3: totient
def phi(n): return sum(1 for k in range(n) if gcd(k,n)==1) if n>1 else 1
ns = range(1,201)
fig, ax = plt.subplots(figsize=(6,3.6))
ax.scatter(ns,[phi(n) for n in ns],s=8,label=r"$\varphi(n)$")
ax.plot(ns,[n-1 for n in ns],"r--",lw=.8,label=r"$n-1$ (primes sit on this line)")
ax.set_xlabel("n"); ax.legend(); ax.set_title("Euler's totient function")
fig.tight_layout(); fig.savefig(IMG+"totient.png", dpi=160); plt.close(fig)
print("phi(10)=",phi(10), " phi(55)=",phi(55))

# Euler / Fermat check
for n in range(2,300):
    for a in range(1,n):
        if gcd(a,n)==1: assert pow(a,phi(n),n)==1

# Fig 4: CRT grid mod 3 x mod 5 -> mod 15
fig, ax = plt.subplots(figsize=(5,4))
G = np.zeros((3,5),int)
for r in range(15): G[r%3, r%5] = r
ax.imshow(G, cmap="Blues")
for i in range(3):
    for j in range(5): ax.text(j,i,G[i,j],ha="center",va="center")
ax.set_xlabel("r mod 5"); ax.set_ylabel("r mod 3"); ax.set_title("CRT: each pair appears exactly once (mod 15)")
fig.tight_layout(); fig.savefig(IMG+"crt_grid.png", dpi=160); plt.close(fig)

# ---------- Day 2: primes ----------
def sieve(n):
    s = bytearray([1])*(n+1); s[0:2]=b"\0\0"
    for i in range(2, isqrt(n)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P = sieve(100000)
print("pi(10^5) =", len(P), " x/ln x =", int(1e5/math.log(1e5)))

xs = np.arange(10, 100001, 50)
import bisect
pix = [bisect.bisect_right(P,x) for x in xs]
fig, ax = plt.subplots(figsize=(6,3.6))
ax.plot(xs,pix,label=r"$\pi(x)$"); ax.plot(xs,xs/np.log(xs),"--",label=r"$x/\ln x$")
ax.set_xlabel("x"); ax.legend(); ax.set_title("Prime Number Theorem")
fig.tight_layout(); fig.savefig(IMG+"pnt.png", dpi=160); plt.close(fig)

# Ulam spiral
N=201; M=np.zeros((N,N)); x=y=N//2; pset=set(sieve(N*N))
dirs=[(1,0),(0,-1),(-1,0),(0,1)]; n=1; step=1; di=0
M[y,x]= 0
while n<N*N:
    for _ in range(2):
        for _ in range(step):
            if n>=N*N: break
            x+=dirs[di][0]; y+=dirs[di][1]; n+=1
            if 0<=x<N and 0<=y<N: M[y,x]= n in pset
        di=(di+1)%4
    step+=1
fig, ax = plt.subplots(figsize=(4.5,4.5)); ax.imshow(M,cmap="gray_r",interpolation="nearest"); ax.axis("off"); ax.set_title("Ulam spiral: primes line up on diagonals")
fig.tight_layout(); fig.savefig(IMG+"ulam.png", dpi=160); plt.close(fig)

# Fermat test / pseudoprimes / Carmichael
def fermat_psp(base, lim):
    return [n for n in range(3,lim) if n not in set(P) and pow(base,n-1,n)==1] if lim<=100000 else None
Pset=set(P)
def is_carm(n):
    if n in Pset or n<3: return False
    return all(pow(a,n-1,n)==1 for a in range(2,n) if gcd(a,n)==1)
print("psp base2 <1000:", [n for n in range(3,1000) if n not in Pset and pow(2,n-1,n)==1])
print("341 base3 ok?", pow(3,340,341), " 91 base2?", pow(2,90,91), " 91 base3?", pow(3,90,91))
carm=[n for n in range(3,3000) if is_carm(n)]
print("Carmichael <3000:", carm)

# Miller-Rabin
def miller_rabin(n, bases=(2,3,5,7,11,13,17,19,23,29,31,37)):
    if n<2: return False
    for p in bases:
        if n%p==0: return n==p
    d,s=n-1,0
    while d%2==0: d//=2; s+=1
    for a in bases:
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True
assert all(miller_rabin(n)==(n in Pset) for n in range(2,5000))
assert not miller_rabin(561) and not miller_rabin(1729)
print("2^127-1 prime (MR):", miller_rabin(2**127-1))

# Lucas-Lehmer
def lucas_lehmer(p):
    if p==2: return True
    s=4; M=2**p-1
    for _ in range(p-2): s=(s*s-2)%M
    return s==0
print("Mersenne exponents <130:", [p for p in P if p<130 and lucas_lehmer(p)])
print("perfect numbers:", [2**(p-1)*(2**p-1) for p in P if p<14 and lucas_lehmer(p)])

# Day 2 assignment Q1: repunits
def repunit_mod(p,n): return (pow(10,n,9*p)-1)//9 % p   # R_n mod p
for p in [3,7,11,13,37]:
    ns=[n for n in range(1,60) if repunit_mod(p,n)==0]
    print("p=",p,"A_p start:",ns[:6])
assert all(repunit_mod(2,n)!=0 and repunit_mod(5,n)!=0 for n in range(1,100))

# Fig: Fermat/Carmichael witnesses
fig, ax = plt.subplots(figsize=(6,3.4))
ns=[n for n in range(3,1800) if n not in Pset]
frac=[sum(1 for a in range(2,n) if gcd(a,n)==1 and pow(a,n-1,n)==1)/max(1,phi(n)) for n in ns]
ax.scatter(ns,frac,s=3); ax.set_xlabel("composite n"); ax.set_ylabel("fraction of liars (coprime a)")
ax.set_title("Fermat liars: 561, 1105 and 1729 hit 1.0")
fig.tight_layout(); fig.savefig(IMG+"fermat_liars.png", dpi=160); plt.close(fig)

# ---------- Crypto ----------
# XOR
print("XOR demo:", format(0b011011 ^ 0b010101, "06b"))

# Caesar + frequency analysis
ENG = {'e':12.7,'t':9.1,'a':8.2,'o':7.5,'i':7.0,'n':6.7,'s':6.3,'h':6.1,'r':6.0,'d':4.3,'l':4.0,'c':2.8,'u':2.8,'m':2.4,'w':2.4,'f':2.2,'g':2.0,'y':2.0,'p':1.9,'b':1.5,'v':1.0,'k':0.8,'j':0.15,'x':0.15,'q':0.1,'z':0.07}
def caesar(t,k): return "".join(chr((ord(c)-97+k)%26+97) if c.isalpha() else c for c in t.lower())
def vigenere(t,key,dec=False,off=0):
    out=[];i=0
    for c in t.lower():
        if c.isalpha():
            s=ord(key.lower()[i%len(key)])-97+off; out.append(chr((ord(c)-97+(-s if dec else s))%26+97)); i+=1
        else: out.append(c)
    return "".join(out)
print(vigenere("THISTHEMESSAGE","ABCD",off=1).upper())  # handout: A=1,B=2,...
def crack_caesar(ct):
    best=None
    for k in range(26):
        pt=caesar(ct,-k)
        score=sum(ENG.get(c,0) for c in pt if c.isalpha())
        if best is None or score>best[0]: best=(score,k,pt)
    return best
text=("cryptography is the study of secure communication in the presence of adversaries the attacker knows the system "
      "except for the key and that is the whole idea behind kerckhoffs principle modern ciphers rely on mathematical problems "
      "that are easy in one direction and hard in the other")
ct=caesar(text,3)
sc,k,pt=crack_caesar(ct); print("Caesar recovered shift",k)
def freq(t):
    letters=[c for c in t if c.isalpha()]; n=len(letters)
    return [100*letters.count(chr(97+i))/n for i in range(26)]
fig, axs = plt.subplots(1,3, figsize=(12,3.2), sharey=True)
axs[0].bar(range(26),[ENG.get(chr(97+i),0) for i in range(26)]); axs[0].set_title("English (reference)")
axs[1].bar(range(26),freq(ct),color="C1"); axs[1].set_title("Caesar, shift 3: fingerprint shifted")
vk="lemon"; vct=vigenere(text*3,vk)
axs[2].bar(range(26),freq(vct),color="C2"); axs[2].set_title("Vigenère: flatter fingerprint")
for a in axs: a.set_xticks(range(26)); a.set_xticklabels([chr(97+i) for i in range(26)],fontsize=6)
axs[0].set_ylabel("% of letters")
fig.tight_layout(); fig.savefig(IMG+"freq.png", dpi=160); plt.close(fig)
# Vigenere key length via index of coincidence
def ic(t):
    n=len(t); return sum(t.count(chr(97+i))*(t.count(chr(97+i))-1) for i in range(26))/(n*(n-1))
vtext="".join(c for c in vct if c.isalpha())
ics={L:np.mean([ic(vtext[i::L]) for i in range(L)]) for L in range(1,13)}
print("IC by key length:", {k:round(v,3) for k,v in ics.items()})
fig, ax = plt.subplots(figsize=(5.5,3.2)); ax.bar(list(ics),list(ics.values())); ax.axhline(0.0667,color="r",ls="--",label="English ≈ 0.067"); ax.axhline(0.0385,color="gray",ls=":",label="random ≈ 0.038")
ax.set_xlabel("guessed key length L"); ax.set_ylabel("mean index of coincidence"); ax.legend(); ax.set_title("Breaking Vigenère: spike at L = 5 (key 'lemon')")
fig.tight_layout(); fig.savefig(IMG+"vigenere_ic.png", dpi=160); plt.close(fig)
assert vigenere(vigenere(text,"lemon"),"lemon",dec=True)==text.lower()

# Many-time pad
def xor_bytes(a,b): return bytes(x^y for x,y in zip(a,b))
m1=b"steal the secret"; m2=b"the cat the bird"; key=b"supersecretverys"
c1,c2=xor_bytes(m1,key),xor_bytes(m2,key)
assert xor_bytes(c1,c2)==xor_bytes(m1,m2)
x12=xor_bytes(c1,c2)
for crib in (b" the ", b"at the "):
    hits=[]
    for pos in range(len(x12)-len(crib)+1):
        guess=xor_bytes(crib,x12[pos:pos+len(crib)])
        if all(chr(g).isalpha() or chr(g)==" " for g in guess): hits.append((pos,guess))
    print("crib",crib,hits)

# Diffie-Hellman
N,g,a,b=11,7,2,4
A,B=pow(g,a,N),pow(g,b,N); print("DH:",A,B,pow(B,a,N),pow(A,b,N))
def order(x,n): return next(k for k in range(1,n+1) if pow(x,k,n)==1)
print("primitive roots mod 5:", [x for x in range(1,5) if order(x,5)==4], " mod 11:", [x for x in range(1,11) if order(x,11)==10])
fig, ax = plt.subplots(figsize=(5.5,3.4))
for x in (2,3,7):
    ax.plot(range(1,11),[pow(x,k,11) for k in range(1,11)],"o-",label=f"{x}^k mod 11")
ax.set_xlabel("k"); ax.legend(); ax.set_title("Powers mod 11 cycle through the group (g=7 is a generator)")
fig.tight_layout(); fig.savefig(IMG+"dh_powers.png", dpi=160); plt.close(fig)

# RSA toy
p,q,e=5,11,7; Nn=p*q; ph=(p-1)*(q-1); d=mod_inverse(e,ph)
print("RSA toy: N,phi,d =",Nn,ph,d, " c=",pow(2,e,Nn)," m=",pow(pow(2,e,Nn),d,Nn))
def to_int(s):
    r=0
    for ch in s: r=r*256+ord(ch)
    return r
def to_str(n):
    r=""
    while n>0: r+=chr(n%256); n//=256
    return r[::-1]
def Encrypt(msg,mod,exp): return fast_pow(to_int(msg),exp,mod)
def Decrypt(c,p,q,e):
    d=mod_inverse(e,(p-1)*(q-1)); return to_str(fast_pow(c,d,p*q))
P1,Q1=1000000007,1000000009
ct=Encrypt("attack",P1*Q1,23917); print("RSA 'attack' ->",ct,"->",Decrypt(ct,P1,Q1,23917))
assert Decrypt(ct,P1,Q1,23917)=="attack"

# Q3 DecipherSimple (small message space -> encrypt every candidate)
def DecipherSimple(c,mod,exp,cands):
    for m in cands:
        if Encrypt(m,mod,exp)==c: return m
    return "don't know"
mod3,exp3=1110014195838866450995043,767549
print("Q3:",DecipherSimple(Encrypt("wait",mod3,exp3),mod3,exp3,["attack","don't attack","wait"]))

# Q4 small prime
def DecipherSmallPrime(c,mod,exp):
    for p in range(2,10**6):
        if mod%p==0: return Decrypt(c,p,mod//p,exp)
    return "don't know"
qbig=1000000007*1000000009*1000003
pp=999983
mod4=pp*(10**30+57) if False else pp*1000000007*1000000009
print("Q4 demo:",DecipherSmallPrime(Encrypt("attack",pp*1000000009,239),pp*1000000009,239))

# Q5 Fermat factorisation
def fermat_factor(n):
    a=isqrt(n)+(0 if isqrt(n)**2==n else 1)
    while True:
        b2=a*a-n; b=isqrt(b2)
        if b*b==b2: return a-b,a+b
        a+=1
f=fermat_factor(1000000007*1000000009); print("Fermat factor:",f)
def DecipherSmallDiff(c,mod,exp):
    p,q=fermat_factor(mod); return Decrypt(c,p,q,exp)
n5=1000000007*1000000009
print("Q5:",DecipherSmallDiff(Encrypt("attack",n5,239),n5,239))

# Q6 common factor
def DecipherCommonDivisor(c1,m1,e1,c2,m2,e2):
    p=gcd(m1,m2); return Decrypt(c1,p,m1//p,e1),Decrypt(c2,p,m2//p,e2)
pc=1000000007; qa=1000000009; qb=998244353
ma,mb=pc*qa,pc*qb
print("Q6:",DecipherCommonDivisor(Encrypt("attack",ma,239),ma,239,Encrypt("wait",mb,65537),mb,65537))

# Q7 Hastad e=2
def DecipherHastad(c1,n1,c2,n2):
    m2=crt(max(n1,n2), c1 if n1>n2 else c2, min(n1,n2), c2 if n1>n2 else c1)
    return to_str(isqrt(m2))
p1,q1=790383132652258876190399065097,662503581792812531719955475509
p2,q2=656917682542437675078478868539,1263581691331332127259083713503
n1,n2=p1*q1,p2*q2
print("Q7:",DecipherHastad(Encrypt("attack",n1,2),n1,Encrypt("attack",n2,2),n2))

# Shor classical part
def shor_classical(N,a):
    assert gcd(a,N)==1
    r=order(a,N); 
    if r%2 or pow(a,r//2,N)==N-1: return None
    return r,gcd(pow(a,r//2)-1,N),gcd(pow(a,r//2)+1,N)
print("Shor(35,a=8):",shor_classical(35,8))
# Fig RSA timing: brute-force factoring vs modular exp
import time
bits=list(range(16,44,3)); tf=[];te=[]
for b in bits:
    pr=[x for x in P if x.bit_length()==b//2]
    if not pr:
        import sympy; pr=[sympy.nextprime(2**(b//2-1)), sympy.nextprime(2**(b//2-1)+1000)]
    n=pr[-1]*pr[-2]
    t=time.perf_counter(); 
    k=2
    while n%k: k+=1
    tf.append(time.perf_counter()-t)
    t=time.perf_counter(); pow(12345,n-1,n); te.append(time.perf_counter()-t)
fig, ax = plt.subplots(figsize=(6,3.6))
ax.semilogy(bits,tf,"o-",label="trial-division factoring"); ax.semilogy(bits,[max(x,1e-7) for x in te],"s-",label="modular exponentiation")
ax.set_xlabel("bits of N"); ax.set_ylabel("seconds"); ax.legend(); ax.set_title("Easy one way, hard the other way")
fig.tight_layout(); fig.savefig(IMG+"one_way.png", dpi=160); plt.close(fig)

# Brute-force table
for name,n in (("DES",56),("AES",128),("3DES",168)):
    t=2**(n-1)/1e18; print(name, f"{2**(n-1):.2e} keys, {t:.2e} s, {t/(3600*24*365):.2e} yr")
# Entropy teasers
print("5 dice entropy bits:",5*math.log2(6))
print("ALL DONE")

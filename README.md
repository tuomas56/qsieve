# QSieve

A simple implementation of the quadratic sieve in Rust for integers `n < (1 << 128)`. 
Currently it is very slow (at least two order of magnitude slower than msieve on all inputs, and very slow in some cases).

Most cases finish fairly fast for inputs of about 20-30 digits (< 0.5 seconds), but occasionally the factor base picking algorithm will get stumped and it can take upwards of ~30 seconds. This has still got to be tuned in multiple ways. For one, instead of solving the congruences modulo prime powers and sieving, we just divide out all the powers. In the best case this is twice as slow.

Optimizations still to do:

* Single or Double prime bounds.
* Log approximations instead of divisions.
* Pollard-Rho or EC exponent vectors.
* Block Lanzcos nullspace algorithm.
* Parallel sieving.
* Multiple polynomial sieving.

### Building

You need a nightly rust compiler.

```
$ cargo build --release
```

Then:

```
$ ./target/release/qsieve N1 N2 ... Nk
```

where `N1` to `Nk` are positive integers of 128 bits or less.


Example usage:

```
$ time ./target/release/qsieve 34028236692093846346337460743176898
trial division up to 20 bits
factor: 2
factor: 31
factor: 131
factor: 2377
factor: 93133
miller-rabin test => false
quadratic sieve on 18925339779841046878049
00000 primes remaining
factor base: [2, 5 ... 5189]
smooth: [(137569399889, 5978546334272), (137569400110, 66784221134051) ... (137570481986, 297734019303626147)]
matrix: 454 x 354
deps: [.. 100 ..]
non-trivial: b3961^2 = b3035^2 (mod 18925339779841046878049) 
trial division up to 20 bits
factor: 1291567
done.
miller-rabin test => false
quadratic sieve on 14653006603483247
00000 primes remaining
factor base: [2, 13 ... 1489]
smooth: [(121049659, 13340533034), (121050349, 180389538554) ... (121202971, 37153575743594)]
matrix: 147 x 119
deps: [.. 20 ..]
non-trivial: b1156^2 = b926^2 (mod 14653006603483247) 
trial division up to 20 bits
factor: 817727017
done.
factor: 17919191
done.
34028236692093846346337460743176898 => [2, 31, 131, 2377, 93133, 1291567, 817727017, 17919191]

real	0m0.348s
user	0m0.316s
sys	    0m0.024s
```




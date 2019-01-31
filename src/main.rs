#![feature(type_ascription)]

extern crate ramp;

use std::ops;
use std::env;
use ramp::Int;

struct PrimeSieve {
    data: Vec<bool>,
    i: usize
}

impl PrimeSieve {
    fn new(n: usize) -> PrimeSieve {
        let mut data = vec![false; n];
        data[0] = true;
        data[1] = true;
        PrimeSieve { data, i: 0 }
    }
}

impl Iterator for PrimeSieve {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.i < self.data.len() && self.data[self.i] {
            self.i += 1;
        }

        if self.i == self.data.len() {
            return None
        }

        self.data[self.i] = true;
        for i in ((2*self.i)..self.data.len()).step_by(self.i) {
            self.data[i] = true;
        }

        Some(self.i)
    }
}

#[derive(Clone, Debug)]
struct FactorBase {
    oprimes: Vec<usize>,
    rprimes: Vec<usize>
}

impl FactorBase {
    fn new(n: u128, s: usize) -> FactorBase {
        let x = (2 * n.isqrt()) as f64;
        let b = ((0.5 - 0.1) * x.ln() * x.ln().ln()).sqrt().exp() as usize + s;
        let mut k = 2 * b;
        let mut primes = Vec::new();
        while primes.len() < b {
            primes = PrimeSieve::new(k).filter(|p| GFp::new(n, *p as u128).is_quadratic_residue()).collect();
            k += b;
        }
        primes.truncate(b);
        FactorBase { oprimes: primes.clone(), rprimes: primes }
    }

    fn pop(&mut self) -> Option<usize> {
        self.rprimes.pop()
    }

    fn len(&self) -> usize {
        self.rprimes.len()
    }
    
    fn exponent_vector(&self, n: u128) -> Vec<u8> {
        let mut res = Vec::new();
        for &prime in &self.oprimes {
            res.push({
                let mut s = 0;
                let mut n = n;
                while n % (prime as u128) == 0 {
                    s += 1;
                    n /= prime as u128;
                }
                s & 1
            });
        }
        res
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
struct GFp {
    n: u128,
    p: u128
}

impl GFp {
    fn new(n: u128, p: u128) -> GFp {
        GFp { n: n % p, p }
    }

    fn pow(mut self, mut n: u128) -> GFp {
        let mut y = GFp::new(1, self.p);

        while n > 1 {
            if n & 1 == 1 {
                y *= self;
            }

            self *= self;
            n >>= 1;
        }

        self * y
    }

    fn legendre_symbol(self) -> GFp {
        self.pow(self.p >> 1)
    }

    fn is_quadratic_residue(self) -> bool {
        self.legendre_symbol().n == 1
    }

    fn sqrt(self) -> (Option<GFp>, Option<GFp>) {
        if self.p == 2 {
            return (Some(self), None)
        }

        if !self.is_quadratic_residue() {
            return (None, None)
        }

        if self.n == 0 {
            return (Some(self), None)
        }

        if self.n == 1 {
            return (Some(self), Some(-self))
        }

        let mut a = GFp::new(0, self.p);
        for i in 2..self.p {
            let i = GFp::new(i, self.p);
            if !(i*i - self).is_quadratic_residue() {
                a = i;
                break;
            }
        }

        let r = a*a - self;
        let xp = GFp2::new(a.n, 1, r.n, self.p);
        let x = xp.pow((self.p >> 1) + 1).a;

        (Some(x), Some(-x))
    }
}

impl ops::Add for GFp {
    type Output = GFp;

    fn add(self, other: GFp) -> GFp {
        if other.n == 0 {
            self
        } else {
            self - (-other)
        }
    }
}

impl ops::AddAssign for GFp {
    fn add_assign(&mut self, other: GFp) {
        if other.n != 0 {
            *self -= -other;
        }
    }
}

impl ops::Neg for GFp {
    type Output = GFp;

    fn neg(self) -> GFp {
        GFp { n: self.p - self.n, p: self.p }
    }
}

impl ops::Sub for GFp {
    type Output = GFp;

    fn sub(self, other: GFp) -> GFp {
        if self.n >= other.n {
            GFp { n: self.n - other.n, p: self.p }
        } else {
            GFp { n: other.p - other.n + self.n, p: self.p }
        }
    }
}

impl ops::SubAssign for GFp {
    fn sub_assign(&mut self, other: GFp) {
        *self = *self - other;
    }
}

impl ops::Div for GFp {
    type Output = GFp;

    fn div(self, other: GFp) -> GFp {
        fn extended_gcd(a: u128, b: u128) -> (u128, u128, u128) {
            let mut s = 0;    
            let mut old_s = 1;
            let mut t = 1;    
            let mut old_t = 0;
            let mut r = b;    
            let mut old_r = a;
            
            while r != 0 {
                let quotient = old_r / r;
                let q = r;
                r = old_r - quotient * r;
                old_r = q;
                let q = s;
                s = old_s - quotient * s;
                old_s = q;
                let q = t;
                t = old_t - quotient * t;
                old_t = q;
            }

            (old_r, old_s, old_t)
        }

        let (_, x, _) = extended_gcd(other.n, other.p);
        self * GFp::new(x, other.p)
    }
}


impl ops::Mul for GFp {
    type Output = GFp;

    fn mul(mut self, mut other: GFp) -> GFp {
        if self.n <= u128::max_value() >> 64 && other.n <= u128::max_value() >> 64 {
            GFp { n: (self.n * other.n) % self.p, p: self.p }
        } else {
            let mut res = GFp::new(0, self.p);

            while other.n > 0 {
                if other.n & 1 == 1 {
                    res += self;
                }

                self += self;
                other.n >>= 1;
            }

            res
        }
    }
}

impl ops::MulAssign for GFp {
    fn mul_assign(&mut self, other: GFp) {
        *self = *self * other;
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
struct GFp2 {
    a: GFp,
    b: GFp,
    r: GFp
}

impl GFp2 {
    fn new(a: u128, b: u128, r: u128, p: u128) -> GFp2 {
        GFp2 { a: GFp::new(a, p), b: GFp::new(b, p), r: GFp::new(r, p) }
    }

    fn pow(mut self, mut n: u128) -> GFp2 {
        let mut y = GFp2::new(1, 0, self.r.n, self.r.p);

        while n > 0 {
            if n & 1 == 1 {
                y *= self;
            }

            self *= self;
            n >>= 1;
        }

        y
    }
}

impl ops::Add for GFp2 {
    type Output = GFp2;

    fn add(self, other: GFp2) -> GFp2 {
        GFp2 { a: self.a + other.a, b: self.b + other.b, r: self.r }
    }
}

impl ops::AddAssign for GFp2 {
    fn add_assign(&mut self, other: GFp2) {
        *self = *self + other;
    }
}

impl ops::Neg for GFp2 {
    type Output = GFp2;

    fn neg(self) -> GFp2 {
        GFp2 { a: -self.a, b: -self.b, r: self.r }
    }
}

impl ops::Sub for GFp2 {
    type Output = GFp2;

    fn sub(self, other: GFp2) -> GFp2 {
        self + (-other)
    }
}

impl ops::SubAssign for GFp2 {
    fn sub_assign(&mut self, other: GFp2) {
        *self = *self - other;
    }
}

impl ops::Mul for GFp2 {
    type Output = GFp2;

    fn mul(self, other: GFp2) -> GFp2 {
        GFp2 { a: self.a * other.a + self.r * self.b * other.b, b: self.a * other.b + self.b * other.a, r: self.r }
    }
}

impl ops::MulAssign for GFp2 {
    fn mul_assign(&mut self, other: GFp2) {
        *self = *self * other;
    }
}

trait NumTools {
    fn isqrt(self) -> Self;
    fn gcd(self, other: Self) -> Self;
    fn clog2(self) -> Self;
}

impl NumTools for u128 {
    fn isqrt(self) -> u128 {
        let mut a = 0;
        let mut b = self;
        let mut m = (a + b + 1) >> 1;
        
        while m*m != self && b - a > 1 {
            if m*m > self {
                b = m;
            } else {
                a = m;
            }

            m = (a + b + 1) >> 1;
        }

        m
    }

    fn gcd(mut self, mut other: u128) -> u128 {
        while other != 0 {
            let t = other;
            other = self % other;
            self = t;
        }

        self
    }

    fn clog2(mut self) -> u128 {
        let mut l = 1;

        while self > 0 {
            self >>= 1;
            l += 1;
        }

        l
    }
}

impl NumTools for Int {
    fn isqrt(self) -> Int {
        let mut a = Int::zero();
        let mut b = self.clone();
        let mut m: Int = (a.clone() + b.clone() + 1) / 2;
        
        while m.clone() * m.clone() != self && b.clone() - a.clone() > 1 {
            if m.clone() * m.clone() > self {
                b = m.clone();
            } else {
                a = m.clone();
            }

            m = (a.clone() + b.clone() + 1) / 2;
        }

        m
    }

    fn gcd(mut self, mut other: Int) -> Int {
        while other != 0 {
            let t = other.clone();
            other = self % other;
            self = t;
        }

        self
    }

    fn clog2(mut self) -> Int {
        let mut l = Int::one();

        while self > 0 {
            self >>= 1;
            l += 1;
        }

        l
    }
}

#[derive(Debug, Clone)]
struct QuadraticSieve {
    factor_base: FactorBase,
    data: Vec<u128>,
    n: u128,
    sqrtn: u128
}

impl QuadraticSieve {
    fn new(n: u128, factor_base: FactorBase) -> QuadraticSieve {
        let sqrtn = n.isqrt();
        let data = Vec::with_capacity((factor_base.len() as f64).powf(2.5) as usize);
        
        QuadraticSieve {
            n, factor_base, data, sqrtn
        }
    }

    fn initialize(&mut self) {
        for i in 0..self.data.capacity() {
            self.data.push((i as u128 + self.sqrtn) * (i as u128 + self.sqrtn) - self.n);
        }
    }

    fn process_prime(&mut self) {
        let p = self.factor_base.pop().unwrap();
        let (sn1, sn2) = GFp::new(self.n, p as u128).sqrt();

        if let Some(sn1) = sn1 {
            let sn1 = (sn1 - GFp::new(self.sqrtn, p as u128)).n as usize;
            for k in 1..5 {
                for i in (sn1..self.data.len()).step_by(p.pow(k)) {
                    self.data[i] /= p as u128;
                }
            }
        }

        if let Some(sn2) = sn2 {
            let sn2 = (sn2 - GFp::new(self.sqrtn, p as u128)).n as usize;
            for k in 1..3 {
                for i in (sn2..self.data.len()).step_by(p.pow(k)) {
                    self.data[i] /= p as u128;
                }
            }
        }
    }

    fn extract_smooth(&self) -> Vec<(u128, u128)> {
        let mut res = Vec::new();

        for (i, &k) in self.data.iter().enumerate() {
            if k == 1 {
                res.push((i as u128 + self.sqrtn, (i as u128 + self.sqrtn) * (i as u128 + self.sqrtn) - self.n));
            }
        }
        
        res
    }

    fn done(&self) -> bool {
        self.factor_base.len() == 0
    }
}

#[derive(PartialEq, Clone, Debug)]
struct Mat2 {
    width: usize,
    height: usize,
    data: Vec<u8>
}

impl Mat2 {
    fn from_columns(columns: &Vec<Vec<u8>>) -> Mat2 {
        Mat2 {
            width: columns.len(),
            height: columns[0].len(),
            data: (0..columns[0].len()).map(|i| {
                columns.iter().map(|column| column[i]).collect::<Vec<_>>()
            }).flatten().collect()
        }
    }

    fn from_rows(rows: &Vec<Vec<u8>>) -> Mat2 {
        Mat2 {
            width: rows[0].len(),
            height: rows.len(),
            data: rows.iter().cloned().flatten().collect()
        }
    }

    fn identity(n: usize) -> Mat2 {
        Mat2::from_columns(&(0..n).map(|i| {
            (0..n).map(|j| if i == j { 1 } else { 0 }).collect()
        }).collect())
    }

    fn columns(&self) -> Vec<Vec<u8>> {
        (0..self.width).map(|col| {
            (0..self.height).map(|row| {
                self.data[self.width*row + col]
            }).collect::<Vec<_>>()
        }).collect()
    }

    fn rows(&self) -> Vec<Vec<u8>> {
        (0..self.height).map(|row| {
            (0..self.width).map(|col| {
                self.data[self.width*row + col]
            }).collect::<Vec<_>>()
        }).collect()
    }

    fn xor(&mut self, i: usize, j: usize) {
        for (ix, jx) in (self.width*i..self.width*(i + 1)).zip(self.width*j..self.width*(j + 1)) {
            self.data[ix] ^= self.data[jx];
        }
    }

    fn swap(&mut self, i: usize, j: usize) {
        if i != j {
            self.xor(i, j);
            self.xor(j, i);
            self.xor(i, j);
        }
    }

    fn augment_rows(&self, other: &Mat2) -> Mat2 {
        let scols = self.columns();
        let mut ocols = other.columns();
        Mat2::from_columns(&(0..self.width).map(|i| {
            let mut c = scols[i].clone();
            c.append(&mut ocols[i]);
            c
        }).collect())
    }

    fn augment_columns(&self, other: &Mat2) -> Mat2 {
        let srows = self.rows();
        let mut orows = other.rows();
        Mat2::from_rows(&(0..self.height).map(|i| {
            let mut c = srows[i].clone();
            c.append(&mut orows[i]);
            c
        }).collect())
    }

    fn gaussian_elimination(mut self) -> Mat2 {
        let mut pivot_row = 0;
        let mut pivot_col = 0;

        while pivot_row < self.width && pivot_col < self.height {
            let mut pivoted = false;
            for row in pivot_row..self.height {
                if self.data[row*self.width + pivot_col] == 1 {
                    self.swap(pivot_row, row);
                    pivoted = true;
                    break;
                }
            }

            if !pivoted {
                pivot_col += 1;
            }

            for row in (pivot_row + 1)..self.height {
                if self.data[row*self.width + pivot_col] == 1 {
                    self.xor(row, pivot_row);
                }
            }

            pivot_row += 1;
            pivot_col += 1;
        }

        pivot_col = 0;
        pivot_row = 0;

        while pivot_col < self.width && pivot_row < self.height {
            if self.data[pivot_row*self.width + pivot_col] == 0 {
                pivot_col += 1;
            } else {
                for row in 0..pivot_row {
                    if self.data[row*self.width + pivot_col] == 1 {
                        self.xor(row, pivot_row);
                    }
                }

                pivot_row += 1;
                pivot_col = 0;
            }
        }

        self
    }

    fn row_null_space(&self) -> Vec<Vec<u8>> {
        let mut res = Vec::new();
        let rn = self.width;
        let matrix = self.augment_columns(&Mat2::identity(self.height));
        let matrix = matrix.gaussian_elimination();
        for row in matrix.rows().iter() {
            let mut zero = true;
            for i in 0..rn {
                if row[i] != 0 {
                    zero = false;
                    break
                }
            }

            if zero {
                res.push(row[rn..row.len()].to_vec());
            }
        }

        res
    }

    fn print(&self) {
        for row in 0..self.height {
            for col in 0..self.width {
                print!("{} ", self.data[self.width*row + col]);
            }
            println!();
        }
    }
}

fn qsieve(n: u128) -> (u128, u128) {
    println!("quadratic sieve on {}", n);
    let mut s = 0;
    let mut sb = 100;
    loop {
        let mut smooth = Vec::new();
        let mut fb = FactorBase::new(n, s);
        while smooth.len() <= fb.len() {
            let mut qs = QuadraticSieve::new(n, fb.clone());
            qs.initialize();
            while !qs.done() {
                qs.process_prime();
                print!("\r{:05} primes remaining", qs.factor_base.len());
            }
            smooth = qs.extract_smooth();
            if smooth.len() <= fb.len() {
                s += 1;
                fb = FactorBase::new(n, s);
            }
        }
        println!();
        smooth.truncate(fb.len() + sb);
        println!("factor base: [{}, {} ... {}]", fb.rprimes[0], fb.rprimes[1], fb.rprimes[fb.rprimes.len() - 1]);
        println!("smooth: [{:?}, {:?} ... {:?}]", smooth[0], smooth[1], smooth[smooth.len() - 1]);
        let matrix = Mat2::from_rows(&smooth.iter().map(|(_, q)| {
            fb.exponent_vector(*q)
        }).collect());
        println!("matrix: {} x {}", matrix.height, matrix.width);
        let mut deps = matrix.row_null_space();
        println!("deps: [.. {} ..]", deps.len());
        deps.sort_by_key(|dep| dep.iter().cloned().sum() : u8);
        for dep in deps.iter() {
            let mut p = Int::one();
            let mut q2 = Int::one();
            for (i, &k) in dep.iter().enumerate()  {
                if k == 1 {
                    p *= Int::from(smooth[i].0);
                    q2 *= Int::from(smooth[i].1);
                }
            }
            assert_eq!((p.clone()*p.clone()) % Int::from(n), q2.clone() % Int::from(n));
            let q = q2.isqrt();
            let f: u128 = (&(p.clone().max(q.clone()) - p.clone().min(q.clone())).gcd(Int::from(n))).into();
            if f != 1 && f != n {
                if 3*p.clone().clog2() > 80 || 3*q.clone().clog2() > 80 {
                    println!("non-trivial: b{}^2 = b{}^2 (mod {}) ", p.clog2(), q.clog2(), n);
                } else {
                    println!("non-trivial: {}^2 = {}^2 (mod {}) ", p, q, n);
                }
                
                return (f, n/f)
            }
        }

        s += 1;
        sb += 100;
    }
}

fn prime(n: u128) -> bool {
    fn strong_fermat_pseudoprime(n: u128, a: u128) -> bool {
        let mut d2s = n - 1;
        let mut s = 0;
        while d2s & 1 == 0 {
            s += 1;
            d2s >>= 1;
        }
        let d = d2s;

        if GFp::new(a, n).pow(d).n == 1 {
            return true
        }

        for r in 0..s {
            if GFp::new(a, n).pow(d * (1 << r)).n == n - 1 {
                return true
            }
        }    

        false
    }

    if n < (1 << 64) {
        for &a in &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] {
            if !strong_fermat_pseudoprime(n, a) {
                println!("miller-rabin test => false");
                return false;
            }
        }

        return true;
    }

    for b in 2..((n.clog2() - 2).pow(2) << 2).min(n - 1) {
        if !strong_fermat_pseudoprime(n, b) {
            println!("miller-rabin test => false");
            return false
        }
    }

    println!("miller-rabin test => true");
    true
}

fn factorize(n: u128) -> Vec<u128> {
    if n < 4 {
        println!("trivial");
        return vec![n]
    }

    let on = n;
    let mut n = n;
    let mut res = Vec::new();
    println!("trial division up to 20 bits");
    let sieve = PrimeSieve::new(1 << 20);
    for p in sieve {
        let p = p as u128;
        while n % p == 0 {
            res.push(p);
            println!("factor: {}", p);
            n /= p;
        }

        if n == 1 {
            break;
        }
    }

    while n > 1 && !prime(n) {
        let (p, q) = qsieve(n);
        n = q;
        res.append(&mut factorize(p));
    }

    if n > 1 {
        println!("factor: {}", n);
        res.push(n);
    }

    assert_eq!(on, res.iter().cloned().product());

    println!("done.");

    res
}

fn main() {
    let mut res = Vec::new();
    for n in env::args().skip(1).map(|x| x.parse().unwrap()) {
        res.push((n, factorize(n)));
    }
    for (n, r) in res.iter() {
        println!("{} => {:?}", n, r);
    }
}
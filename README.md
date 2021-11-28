
<a id='Cyclotomics'></a>

<a id='Cyclotomics-1'></a>

# Cyclotomics

- [Cyclotomics](index.md#Cyclotomics)

<a id='Cyclotomics' href='#Cyclotomics'>#</a>
**`Cyclotomics`** &mdash; *Module*.



This  package deals with cyclotomic numbers,  the complex numbers which are linear  combinations  of  roots  of  unity  with  rational coefficients. It depends on the packages `ModuleElt` and `Primes`.

The  cyclotomic numbers form a field, the  cyclotomic field. It is also the maximal  extension of the rationals which  has an abelian Galois group. Its ring of integers are the sums of unity with integral coefficients.

Cyclotomics are very important for finite groups, since character values of finite groups are cyclotomic integers.

This package is a port of the GAP implementation of cyclotomics, which uses a  normal form  given by  writing them  in the  Zumbroich basis.  This form allows to find the smallest Cyclotomic field which contains a given number, and  decide in particular  if a cyclotomic  is zero. Let ζₙ=exp(2iπ/n). The Zumbroich  basis is  a particular  subset of  size φ(n) of 1,ζₙ,ζₙ²,…,ζₙⁿ⁻¹ which forms a basis of ℚ (ζₙ).

I  started  this  file  by  porting  Christian  Stump's Sage code, which is simpler to understand than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

As  does  GAP,  I  lower  automatically  numbers  after  each  computation; currently  the code is  somewhat slower (depending  on the operation it has the same speed or is slower up to 50%) than the C code in GAP but there are probably many opportunities to optimize that I missed.

What GAP does which I do not do is convert automatically a Cyclotomic which is  rational to a Rational,  a Rational which is  integral to an Integer, a BigInt  which is  small to  an Int,  etc… This  is a tremendously important optimization  but because of type stability in Julia it needs a new type of number to be added to Julia, which I am not competent enough to try.

We  define two types in  this package: `Root1` represents  a root of unity, and `Cyc` a cyclotomic number. The main way to build a Cyclotomic number is to  use the function `E(n,k=1)` which  constructs the `Root1` `ζₙᵏ`, and to make   linear  combinations  of  such  numbers  with  integer  or  rational coefficients.

**Examples**

```julia-repl
julia> E(3,2) # a root of unity
Root1: ζ₃²

julia> E(3)+E(4) # nice display at the repl
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹
```

```julia-rep1
julia> print(E(3)+E(4)) # otherwise give output which can be read back in
E(12,4)-E(12,7)-E(12,11)
```

```julia-repl
julia> E(12,11)-E(12,7) # square roots of integers are recognized
Cyc{Int64}: √3

# but you can prevent that (we assume julia 1.7)
julia> repr(E(12,11)-E(12,7),context=(:limit=>true,:quadratic=>false))
"-ζ₁₂⁷+ζ₁₂¹¹"

julia> a=E(3)+E(3,2)
Cyc{Int64}: -1

julia> conductor(a) # a has been lowered to ℚ (ζ₁)=ℚ 
1

julia> typeof(Int(a))
Int64
```

```julia-rep1
julia> Int(E(4))
ERROR: InexactError: convert(Int64, E(4))
```

```julia-repl
julia> inv(1+E(4)) # inverses often need Rational coefficients
Cyc{Rational{Int64}}: (1-ζ₄)/2

julia> inv(E(5)+E(5,4)) # but not always (we have here a unit)
Cyc{Int64}: -ζ₅²-ζ₅³

julia> Cyc(1//2+im) # we can convert to Cyclotomics Gaussian rationals
Cyc{Rational{Int64}}: (1+2ζ₄)/2

julia> conj(1+E(4)) # complex conjugate
Cyc{Int64}: 1-ζ₄

julia> real(E(3))  # real part
Cyc{Rational{Int64}}: -1/2

julia> Rational{Int}(real(E(3)))
-1//2

julia> imag(E(3))  # imaginary part
Cyc{Rational{Int64}}: √-3/2

julia> c=Cyc(E(9))   # an effect of the Zumbroich basis
Cyc{Int64}: -ζ₉⁴-ζ₉⁷

julia> Root1(c) # you can convert back to Root1 if possible
Root1: ζ₉

julia> Root1(1+E(4)) # the constructor returns nothing for a non-root
```

The  group of  roots of  unity is  isomorphic to  ℚ /ℤ  , thus  `Root1` are represented internally by a rational number in `[0,1[`.

```julia-repl
julia> Root1(;r=1//4) # this contructor ensures the fraction is in [0,1[
Root1: ζ₄

julia> c=E(4)*E(3) # faster computation if staying indside roots of unity
Root1: ζ₁₂⁷

julia> c=Complex{Float64}(E(3))  # convert to Complex{float} is sometimes useful
-0.4999999999999999 + 0.8660254037844387im
```

In  presence of a float of type  `T`, a `Cyc` is converted to `Complex{T}`. In  presence of a `Cyc`,  an integer, a rational,  or a complex of these is converted to a `Cyc`.

```julia-repl
julia> 0.0+E(3)
-0.4999999999999999 + 0.8660254037844387im

julia> E(3)+1//2
Cyc{Rational{Int64}}: √-3/2

julia> E(3)+im
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

```

The  function  `complex`  converts  a  `Cyc{T}`  to  a  `Complex{T}` if the conductor is 1 or 4, to a `Complex{float(T)}` otherwise.

```julia-repl
julia> complex(E(4))
0 + 1im

julia> complex(E(3))
-0.4999999999999999 + 0.8660254037844387im
```

`Cyc`s have methods `copy, hash, ==, cmp, isless` (total order) so they can be  keys in hashes or  elements of sets. Cyclotomics  which are integers or rationals  compare  correctly  to  `Real`s  (contrary  to  irrational  real `Cyc`s):

```julia-repl
julia> -1<Cyc(0)<1
true
```

You  can pick apart a cyclotomic in various ways. The fastest is to use the iterator  `pairs` which, for a cyclotomic  `a` of conductor `e` iterates on the  pairs `(i,c)` such that  `a` has a non-zero  coefficient `c` on `ζₑⁱ`. You  can also get the coefficient `ζₑⁱ` by indexing `a[i]` but it is slower than  `pairs` to iterate on coefficients this  way. Finally you can get the vector of all coefficients by `coefficients`.

```julia-repl
julia> a=E(3)+E(4)
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

julia> collect(pairs(a))
3-element Vector{Tuple{Int64, Int64}}:
 (4, 1)
 (7, -1)
 (11, -1)

julia> a[6],a[7]
(0, -1)

julia> coefficients(a)
12-element Vector{Int64}:
  0
  0
  0
  0
  1
  0
  0
 -1
  0
  0
  0
 -1
```

For more information see the methods denominator, Quadratic, galois, root. 

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime Cyclotomics.testmat(12)^2;  # on Julia 1.7
  315.271 ms (4331137 allocations: 302.51 MiB)
```

The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```

testmat(12)^2 takes 0.35s in GAP3, 0.24s in GAP4


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L1-L215' class='documenter-source'>source</a><br>

<a id='Cyclotomics.E' href='#Cyclotomics.E'>#</a>
**`Cyclotomics.E`** &mdash; *Function*.



`E(n,p=1)` makes the `Root1` equal to `ζₙᵖ`


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L290' class='documenter-source'>source</a><br>

<a id='Cyclotomics.Root1' href='#Cyclotomics.Root1'>#</a>
**`Cyclotomics.Root1`** &mdash; *Type*.



`Root1(c)`

`c` should be a cyclotomic number (a `Cyc`), or a `Real`. `Root1` returns a `Root1` object containing the rational `e/n` with `0≤e<n` (that is, `e/n∈ ℚ /ℤ`) if `c==E(n,e)`, and `nothing` if `c` is not a root of unity.

```julia-repl
julia> r=Root1(-E(9,2)-E(9,5))
Root1: ζ₉⁸

julia> conductor(r)
9

julia> exponent(r)
8

julia> Cyc(r)
Cyc{Int64}: -ζ₉²-ζ₉⁵

julia> Root1(-E(9,4)-E(9,5)) # nothing
```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L1028-L1050' class='documenter-source'>source</a><br>

<a id='Cyclotomics.conductor' href='#Cyclotomics.conductor'>#</a>
**`Cyclotomics.conductor`** &mdash; *Function*.



`conductor(c::Cyc)`    `conductor(a::AbstractArray)`

returns the smallest positive integer  n such that `c∈ ℚ (ζₙ)` (resp. all elements of `a` are in `ℚ (ζₙ)`).

```julia-repl
julia> conductor(E(9))
9

julia> conductor([E(3),1//2,E(4)])
12
```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L401-L415' class='documenter-source'>source</a><br>

<a id='Cyclotomics.coefficients' href='#Cyclotomics.coefficients'>#</a>
**`Cyclotomics.coefficients`** &mdash; *Function*.



`coefficients(c::Cyc)`

for  a cyclotomic `c` of conductor `n`,  returns a vector `v` of length `n` such that `c==∑ᵢ vᵢ₋₁ ζⁱ`.

```julia-repl
julia> coefficients(Cyc(E(9)))
9-element Vector{Int64}:
  0
  0
  0
  0
 -1
  0
  0
 -1
  0
```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L446-L465' class='documenter-source'>source</a><br>

<a id='Base.denominator' href='#Base.denominator'>#</a>
**`Base.denominator`** &mdash; *Function*.



`denominator(c::Cyc{Rational})`

returns the smallest `d` such that `d*c` has integral coefficients (thus is an algebraic integer).


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L476-L481' class='documenter-source'>source</a><br>

<a id='Cyclotomics.galois' href='#Cyclotomics.galois'>#</a>
**`Cyclotomics.galois`** &mdash; *Function*.



galois(c::Cyc,n::Int) applies to c the galois automorphism   of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.   n should be prime to conductor(c).

**Examples**

```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
Cyc{Int64}: 1-ζ₄

julia> galois(root(5),2)==-root(5)
true
```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L979-L991' class='documenter-source'>source</a><br>

<a id='Cyclotomics.Quadratic' href='#Cyclotomics.Quadratic'>#</a>
**`Cyclotomics.Quadratic`** &mdash; *Type*.



`Quadratic(c::Cyc)` 

determines  if  `c`  lives  in  a  quadratic  extension  of  `ℚ`. The call `q=Quadratic(c)`  returns a  struct `Quadratic`  with fields  `q.a`, `q.b`, `q.root`,  `q.den` representing `c` as `(q.a + q.b root(q.root))//q.den` if such a representation is possible or returns `q===nothing` otherwise.

**Examples**

```julia-repl
julia> Quadratic(1+E(3))
(1+√-3)/2

julia> Quadratic(1+E(5))

```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L1094-L1110' class='documenter-source'>source</a><br>

<a id='Cyclotomics.root' href='#Cyclotomics.root'>#</a>
**`Cyclotomics.root`** &mdash; *Function*.



`root(x,n=2)`

computes  the `n`-th root of `x` when we know  how to do it. We know how to compute  `n`-th  roots  for  roots  of  unity, square roots of integers and `n`-th  roots  of  perfect  `n`-th  powers  of  integers or square roots of integers.

```julia-repl
julia> root(-1)
Cyc{Int64}: ζ₄

julia> root(E(4))
Root1: ζ₈

julia> root(27,6)
Cyc{Int64}: √3
```


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L1189-L1207' class='documenter-source'>source</a><br>

<a id='Cyclotomics.Elist' href='#Cyclotomics.Elist'>#</a>
**`Cyclotomics.Elist`** &mdash; *Function*.



Cyclotomics.Elist(n,i)  

expresses  ζₙⁱ  in  zumbroich_basis(n):  it  is  a  sum  of some ζₙʲ with   coefficients all 1 or all -1. The result is a Pair sgn=>inds where sgn is   true  if coefficients are all 1 and false otherwise, and inds is the list   of i in 0:n-1 such that ζₙⁱ occurs with a non-zero coefficient (the i in   1:n such that ζₙⁱ⁻¹.. for :vec and :svec)


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L488-L496' class='documenter-source'>source</a><br>

<a id='Cyclotomics.zumbroich_basis' href='#Cyclotomics.zumbroich_basis'>#</a>
**`Cyclotomics.zumbroich_basis`** &mdash; *Function*.



Cyclotomics.zumbroich_basis(n::Int) 

returns  the Zumbroich basis of  ℚ (ζₙ) as the  vector of i in 0:n-1 such   that `ζₙⁱ` is in the basis


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L421-L426' class='documenter-source'>source</a><br>

<a id='Cyclotomics.prime_residues' href='#Cyclotomics.prime_residues'>#</a>
**`Cyclotomics.prime_residues`** &mdash; *Function*.



the numbers less than n and prime to n 


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L267' class='documenter-source'>source</a><br>

<a id='Cyclotomics.factor' href='#Cyclotomics.factor'>#</a>
**`Cyclotomics.factor`** &mdash; *Function*.



`factor(n::Integer)` make `Primes.factor` fast for small Ints by memoizing it


<a target='_blank' href='https://github.com/jmichel7/Cyclotomics.jl/blob/2e6718f5a928a88933474195c81360db75083126/src/Cyclotomics.jl#L275-L278' class='documenter-source'>source</a><br>


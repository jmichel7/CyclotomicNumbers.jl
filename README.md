
<a id='CyclotomicNumbers'></a>

<a id='CyclotomicNumbers-1'></a>

# CyclotomicNumbers

- [CyclotomicNumbers](index.md#CyclotomicNumbers)

<a id='CyclotomicNumbers' href='#CyclotomicNumbers'>#</a>
**`CyclotomicNumbers`** &mdash; *Module*.



This  package deals with cyclotomic numbers,  the complex numbers which are linear  combinations  of  roots  of  unity  with  rational coefficients. It depends on the packages `ModuleElt` and `Primes`. It can also handle linear combinations  of roots of  unity with coefficients  in an arbitrary ring of real  numbers like  `Floats`. The  cyclotomic numbers  can be  converted to `Complex` numbers.

The  cyclotomic numbers form a field, the  cyclotomic field. It is also the maximal  extension of the rationals which  has an abelian Galois group. Its ring of integers are the sums of roots of unity with integer coefficients.

Cyclotomic  numbers are very  important for finite  groups, since character values of finite groups are cyclotomic integers.

This package is a port of the GAP implementation of cyclotomics, which uses a  normal form  given by  writing them  in the  Zumbroich basis.  This form allows to find the smallest Cyclotomic field which contains a given number, and  decide in particular  if a cyclotomic  is zero. Let ζₙ=exp(2iπ/n). The Zumbroich  basis is  a particular  subset of  size φ(n) of 1,ζₙ,ζₙ²,…,ζₙⁿ⁻¹ which forms a basis of ℚ (ζₙ).  The reference is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

I  started  this  file  by  porting  Christian  Stump's Sage code, which is simpler  to understand than GAP's C  code. 

As GAP does, I lower automatically numbers after each computation (that is, reduce  them to the smallest cyclotomic field where they belong). Currently the  code is somewhat  slower (depending on  the operation it  has the same speed or is slower up to 50%) than the C code in GAP but there are probably opportunities to optimize that I missed.

What GAP does which I do not do is convert automatically a Cyclotomic which is  rational to a Rational,  a Rational which is  integral to an Integer, a BigInt  which is  small to  an Int,  etc… This  is a tremendously important optimization  but because of type stability in Julia it needs a new type of number to be added to Julia, which I am not competent enough to try.

This package is similar (and mostly compatible) with Marek Kaluba's package `Cyclotomics`,  whose existence I discovered after writing this package. We discussed merging them but concluded it would be a lot of work for benefits which are not clear currently. Some differences are:

  * I define two types in this  package: `Root1` represents a root of unity, and  `Cyc` a cyclotomic number. The advantage of having a separate type for  roots of  unity is  that computations  are very  fast for them, of which  I take advantage  in the package  `CycPol` for polynomials whose zeros are roots of unity.
  * In Kaluba's package  numbers are not  systematically lowered but only on demand  (like  for  printing).  this  speeds  up some computations by a factor  approaching 2, but it also makes  some computations I had to do infeasible,  like the following one (which if not lowering involves too large fields); the answer is `-36ζ₃²`:

```julia-rep1
julia> prod(x->1-x,[E(3),E(3,2),E(6),E(6,5),E(8),E(8),E(8,5),E(8,7),E(9,2),
E(9,5),E(9,8),E(12,7),E(12,11),E(16),E(16,3),E(16,5),E(16,9),E(16,11),E(16,13),
E(18,5),E(18,5),E(18,11),E(18,11),E(18,17),E(18,17),E(21,2),E(21,5),E(21,8),
E(21,11),E(21,17),E(21,20),
E(27,2),E(27,5),E(27,8),E(27,11),E(27,14),E(27,17),E(27,20),E(27,23),E(27,26),
E(32,7),E(32,15),E(32,23),E(32,31),E(39),E(39,4),E(39,7),E(39,10),E(39,16),
E(39,19),E(39,22),E(39,25),E(39,28),E(39,31),E(39,34),E(39,37),E(42),E(42,13),
E(42,19),E(42,25),E(42,31),E(42,37),E(48,11),E(48,19),E(48,35),E(48,43),E(60,7),
E(60,19),E(60,31),E(60,43),E(78,5),E(78,11),E(78,17),E(78,23),E(78,29),E(78,35),
E(78,41),E(78,47),E(78,53),E(78,59),E(78,71),E(78,77),E(80,7),E(80,23),E(80,31),
E(80,39),E(80,47),E(80,63),E(80,71),E(80,79),E(88,3),E(88,19),E(88,27),E(88,35),
E(88,43),E(88,51),E(88,59),E(88,67),E(88,75),E(88,83),E(90),E(90,7),E(90,13),
E(90,19),E(90,31),E(90,37),E(90,43),E(90,49),E(90,61),E(90,67),E(90,73),
E(90,79),E(96,5),E(96,13),E(96,29),E(96,37),E(96,53),E(96,61),E(96,77),E(96,85),
E(104),E(104,9),E(104,17),E(104,25),E(104,33),E(104,41),E(104,49),E(104,57),
E(104,73),E(104,81),E(104,89),E(104,97),E(144),E(144,17),E(144,25),E(144,41),
E(144,49),E(144,65),E(144,73),E(144,89),E(144,97),E(144,113),E(144,121),
E(144,137),E(152,5),E(152,13),E(152,21),E(152,29),E(152,37),E(152,45),
E(152,53),E(152,61),E(152,69),E(152,77),E(152,85),E(152,93),E(152,101),
E(152,109),E(152,117),E(152,125),E(152,141),E(152,149),E(204,11),E(204,23),
E(204,35),E(204,47),E(204,59),E(204,71),E(204,83),E(204,95),E(204,107),
E(204,131),E(204,143),E(204,155),E(204,167),E(204,179),E(204,191),E(204,203)])
```

If you `develop` my package it is easy to use the strategy of lowering only on  demand or to use alternate implementations like dense vectors or sparse vectors (like `Cyclotomics`) –- I prepared boolean flags to choose various implementations in the code. I choosed the implementation with `ModuleElts` and systematic lowering as giving the best results.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)` which   constructs  the  `Root1`  equal  to   `ζₙᵏ`,  and  to  make  linear combinations of such numbers.

**Examples**

```julia-repl
julia> E(3,2) # a root of unity
Root1: ζ₃²

julia> E(3)+E(4) # nice display at the repl
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹
```

```julia-rep1
julia> print(E(3)+E(4)) # otherwise give output which can be parsed back
E(12,4)-E(12,7)-E(12,11)
```

```julia-repl
julia> E(12,11)-E(12,7) # square roots of integers are recognized on output
Cyc{Int64}: √3

# but you can prevent that recognition
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
julia> inv(1+E(4)) # like for Integers inverse involves floats
Cyc{Float64}: 0.5-0.5ζ₄

julia> 1//(1+E(4))  # but not written this way
Cyc{Rational{Int64}}: (1-ζ₄)/2

julia> Cyc(1//2+im) # we can convert Gaussian rationals to cyclotomics
Cyc{Rational{Int64}}: (1+2ζ₄)/2

julia> conj(1+E(4)) # complex conjugate
Cyc{Int64}: 1-ζ₄

julia> real(E(3))  # real part
Cyc{Rational{Int64}}: -1/2

julia> Rational{Int}(real(E(3)))
-1//2

julia> imag(E(3))  # imaginary part
Cyc{Rational{Int64}}: √3/2

julia> c=Cyc(E(9))   # the normal form in the Zumbroich basis has two terms
Cyc{Int64}: -ζ₉⁴-ζ₉⁷

julia> Root1(c) #  but you can convert back to Root1 if possible
Root1: ζ₉

julia> Root1(1+E(4)) # the constructor Root1 returns nothing for a non-root
```

The  group of  roots of  unity is  isomorphic to  ℚ /ℤ  , thus  `Root1` are represented internally by a rational number in `[0,1[`.

```julia-repl
julia> Root1(;r=1//4) # this constructor ensures the fraction is in [0,1[
Root1: ζ₄

julia> c=E(4)*E(3) # fast computation if staying inside roots of unity
Root1: ζ₁₂⁷
```

`Root1` have the same operations as `Cyc`, but are first converted to `Cyc` for  any  operation  other  than  `one,  isone,  *, ^, inv, conj, /, //`. A `Root1`  can be  raised to  a `Rational`  power, which  extracts a  root if needed.  A `Root1` can be  taken apart using `order`  and `exponent` –- if `a`  and `b` are prime to each other,  `a` is the order of `E(a,b)` and `b` (taken  mod `a`) is the exponent. Note  that the order is not the conductor since `E(6)==-E(3)` has order 6 and conductor 3.

```julia-repl
julia> c=Complex{Float64}(E(3))  # convert to Complex{float} is sometimes useful
-0.4999999999999999 + 0.8660254037844387im
```

In  presence  of  a  `Cyc`  a  number  `<:Real`  or  `<:Complex{<:Real}` is converted to a `Cyc`.

```julia-repl
julia> 0.0+E(3)
Cyc{Float64}: 1.0ζ₃

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

`Cyc`s have methods `copy, hash, ==, cmp, isless` (total order) so they can be  keys in hashes, elements of sets,  and can be sorted. Cyclotomics which are  integers  or  rationals  compare  correctly  to `Real`s (this does not extend to irrational real `Cyc`s):

```julia-repl
julia> -1<Cyc(0)<1
true
```

`Cyc`s have the operations `+, -, *, /, //, inv, ^, conj, abs2, abs, image, real, reim, isinteger, isreal, one, isone, zero, iszero, complex, adjoint`. Cyclotomics   with  rational  or  integer  coefficients  have  the  methods `numerator`  and  `denominator`:  a  `Cyc`  `x`  is a Cyclotomic integer if `denominator(x)==1`   and  then  `numerator(x)`   gives  the  corresponding `Cyc{<:Integer}`.

You  can pick apart a cyclotomic in various ways. The fastest is to use the iterator  `pairs` which, for a cyclotomic  `a` of conductor `e` iterates on the  pairs `(i,c)` such that  `a` has a non-zero  coefficient `c` on `ζₑⁱ`. You  can also get  the coefficient of  `ζₑⁱ` as `a[i]`  but it is slower to iterate  on coefficients  this way.  Finally you  can get (efficiently) the vector of all coefficients by `coefficients`.

```julia-repl
julia> a=E(3)+E(4)
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹

julia> collect(pairs(a))
3-element Vector{Pair{Int64, Int64}}:
  4 => 1
  7 => -1
 11 => -1

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

julia> valtype(a) # the type of the coefficients of a
Int64
```

For more information see the docstring for the methods Quadratic,  galois, root, conjugates.

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime CyclotomicNumbers.testmat(12)^2;  # on Julia 1.7.2
  287.781 ms (3774785 allocations: 279.66 MiB)
```

The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```

testmat(12)^2 takes 0.31s in GAP3, 0.22s in GAP4


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1-L279' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.E' href='#CyclotomicNumbers.E'>#</a>
**`CyclotomicNumbers.E`** &mdash; *Function*.



`E(n,p=1)` returns the `Root1` equal to `ζₙᵖ`


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L359' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.Root1' href='#CyclotomicNumbers.Root1'>#</a>
**`CyclotomicNumbers.Root1`** &mdash; *Type*.



`Root1(c)`

`c`  should be a cyclotomic number (a  `Cyc`), or a `Real`. `Root1` returns `E(n,e)` if `c==E(n,e)`, and `nothing` if `c` is not a root of unity.

```julia-repl
julia> r=Root1(-E(9,2)-E(9,5))
Root1: ζ₉⁸

julia> order(r)
9

julia> exponent(r)
8

julia> Cyc(r) # the Zumbroich normal form is a sum of 2 roots of unity
Cyc{Int64}: -ζ₉²-ζ₉⁵

julia> Root1(-E(9,4)-E(9,5)) # nothing
```


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1149-L1170' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.conductor' href='#CyclotomicNumbers.conductor'>#</a>
**`CyclotomicNumbers.conductor`** &mdash; *Function*.



`conductor(c::Cyc)`    `conductor(a::AbstractArray)`

returns the smallest positive integer  n such that `c∈ ℚ (ζₙ)` (resp. all elements of `a` are in `ℚ (ζₙ)`).

```julia-repl
julia> conductor(E(9))
9

julia> conductor([E(3),1//2,E(4)])
12
```


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L477-L491' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.coefficients' href='#CyclotomicNumbers.coefficients'>#</a>
**`CyclotomicNumbers.coefficients`** &mdash; *Function*.



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


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L520-L539' class='documenter-source'>source</a><br>

<a id='Base.denominator' href='#Base.denominator'>#</a>
**`Base.denominator`** &mdash; *Function*.



`denominator(c::Cyc{Rational})`

returns the smallest `d` such that `d*c` has integral coefficients (thus is an algebraic integer).


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L550-L555' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.galois' href='#CyclotomicNumbers.galois'>#</a>
**`CyclotomicNumbers.galois`** &mdash; *Function*.



galois(c::Cyc,n::Int) applies to c the galois automorphism   of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.   n should be prime to conductor(c).

**Examples**

```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
Cyc{Int64}: 1-ζ₄

julia> galois(root(5),2)==-root(5)
true
```


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1081-L1093' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.conjugates' href='#CyclotomicNumbers.conjugates'>#</a>
**`CyclotomicNumbers.conjugates`** &mdash; *Function*.



`conjugates(c)`

returns the list of distinct galois conjugates of `c` (over the Rationals), starting with c

```julia-repl
julia> conjugates(1+root(5))
2-element Vector{Cyc{Int64}}:
 1+√5
 1-√5
```


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1110-L1122' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.Quadratic' href='#CyclotomicNumbers.Quadratic'>#</a>
**`CyclotomicNumbers.Quadratic`** &mdash; *Type*.



`Quadratic(c::Cyc)` 

determines  if  `c`  lives  in  a  quadratic  extension  of  `ℚ`. The call `q=Quadratic(c)`  returns a  struct `Quadratic`  with fields  `q.a`, `q.b`, `q.root`,  `q.den` representing `c` as `(q.a + q.b root(q.root))//q.den` if such a representation is possible or returns `q===nothing` otherwise.

**Examples**

```julia-repl
julia> Quadratic(E(3,2)-2E(3))
(1-3√-3)/2

julia> Quadratic(1+E(5))

```


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1214-L1230' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.root' href='#CyclotomicNumbers.root'>#</a>
**`CyclotomicNumbers.root`** &mdash; *Function*.



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


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L1308-L1326' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.Elist' href='#CyclotomicNumbers.Elist'>#</a>
**`CyclotomicNumbers.Elist`** &mdash; *Function*.



CyclotomicNumbers.Elist(n,i)  

expresses  ζₙⁱ  in  zumbroich_basis(n):  it  is  a  sum  of some ζₙʲ with   coefficients all 1 or all -1. The result is a Pair sgn=>inds where sgn is   true  if coefficients are all 1 and false otherwise, and inds is the list   of i in 0:n-1 such that ζₙⁱ occurs with a non-zero coefficient.

Should only be called with i∈ 0:n-1


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L690-L699' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.zumbroich_basis' href='#CyclotomicNumbers.zumbroich_basis'>#</a>
**`CyclotomicNumbers.zumbroich_basis`** &mdash; *Function*.



CyclotomicNumbers.zumbroich_basis(n::Int) 

returns  the Zumbroich basis of  ℚ (ζₙ) as the  vector of i in 0:n-1 such   that `ζₙⁱ` is in the basis.


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L496-L501' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.prime_residues' href='#CyclotomicNumbers.prime_residues'>#</a>
**`CyclotomicNumbers.prime_residues`** &mdash; *Function*.



`CyclotomicNumbers.prime_residues(n)` the numbers less than `n` and prime to `n` 


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L330' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.factor' href='#CyclotomicNumbers.factor'>#</a>
**`CyclotomicNumbers.factor`** &mdash; *Function*.



`CyclotomicNumbers.factor(n::Integer)` make `Primes.factor` fast for small Ints by memoizing it


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L343-L346' class='documenter-source'>source</a><br>

<a id='CyclotomicNumbers.modZ' href='#CyclotomicNumbers.modZ'>#</a>
**`CyclotomicNumbers.modZ`** &mdash; *Function*.



`modZ(x::Rational{<:Integer})` returns `x mod ℤ` as a rational in `[0,1[`


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L355' class='documenter-source'>source</a><br>

<a id='Base.numerator-Union{Tuple{Cyc{<:Union{Rational{T}, T}}}, Tuple{T}} where T<:Integer' href='#Base.numerator-Union{Tuple{Cyc{<:Union{Rational{T}, T}}}, Tuple{T}} where T<:Integer'>#</a>
**`Base.numerator`** &mdash; *Method*.



`numerator(c::Cyc{Rational})`

returns `denominator(c)*c` as a cyclotomic over the integers.


<a target='_blank' href='https://github.com/jmichel7/CyclotomicNumbers.jl/blob/860a466d74f5885448cee6e47684c0f9ba7d1e32/src/CyclotomicNumbers.jl#L558-L562' class='documenter-source'>source</a><br>


"""
This  package deals with cyclotomic numbers,  the complex numbers which are
linear  combinations  of  roots  of  unity,  usually  with with rational or
integer  coefficients; but it is possible to use any coefficients `<:Real`.
Cyclotomic numbers can be converted to `Complex` numbers.

The  cyclotomic  numbers  with  rational  coefficients  form  a  field, the
cyclotomic  field. It  is the  maximal extension  of the  rationals with an
abelian  Galois group. Its ring of  integers is the cyclotomic numbers with
integer coefficients (called cyclotomic integers).

Cyclotomic  numbers are very  important for finite  groups, since character
values of finite groups are cyclotomic integers.

This package depends only on the packages `ModuleElt` and `Primes`. It is a
port  of the  GAP implementation  of cyclotomics,  which uses a normal form
given  by writing them in the Zumbroich basis. This form allows to find the
smallest  Cyclotomic field  which contains  a given  number, and  decide in
particular  if a cyclotomic  is zero. Let  `ζₙ=exp(2im*π/n)`. The Zumbroich
basis is a particular subset of size φ(n) of 1,ζₙ,ζₙ²,…,ζₙⁿ⁻¹ which forms a
basis of ℚ (ζₙ). The reference is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

I  started  this  file  by  porting  Christian  Stump's Sage code, which is
simpler  to understand than GAP's C  code. 

As GAP does, I lower automatically numbers after each computation, that is,
reduce  them to the smallest cyclotomic  field where they belong. Currently
the  code is somewhat  slower (depending on  the operation it  has the same
speed or is slower up to 50%) than the C code in GAP but there are probably
opportunities to optimize that I missed.

What GAP does which I do not do is convert automatically a Cyclotomic which
is  rational  to  a  `Rational`,  a  `Rational`  which  is  integral  to an
`Integer`,  a  `BigInt`  which  is  small  to  an  `Int`,  etc…  This  is a
tremendously  important optimization but because of type stability in Julia
it  needs  a  new  type  of  number  to  be  added to Julia, which I am not
competent enough to try.

This package is similar (and mostly compatible) with Marek Kaluba's package
`Cyclotomics`,  whose existence I discovered after writing this package. We
discussed merging them but concluded it would be a lot of work for benefits
which are not clear currently. Some differences are:

  - I define two types in this  package: `Root1` represents a root of unity,
    and  `Cyc` a cyclotomic number. The advantage of having a separate type
    for  roots of  unity is  that computations  are very  fast for them, of
    which  I take advantage  in the package  `CycPol` for polynomials whose
    zeros are roots of unity.

  - In Kaluba's package  numbers are not  systematically lowered but only on
    demand  (like  for  printing).  this  speeds  up some computations by a
    factor  approaching 2, but it also makes some computations I have to do
    infeasible,  like the following one, which if not lowering involves too
    large fields; the answer is `-36ζ₃²`:

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
If you `develop` my package it is easy to use the strategy of lowering only
on  demand or to use alternate implementations like dense vectors or sparse
vectors (like `Cyclotomics`) --- I prepared boolean flags to choose various
implementations  in the  code. I  have currently  chosen the implementation
with `ModuleElts` and systematic lowering as giving the best results.

The main way to build a Cyclotomic number is to use the function `E(n,k=1)`
which  constructs  the  `Root1`  equal  to  `ζₙᵏ`,  and  then  make  linear
combinations of such numbers.

# Examples
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
Cyc{Rational{Int64}}: -1//2

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

The  group of  roots of  unity is  isomorphic to  ℚ /ℤ  , thus  `Root1` are
represented internally by a rational number in `[0,1[`.

```julia-repl
julia> Root1(;r=1//4) # this constructor ensures the fraction is in [0,1[
Root1: ζ₄

julia> c=E(4)*E(3) # fast computation if staying inside roots of unity
Root1: ζ₁₂⁷
```
`Root1` have the same operations as `Cyc`, but are first converted to `Cyc`
for  any  operation  other  than  `one,  isone,  *, ^, inv, conj, /, //`. A
`Root1`  can be  raised to  a `Rational`  power, which  extracts a  root if
needed.  A `Root1` can be  taken apart using `order`  and `exponent` --- if
`a`  and `b` are prime to each other,  `a` is the order of `E(a,b)` and `b`
(taken  mod `a`) is the exponent. Note  that the order is not the conductor
since `E(6)==-E(3)` has order 6 and conductor 3.

```julia-repl
julia> c=Complex{Float64}(E(3))  # convert to Complex{float} is sometimes useful
-0.4999999999999999 + 0.8660254037844387im
```
In  presence  of  a  `Cyc`  a  number  `<:Real`  or  `<:Complex{<:Real}` is
converted to a `Cyc`.

```julia-repl
julia> 0.0+E(3)
Cyc{Float64}: 1.0ζ₃

julia> E(3)+1//2
Cyc{Rational{Int64}}: √-3/2

julia> E(3)+im
Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹
```
The  function  `complex`  converts  a  `Cyc{T}`  to  a  `Complex{T}` if the
conductor is 1 or 4, to a `Complex{float(T)}` otherwise.

```julia-repl
julia> complex(E(4))
0 + 1im

julia> complex(E(3))
-0.4999999999999999 + 0.8660254037844387im
```
`Cyc`s have methods `copy, hash, ==, cmp, isless` (total order) so they can
be  keys in hashes, elements of sets,  and can be sorted. Cyclotomics which
are  integers  or  rationals  compare  correctly  to `Real`s (this does not
extend to irrational real `Cyc`s):

```julia-repl
julia> -1<Cyc(0)<1
true
```

`Cyc`s have the operations `+, -, *, /, //, inv, ^, conj, abs2, abs, image,
real, reim, isinteger, isreal, one, isone, zero, iszero, complex, adjoint`.
Cyclotomics   with  rational  or  integer  coefficients  have  the  methods
`numerator`  and  `denominator`:  a  `Cyc`  `x`  is a Cyclotomic integer if
`denominator(x)==1`   and  then  `numerator(x)`   gives  the  corresponding
`Cyc{<:Integer}`.

You  can pick apart a cyclotomic in various ways. The fastest is to use the
iterator  `pairs` which, for a cyclotomic  `a` of conductor `e` iterates on
the  pairs `(i,c)` such that  `a` has a non-zero  coefficient `c` on `ζₑⁱ`.
You  can also get  the coefficient of  `ζₑⁱ` as `a[i]`  but it is slower to
iterate  on coefficients  this way.  Finally you  can get (efficiently) the
vector of all coefficients by `coefficients`.

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

For more information see the docstring for the methods Quadratic, 
galois, root, conjugates.

Finally, a benchmark:

```benchmark
julia> function testmat(p) 
         ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
         [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
       end
testmat (generic function with 1 method)

julia> @btime CyclotomicNumbers.testmat(12)^2;  # on Julia 1.8.5
  182.929 ms (2032503 allocations: 198.22 MiB)
```
The equivalent in GAP:

```
testmat:=function(p)local ss;ss:=Combinations([0..p-1],2);
  return List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p));
end; 
```
testmat(12)^2 takes 0.31s in GAP3, 0.22s in GAP4
"""
module CyclotomicNumbers
export coefficients, root, E, Cyc, conductor, galois, Root1, Quadratic, 
       order, conjugates, modZ

#---- formatting utilities duplicated here to avoid dependency ----------
const ok="^[-+]?([^-+*/]|√-|{-)*"
const par=Regex("(\\([^()]*\\))")
const nobf=Regex("$ok(/+)?[0-9]*\$")
const nob=Regex("$ok\$")

# `c` is a coefficient. Bracket it for clarity if needed
function bracket_if_needed(c::String;allow_frac=false)
  u=c
  while(match(par,u)!=nothing) u=replace(u,par=>"") end
  e=allow_frac ? nobf : nob
  if match(e,u)!==nothing c
  else "("*c*")" 
  end
end

# `c` is a coefficient. Bracket it if needed, simplify it if ±1
function format_coefficient(c::String;allow_frac=false)
  if c=="1" ""
  elseif c=="-1" "-"
  else bracket_if_needed(c;allow_frac)
  end
end

# format exponent as unicode or TeX
function stringexp(io::IO,n::Integer)
  if isone(n) ""
  elseif get(io,:TeX,false) 
    "^"*(n in 0:9 ? string(n) : "{"*string(n)*"}")
  elseif get(io,:limit,false)
    if n<0 res=['⁻']; n=-n else res=Char[] end
    for i in reverse(digits(n)) 
      push!(res,['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'][i+1])
    end
    String(res)
  else "^"*string(n)
  end
end

# format index as unicode or TeX
function stringind(io::IO,n::Integer)
  if get(io,:TeX,false) 
    n in 0:9 ? "_"*string(n) : "_{"*string(n)*"}"
  elseif get(io,:limit,false)
    if n<0 res=['₋']; n=-n else res=Char[] end
    for i in reverse(digits(n)) push!(res,Char(0x2080+i)) end
    String(res)
  else "_"*string(n)
  end
end

#---- next function duplicated from Combinat to avoid dependency ----------
"`CyclotomicNumbers.prime_residues(n)` the numbers less than `n` and prime to `n` "
function prime_residues(n)
  if n==1 return [0] end
  pp=trues(n-1)
  for i in 2:div(n,2)
    if !pp[i] || n%i!=0 continue end
    pp[i:i:n-1].=false
  end
  (1:n-1)[pp]
end

using Primes: factor, eachfactor

#------------------------ type Root1 ----------------------------------
"""
`Root1`  is a type representing roots of unity. The internal representation
is  by a `Rational{Int}` of the form  `y//x` where the integers `y` and `x`
satisfy  `0≤y<x`, representing the  root of unity  `ζₓʸ` (where `ζₓ` is the
root  of  unity  whose  approximate  value  is  `exp(2im*π/x)`).  Efficient
constructors are `Root1(;r=y//x)` and `E`.
"""
struct Root1 <: Number # E(c,n)
  r::Rational{Int}
  global Root1_(x)=new(x)
end

"`modZ(x::Rational{<:Integer})` returns `x mod ℤ ` as a rational in `[0,1[`"
modZ(x::Rational{<:Integer})=Base.unsafe_rational(mod(numerator(x),
                                    denominator(x)),denominator(x))

"""
`E(n,p=1)`  returns the `Root1`  equal to `ζₙᵖ` (where  `ζₙ` is the root of
unity whose approximate value is `exp(2im*π/n)`)
"""
E(c,n=1)=Root1_(modZ(n//c))

Base.exponent(a::Root1)=numerator(a.r)
order(a::Root1)=denominator(a.r)
conductor(a::Root1)=order(a)%4==2 ? div(order(a),2) : order(a)

Root1(;r=0//1)=Root1_(modZ(Rational{Int}(r)))

function Root1(c::Real)
  if isone(c) Root1_(0//1)
# elseif iszero(c) Cyc(0)
  elseif c==-1 Root1_(1//2)
  else nothing
  end
end

function Root1(c::Complex)
  if isone(c) Root1_(0//1)
  elseif c==-1 Root1_(1//2)
  elseif c==im Root1_(1//4)
  elseif c==-im Root1_(3//4)
  else nothing
  end
end

Base.broadcastable(r::Root1)=Ref(r)

function Base.show(io::IO, ::MIME"text/plain", r::Root1)
  if !haskey(io,:typeinfo) print(io,"Root1: ") end
  show(io,r)
end

function Base.show(io::IO, r::Root1)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  d=exponent(r)
  c=order(r)
  if repl || TeX
    if c==1 print(io,"1")
    elseif c==2 print(io,"-1")
    else print(io,TeX ? "\\zeta" : "ζ",stringind(io,c),stringexp(io,d))
    end
  else
    print(io,(d==1 || d==0) ? "E($c)" : "E($c,$d)")
  end
end

function Base.cmp(a::Root1,b::Root1)
  r=cmp(conductor(a),conductor(b))
  if !iszero(r) return r end
  cmp(Cyc(a),Cyc(b))
end

Base.isless(a::Root1,b::Root1)=cmp(a,b)==-1
Base.isless(c::Root1,d::Real)=Cyc(c)<Cyc(d)
Base.isless(d::Real,c::Root1)=Cyc(d)<Cyc(c)
Base.:(==)(a::Root1,b::Root1)=a.r==b.r
Base.one(a::Root1)=Root1_(0//1)
Base.zero(a::Root1)=zero(Cyc{Int})
Base.isone(a::Root1)=iszero(a.r)
Base.iszero(a::Root1)=false
Base.:*(a::Root1,b::Root1)=Root1_(modZ(a.r+b.r))

Base.:^(a::Root1,n::Integer)=Root1_(modZ(n*a.r))
Base.:^(a::Root1,r::Rational{<:Integer})=root(a,denominator(r))^numerator(r)
Base.inv(a::Root1)=Root1_(modZ(-a.r))
Base.conj(a::Root1)=inv(a)
Base.:/(a::Root1,b::Root1)=a*inv(b)
Base.://(a::Root1,b::Root1)=a/b

#------------------------ type Cyc ----------------------------------
const impl=:MM # I tried 4 different implementations. 
# For testmat(12)^2
# :MM,ModuleElt is fastest
# :MM,HModuleElt is 70% slower than ModuleElt
# :svec is 60% slower than ModuleElt
# :vec is 50% slower than ModuleElt
const lazy=false # whether to lower all the time or on demand

if impl==:vec
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  d::Vector{T} # length(d)==conductor; i-th element is the coefficient on ζⁱ⁻¹
  global function Cyc_(d::Vector{T}) where T<:Real 
    new{T}(d)
  end
end
function Cyc(c::Integer,v::AbstractVector)
  if c!=length(v) error("c=$c but length(v)=$(length(v))\n") end
  Cyc_(v)
end
Base.pairs(c::Cyc)=Iterators.map(x->x[1]-1=>x[2],
                              Iterators.filter(x->x[2]!=0,enumerate(c.d)))
elseif impl==:svec
using SparseArrays
mutable struct Cyc{T <: Real}<:Number    # a cyclotomic number
  d::SparseVector{T,Int} # d[i]==coeff on ζⁱ⁻¹ (where i∈ zumbroich_basis(n))
  global function Cyc_(d::SparseVector{T}) where T<:Real 
#   for debugging you may uncomment the following line
#   if !issorted(d.nzind) || any(iszero,d.nzval) error(d) end
    new{T}(d)
  end
end
Base.pairs(c::Cyc)=Iterators.map(x->x[1]-1=>x[2],zip(c.d.nzind,c.d.nzval))
function Cyc(c::Integer,v::SparseVector)
  if c!=length(v) error("c=$c but length(v)=$(length(v))\n") end
  Cyc_(v)
end

elseif impl==:MM
using ModuleElts
const MM=ModuleElt # you can try with HModuleElt
mutable struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int
  d::MM{Int,T} # pairs: i=>coeff on ζⁱ (where i∈ zumbroich_basis(n))
end
Base.pairs(c::Cyc)=c.d
end

if impl==:MM
  conductor(c::Cyc)=c.n
  Base.getindex(c::Cyc,i::Integer)=c.d[mod(i,conductor(c))]
else
  conductor(c::Cyc)=length(c.d)
  Base.getindex(c::Cyc,i::Integer)=c.d[mod(i,conductor(c))+1]
end

Base.valtype(c::Cyc{T}) where T =T # how to recover T from c

"""
   `conductor(c::Cyc)`

returns the smallest positive integer  n such that `c∈ ℚ (ζₙ)`

   `conductor(a::AbstractArray)`

smallest positive integer  n such that all elements of `a` are in `ℚ (ζₙ)`

```julia-repl
julia> conductor(E(6))
3

julia> conductor([E(3),1//2,E(4)])
12
```
"""
conductor(a::AbstractArray)=lcm(conductor.(a))
conductor(i::Integer)=1 # for convenience
conductor(i::Rational{<:Integer})=1 # for convenience

"""
  CyclotomicNumbers.zumbroich_basis(n::Int) 

  returns  the Zumbroich basis of ℚ (ζₙ) as the sorted vector of i in 0:n-1
  such that `ζₙⁱ` is in the basis.
"""
function zumbroich_basis(n::Int)
# This function is not used in the rest of this module. We use Elist.
  if n==1 return [0] end
  function J(k::Int, p::Int) # see [Breuer] Rem. 1 p. 283
    k==0 ? (p==2 ? (0:0) : (1:p-1)) : p==2 ? (0:1) : (div(1-p,2):div(p-1,2))
  end
  sort(vec(sum.(Iterators.product((div(n,p^k)*J(k-1,p) 
       for (p,np) in eachfactor(n) for k in 1:np)...))).%n)
end

"""
`coefficients(c::Cyc)`

for  a cyclotomic `c` of conductor `n`,  returns a vector `v` of length `n`
such that `c==∑ᵢ vᵢ₊₁ ζⁱ` for `i∈ 0:n-1`.

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
"""
function coefficients(c::Cyc)
if impl==:svec return Array(c.d)
elseif impl==:vec return copy(c.d)
else
  res=zeros(valtype(c),conductor(c))
  @inbounds for (i,v) in pairs(c) res[i+1]=v end
  res
end
end
  
"""
`denominator(c::Cyc{<:Rational})`

returns   the   smallest   integer   `d`   such  that  `d*c`  has  integral
`coefficients` (thus is an algebraic integer).
"""
Base.denominator(c::Cyc)=lcm(denominator.(values(c.d)))

"""
`numerator(c::Cyc{<:Rational})`

returns `denominator(c)*c` as a cyclotomic over the integers.
"""
Base.numerator(c::Cyc{<:Union{T,Rational{T}}}) where T<:Integer=Cyc{T}(c*denominator(c))

if impl==:vec
Cyc(i::Real)=Cyc_([i])
elseif impl==:svec
Cyc(i::Real)=Cyc_(iszero(i) ? spzeros(typeof(i),1) : SparseVector(1,[1],[i]))
else
Cyc(i::Real)=Cyc(1,MM(iszero(i) ? Pair{Int,typeof(i)}[] : [0=>i];check=false))
end
Cyc{T}(i::Real) where T<:Real=Cyc(T(i))

Cyc{T}(a::Root1) where T<:Real=Cyc{T}(Cyc(a))

Base.zero(c::Cyc{T}) where T=Cyc{T}(0)

if impl==:svec obviouslyzero(c::Cyc)=nnz(c.d)==0 # faster than iszero(c.d)
else           obviouslyzero(c::Cyc)=iszero(c.d)
end
if lazy Base.iszero(c::Cyc)=obviouslyzero(lower!(c))
        Base.isone(c::Cyc)=lower!(c);isone(conductor(c)) && isone(num(c))
else    Base.iszero(c::Cyc)=obviouslyzero(c)
        Base.isone(c::Cyc)=isone(conductor(c)) && isone(num(c))
end
Base.zero(::Type{Cyc{T}}) where T=Cyc{T}(0)
Base.one(c::Cyc{T}) where T =Cyc{T}(1)

function Cyc(c::Complex)
  if iszero(imag(c)) return Cyc(real(c)) end
if impl==:vec
  Cyc(4,[real(c),imag(c),zero(real(c)),zero(real(c))])
elseif impl==:svec
  Cyc_(dropzeros!(SparseVector(4,[1,2],[real(c),imag(c)])))
else
  if iszero(real(c)) Cyc(4,MM(1=>imag(c);check=false))
  else Cyc(4,MM(0=>real(c),1=>imag(c);check=false))
  end
end
end

Cyc{T}(c::Complex) where T=T(real(c))+E(4)*T(imag(c))

function Cyc{T}(c::Cyc{T1}) where {T,T1}
  if T==T1 return c end
  if impl==:MM Cyc(conductor(c),convert(MM{Int,T},c.d))
  else Cyc(conductor(c),T.(c.d))
  end
end

# num(c): value of c as a real when conductor(c)==1
if impl==:MM
  num(c::Cyc{T}) where T =iszero(c) ? zero(T) : last(first(c.d))
else
  num(c::Cyc)=c.d[1]
end

function (::Type{T})(c::Cyc)where T<:Union{Integer,Rational}
  @static if lazy lower!(c) end
  if conductor(c)==1 return T(num(c)) end
  throw(InexactError(:convert,T,c))
end

function Base.convert(::Type{T},c::Cyc;check=true)where T<:AbstractFloat
  if check && !isreal(c) throw(InexactError(:convert,T,c)) end
  real(convert(Complex{T},c))
end

function (::Type{T})(c::Cyc)where T<:AbstractFloat 
  convert(T,c;check=true)
end

(::Type{T})(a::Root1) where T<:Number = T(Cyc(a))

function Complex{T}(c::Cyc)where T<:AbstractFloat
  sum(((k,v),)->T(v)*cispi(2*T(k)/conductor(c)),pairs(c);init=Complex{T}(0.0))
end

function Complex{T}(c::Cyc)where T<:Union{Integer,Rational}
  @static if lazy lower!(c) end
  if conductor(c)==1 return Complex{T}(num(c)) end
  if conductor(c)==4 
    res=Complex{T}(0)
    for (k,v) in pairs(c)
      res+=k==0 ? v : im*v
    end
    return res
  end
  throw(InexactError(:convert,Complex{T},c))
end

Complex{T}(a::Root1) where T=Complex{T}(Cyc(a))
function Base.complex(c::Cyc{T}) where T 
  @static if lazy lower!(c) end
  (conductor(c)==1||conductor(c)==4) ? Complex{T}(c) : Complex{float(T)}(c)
end

Base.complex(a::Root1)=complex(Cyc(a))

function Base.isinteger(c::Cyc)
  @static if lazy lower!(c) end
  conductor(c)==1 && isinteger(num(c))
end

Base.isinteger(a::Root1)=a.r==0//1 || a.r==1//2

Base.isreal(c::Cyc)=conductor(c)==1 || c==conj(c)
Base.isreal(a::Root1)=isinteger(a)

function Base.real(c::Cyc)
  @static if lazy lower!(c) end
  if conductor(c)==1 return num(c) end
  (c+conj(c))*1//2
end

Base.real(a::Root1)=real(Cyc(a))

function Base.imag(c::Cyc)
  @static if lazy lower!(c) end
  if conductor(c)==1 return 0 end
  E(4)*(conj(c)-c)*1//2
end

Base.imag(a::Root1)=imag(Cyc(a))
Base.reim(c::Cyc)=(real(c),imag(c))
Base.abs2(a::Root1)=abs2(Cyc(a))
Base.abs(a::Root1)=one(a)

# memoize Elist
const Elist_dict=Dict{Tuple{Int,Int},Pair{Bool,Vector{Int}}}((1,0)=>(true=>[0])) 
"""
`CyclotomicNumbers.Elist(n,i)`

expresses  `ζₙⁱ` in  `zumbroich_basis(n)`: it  is a  sum of some `ζₙʲ` with
coefficients  all 1  or all  -1. The  result is  a `Pair` `sgn=>inds` where
`sgn` is `true` if coefficients are all 1 and `false` otherwise, and `inds`
is  the  list  of  `i`  in  `0:n-1`  such that `ζₙⁱ` occurs with a non-zero
coefficient.

Should only be called with `i∈ 0:n-1`
"""
function Elist(n::Int,i::Int)
  get!(Elist_dict,(n,i)) do
    mp=Int[]
    j=i
    for (p,np) in eachfactor(n)
      f=p^np
      m=div(n,f)
      cnt=mod(j*invmod(m,f),f)
      j-=cnt*m
      if p==2
        if 1==div(cnt,div(f,p)) push!(mp,p) end
      else
        tmp=zeros(Int,np)
        for k in 1:np
          f=div(f,p)
          tmp[k],cnt=divrem(cnt,f)
        end
        for k in np-1:-1:1
          if tmp[k+1]>div(p-1,2)
            tmp[k]+=1
            if k==1 && tmp[k]==p tmp[k]=0 end
          end
        end
        if tmp[1]==0 push!(mp,p) end
      end
    end
    if isempty(mp) return true=>[i] end
    v=vec(sum.(Iterators.product((div(n,p)*(1:p-1) for p in mp)...)))
    v=sort((i.+v).%n)
    iseven(length(mp))=>v
  end
end

# p iterator of pairs i=>c meaning c*E(n,i)
# This constructor is not in API since the result may or may not need lowering
if true
function Cyc(n::Integer,::Type{T},p) where T 
  res=if     impl==:vec  zeros(T,n)
      elseif impl==:svec spzeros(T,n)
      elseif impl==:MM   Pair{Int,T}[]
      end
  for (i,c) in p
    (s,v)=Elist(n,mod(i,n))
    if !s c=-c end
    if impl==:MM for k in v push!(res,k=>c) end
    else @inbounds view(res,v.+1).+=c
    end
  end
  if impl==:MM Cyc(n,MM(res))
  elseif impl==:svec Cyc(n,dropzeros!(res))
  else Cyc(n,res)
  end
end
else # attempt to be revisited
function Cyc(n::Integer,::Type{T},p) where T 
  if impl==:MM && length(p)==1
    (i,c)=first(p)
    (s,v)=Elist(n,mod(i,n))
    if !s c=-c end
    return Cyc(n,MM(v.=>c;check=iszero(c)))
  end
  res=if     impl==:vec  zeros(T,n)
      elseif impl==:svec spzeros(T,n)
      elseif impl==:MM 
        (0:n-1).=>T(0)
      end
  for (i,c) in p
    (s,v)=Elist(n,mod(i,n))
    if !s c=-c end
    if impl==:MM @inbounds for k in v res[k+1]=k=>last(res[k+1])+c end
    else @inbounds view(res,v.+1).+=c
    end
  end
  if impl==:MM Cyc(n,MM(filter!(r->!iszero(last(r)),res);check=false))
  elseif impl==:svec Cyc(n,dropzeros!(res))
  else Cyc(n,res)
  end
end
end

function Cyc(a::Root1) # the result is guaranteed lowered
  n=order(a)
  e=exponent(a)
  n%4==2 ? Cyc(div(n,2),Int,(div(e,2)+div(n+2,4)=>-1,)) : Cyc(n,Int,(e=>1,))
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:Real}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:Real}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Cyc{T2}})where {T1,T2}
  Cyc{promote_type(T1,T2)}
end

# total order is necessary to put Cycs in a sorted list
# for conductor(c)==1  a<b is as expected
function Base.cmp(a::Cyc,b::Cyc)
  @static if lazy lower!(a);lower!(b) end
  t=cmp(conductor(a),conductor(b))
  if !iszero(t) return t end
  conductor(a)==1 ? cmp(num(a),num(b)) : cmp(a.d,b.d)
end

if lazy
function Base.:(==)(a::Cyc,b::Cyc)
  if conductor(a)==conductor(b) return a.d==b.d end
  lower!(a);lower!(b)
  conductor(a)==conductor(b) && a.d==b.d
end
else
Base.:(==)(a::Cyc,b::Cyc)=conductor(a)==conductor(b) && a.d==b.d
end
Base.isless(a::Cyc,b::Cyc)=cmp(a,b)==-1
Base.isless(c::Cyc,d::Real)=c<Cyc(d)
Base.isless(d::Real,c::Cyc)=Cyc(d)<c
Base.isless(a::Root1,b::Cyc)=Cyc(a)<b
Base.isless(a::Cyc,b::Root1)=a<Cyc(b)

# hash is necessary to put Cycs as keys of a Dict or make a Set
function Base.hash(a::Cyc, h::UInt)
  @static if lazy lower!(a) end
  hash(a.d, hash(conductor(a), h))
end

function Base.show(io::IO, ::MIME"text/html", a::Cyc)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", r::Cyc)
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function normal_show(io::IO,p::Cyc)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if valtype(p)<:Rational{<:Integer}
    den=denominator(p)
    p=Cyc{typeof(den)}(p*den)
  else
    den=1
  end
  res=""
  for (deg,v) in pairs(p)
    if deg==0 t=string(v)
    else 
      t=format_coefficient(string(v))
      if repl || TeX
        r=(TeX ? "\\zeta" : "ζ")*stringind(io,conductor(p))*stringexp(io,deg)
      else
        r=(deg==1 ? "E($(conductor(p)))" : "E($(conductor(p)),$deg)")
      end
      t*=r
    end
    if t[1]!='-' t="+"*t end
    res*=t
  end
  if res[1]=='+' res=res[2:end] end
  if !isone(den) 
    res=bracket_if_needed(res)
    res*=repl||TeX ? "/$den" : "//$den"
  end
  res
end

function Base.show(io::IO, p::Cyc)
  quadratic=get(io,:quadratic,true)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  function printtyped(n)
    if repl||TeX|| get(io,:typeinfo,Any)==typeof(p) print(io,n)
    else print(io,typeof(p),"(",n,")")
    end
  end
  @static if lazy lower!(p) end
  if conductor(p)==1
    n=num(p)
    if n isa Integer || n isa Rational{<:Integer}
      if denominator(n)==1 printtyped(numerator(n))
      else printtyped(n)
      end
      return
    else
      printtyped(n)
      return
    end
  end
  rqq=[normal_show(io,p)]
  if quadratic && (valtype(p)<:Integer || valtype(p)<:Rational{<:Integer})
    q=Quadratic(p)
    if !isnothing(q) pushfirst!(rqq,repr(q;context=io)) end
    for test in [1-E(4),1+E(4),E(3),E(3,2),1-E(3),1-E(3,2),1+E(3),1+E(3,2),root(-3)]
      if !iszero(conductor(p)%conductor(test)) continue end
      q=Quadratic(p*1//test)
      if isnothing(q) continue end
      rq=repr(q;context=io)
      rq=format_coefficient(rq;allow_frac=true)
      t=format_coefficient(normal_show(io,test))
      if !(repl||TeX) t*="*" end
      if !isempty(rq) && rq[1]=='-' rq="-"*t*rq[2:end] else rq=t*rq end
      push!(rqq,rq)
    end
  end
  n=rqq[argmin(length.(rqq))]
  if repl||TeX|| get(io,:typeinfo,Any)==typeof(p) print(io,n)
  else print(io,typeof(p),"(",n,")")
  end
end

# write a,b in common field Q(ζ_lcm(conductor(a),conductor(b)))
function promote_conductor(a::Cyc{T},b::Cyc{T})where T
  if conductor(a)==conductor(b) return (a, b) end
  n=lcm(conductor(a),conductor(b))
  if n!=conductor(a)
    m=div(n,conductor(a))
    let m=m
      a=Cyc(n,valtype(a),i*m=>v for (i,v) in pairs(a))
    end
  end
  if n!=conductor(b)
    m=div(n,conductor(b))
    let m=m # the eltype of next iterator is Any!
      b=Cyc(n,valtype(b),i*m=>v for (i,v) in pairs(b))
    end
  end
  (a,b)
end

function Base.:+(x::Cyc,y::Cyc)
  a,b=promote(x,y)
  if obviouslyzero(a) return b end
  if obviouslyzero(b) return a end
  a,b=promote_conductor(a,b)
if impl==:svec res=Cyc_(dropzeros!(a.d+b.d))
else           res=Cyc(conductor(a),a.d+b.d)
end
  @static if !lazy lower!(res) end
  res
end

Base.:-(a::Cyc)=Cyc(conductor(a),-a.d)

function Base.:-(x::Cyc,y::Cyc)
  a,b=promote(x,y)
  if obviouslyzero(a) return -b end
  if obviouslyzero(b) return a end
  a,b=promote_conductor(a,b)
if impl==:svec res=Cyc_(dropzeros!(a.d-b.d))
else           res=Cyc(conductor(a),a.d-b.d)
end
  @static if !lazy lower!(res) end
  res
end

function Base.div(c::Cyc,a::Real)
if impl==:MM
  n=merge(div,c.d,a)
  Cyc(iszero(n) ? 1 : conductor(c),n)
else
  res=div.(c.d,a)
  iszero(res) ? Cyc{eltype(res)}(0) : Cyc(conductor(c),res)
end
end

function Base.:*(c::Cyc,a::Real)
  # the check for one saves a lot of time in some applications
  if isone(a) return Cyc{promote_type(valtype(c),typeof(a))}(c) end
  if iszero(a) return zero(Cyc{promote_type(valtype(c),typeof(a))}) end
  Cyc(conductor(c),c.d*a)
end
Base.:*(a::Real,c::Cyc)=c*a

function Base.:*(a::Cyc,b::Cyc;reduce=!lazy)
  a,b=promote(a,b)
  if obviouslyzero(a) return a end
  if obviouslyzero(b) return b end
  if conductor(a)==1 return b*num(a)
  elseif conductor(b)==1 return a*num(b)
  end
  n=lcm(conductor(a),conductor(b))
  na=div(n,conductor(a))
  nb=div(n,conductor(b))
  c=Cyc(n,valtype(a),na*i+nb*j=>ai*bj for (i,ai) in pairs(a),(j,bj) in pairs(b))
  reduce ? lower!(c) : c
end

Base.://(c::Cyc,a::Real)=Cyc(conductor(c),impl==:MM ? c.d//a : c.d.//a)
Base.://(a::Union{Cyc,Real,Root1},c::Cyc)=a*inv(
                               Cyc{promote_type(valtype(c),Rational{Int})}(c))
Base.://(c::Root1,a::Real)=Cyc(c)//a
Base.://(a::Union{Real,Cyc},c::Root1)=a//Cyc(c)
Base.:/(c::Cyc,a::Real)=c*inv(a)
Base.:/(a::Cyc,c::Cyc)=a*inv(c)
Base.div(c::Root1,a)=div(Cyc(c),a)

# change c to have data n,v
function Cyc!(c,n,v) 
  if impl==:svec
    c.d=v
  elseif impl==:vec
    resize!(c.d,length(v))
    c.d.=v
  else
    c.n=n
    if MM==ModuleElt 
      resize!(c.d.d,length(v))
      c.d.d.=v
    else 
      dv=Dict(v)
      empty!(c.d.d)
      merge!(c.d.d,dv)
    end
  end
  c
end

function lower!(c::Cyc) # write c in smallest Q(ζ_n) where it sits
  n=conductor(c)
#  print("lowering conductor=",conductor(c));@show c.d
  if n==1 return c end
  if obviouslyzero(c)
    return Cyc!(c,1,impl==:MM ? empty(c.d.d) : impl==:vec ? zeros(valtype(c),1)
                 : spzeros(valtype(c),1))
  end
  for (p,np) in eachfactor(n)
    m=div(n,p)
if impl==:vec
    kk=Iterators.filter(i->c.d[i]!=0,eachindex(c.d))
    if np>1 
      if all(k->(k-1)%p==0,kk) 
        v=zeros(valtype(c),m)
        for x in kk v[1+div(x-1,p)]=c.d[x] end
        return lower!(Cyc!(c,m,v))
      end
    elseif iszero(count(!iszero,c.d)%(p-1))
      u=zeros(Int,m)
      for k in kk u[1+((k-1)%m)]+=1 end
      if !all(x->iszero(x) || x==p-1,u) continue end
      i=0; for j in eachindex(u) if !iszero(u[j]) i+=1;u[i]=j-1 end end
      resize!(u,i)
      @. u=div(u+m*mod(-u,p)*invmod(m,p),p)%m
      v=zeros(valtype(c),m)
      if p==2  
        @views @. v[u+1].=c.d[1+(u*p)%n]
        return lower!(Cyc!(c,m,v))
      elseif all(k->allequal(c.d[1+(m*i+k*p)%n] for i in 1:p-1),u)
        @views @. v[u+1].=-c.d[1+(m+u*p)%n]
        return lower!(Cyc!(c,m,v))
      end
    end
elseif impl==:svec
    kk=c.d.nzind
    val=c.d.nzval
    if np>1 
      if all(k->(k-1)%p==0,kk) 
        return lower!(Cyc!(c,m,SparseVector(m,map(x->1+div(x-1,p),kk),val)))
      end
    elseif iszero(length(kk)%(p-1))
      u=zeros(Int,m)
      for k in kk u[1+(k-1)%m]+=1 end
      if !all(x->iszero(x) || x==p-1,u) continue end
      i=0; for j in eachindex(u) if !iszero(u[j]) i+=1;u[i]=j-1 end end
      resize!(u,i)
      @. u=div(u+m*mod(-u,p)*invmod(m,p),p)%m
      if !issorted(u) sort!(u) end
      if p==2  
        return lower!(Cyc!(c,m,SparseVector(m,u.+1,[c.d[1+(k*p)%n] for k in u])))
      elseif all(k->allequal(c.d[1+(m*i+k*p)%n] for i in 1:p-1),u)
        return lower!(Cyc!(c,m,SparseVector(m,u.+1,[-c.d[1+(m+k*p)%n] for k in u])))
      end
    end
elseif impl==:MM
    if np>1 
      if all(k->first(k)%p==0,c.d) 
        return lower!(Cyc!(c,m,div(k,p)=>v for (k,v) in c.d))
      end
    elseif iszero(length(c.d)%(p-1))
      u=zeros(Int,m)
      for (k,v) in c.d u[1+(k%m)]+=1 end
      if !all(x->iszero(x) || x==p-1,u) continue end
      i=0; for j in eachindex(u) if !iszero(u[j]) i+=1;u[i]=j-1 end end
      resize!(u,i)
      @. u=div(u+m*mod(-u,p)*invmod(m,p),p)%m
      if p==2  
        v=[k=>c.d[(k*p)%n] for k in u];if !issorted(v) sort!(v) end
        return lower!(Cyc!(c,m,v))
      elseif all(k->allequal(c.d[(m*i+k*p)%n] for i in 1:p-1),u)
        v=[k=>-c.d[(m+k*p)%n] for k in u];if !issorted(v) sort!(v) end
        return lower!(Cyc!(c,m,v))
      end
    end
end
  end
  c
end

galois(c::Rational,n::Integer)=c

galois(c::Integer,n::Integer)=c

"""
`galois(c::Cyc,n::Integer)`  applies  to  `c`  the  galois  automorphism  of `ℚ
(ζ_conductor(c))`  raising  all  roots  of  unity  to the `n`-th power. `n`
should be prime to `conductor(c)`.
# Examples
```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
Cyc{Int64}: 1-ζ₄

julia> galois(root(5),2)==-root(5)
true
```
"""
function galois(c::Cyc,n::Integer)
  if gcd(n,conductor(c))!=1 
    throw(DomainError(n,"should be prime to conductor $(conductor(c))"))
  end
  if obviouslyzero(c) return c end
  Cyc(conductor(c),valtype(c),e*n=>p for (e,p) in pairs(c))
end

function galois(c::Root1,n::Integer)
  d=order(c)
  if gcd(n,d)!=1 error("$n should be prime to order($c)=$d") end
  Root1(;r=n*c.r)
end

Base.conj(c::Cyc)=galois(c,-1)

"""
`conjugates(c)`

returns the list of distinct galois conjugates of `c` (over the Rationals),
starting with `c`

```julia-repl
julia> conjugates(1+root(5))
2-element Vector{Cyc{Int64}}:
 1+√5
 1-√5
```
"""
function conjugates(c::Union{Cyc,Root1}) # Root1 or Cyc
  res=[c]
  for i in prime_residues(c isa Cyc ? conductor(c) : order(c))[2:end]
    c1=galois(c,i)
    if !(c1 in res) push!(res,c1) end
  end
  res
end

conjugates(c::Union{Rational{<:Integer},Integer})=[c]

function Base.inv(c::Cyc{<:Integer})
  if conductor(c)==1 return Cyc(inv(num(c))) end
  l=conjugates(c)
  r=l[2]
  for t in l[3:end] r=*(r,t;reduce=false) end
  a=num(*(c,r;reduce=true))
# if isone(a) return r end
# T<:Integer ? r//a : r/a
  r//a
end

function Base.inv(c::Cyc)
  if conductor(c)==1 return Cyc(inv(num(c))) end
  l=conjugates(c)
  r=l[2]
  for t in l[3:end] r=*(r,t;reduce=false) end
  a=num(*(c,r;reduce=true))
  r/a
end

Base.:^(a::Cyc, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                   Base.power_by_squaring(inv(a),-n)

Base.abs2(c::Cyc)=c*conj(c)
Base.abs(c::Cyc)=conductor(c)==1 ? Cyc(abs(num(c))) : abs(complex(c))
Base.adjoint(c::Cyc)=conj(c)

"""
`Root1(c)`
    
`c`  should  be  a  `Cyc`,  a  `Real`,  or  a `Complex`. `Root1(c)` returns
`E(n,e)`  if `c==E(n,e)`, and  `nothing` if `c`  is not equal  to a root of
unity.

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
""" 
function Root1(c::Cyc)
  @static if lazy lower!(c) end
  if !(all(x->last(x)==1,pairs(c)) || all(x->last(x)==-1,pairs(c)))
    return nothing
  end
  n=conductor(c)
  for i in prime_residues(n)
    if c==E(n,i) return Root1_(i//n) end
    if -c==E(n,i)
      if iseven(n) return E(n,div(n,2)+i)
      else return E(2n,n+2*i)
      end
    end
  end
  # return nothing
end

Base.:(==)(a::Root1,b::Number)=Cyc(a)==b # too expensive in lazy case
Base.:+(a::Root1,b::Root1)=Cyc(a)+Cyc(b)
Base.:-(a::Root1,b::Root1)=Cyc(a)-Cyc(b)
Base.:-(r::Root1)=-Cyc(r)
Base.promote_rule(a::Type{Root1},b::Type{Cyc{T}}) where T =b
Base.promote_rule(a::Type{Root1},b::Type{<:Real})=Cyc{b}
Base.promote_rule(a::Type{Root1},b::Type{Complex{T}}) where T =Cyc{promote_type(T,Int)}

struct Quadratic
  a
  b
  root
  den
end

"""
  `Quadratic(c::Cyc)` 
  
determines  if  `c`  lives  in  a  quadratic  extension  of  `ℚ `. The call
`q=Quadratic(c)`  returns a  struct `Quadratic`  with fields  `q.a`, `q.b`,
`q.root`,  `q.den` such that `c==(q.a + q.b root(q.root))//q.den` if such a
representation is possible or returns `nothing` otherwise.

# Examples
```julia-repl
julia> Quadratic(E(3,2)-2E(3))
(1-3√-3)/2

julia> Quadratic(1+E(5))

```
"""
function Quadratic(c::Cyc)
  den=denominator(c)
  c=numerator(c)
  if conductor(c)==1 return Quadratic(num(c),0,1,den) end
  f=factor(conductor(c))
  v2=get(f,2,0)
  if v2>3 || (v2==2 && any(p->p[1]!=2 && p[2]!=1,f)) ||
     (v2<2 && any(x->x!=1,values(f)))
    return nothing
  end
  f=keys(f)
  if v2==0
    sqr=conductor(c)
    if sqr%4==3 sqr=-sqr end
    gal=conjugates(c)
    if length(gal)!=2 return nothing end
    a=numerator(convert(valtype(c),sum(gal)))  # trace of 'c' over the rationals
    if iseven(length(f)) b=2*c[1]-a
    else b=2*c[1]+a
    end
    if iseven(a) && iseven(b) a>>=1; b>>=1; d=1
    else d=2
    end
  elseif v2==2
    sqr=conductor(c)>>2
    if sqr==1 a=c[0];b=-c[1]
    else
      a=c[4]
      if iseven(length(f)) a=-a end
      b=-c[sqr+4]
    end
    if sqr%4==1 sqr=-sqr; b=-b end
    d=1
  else		# v2 = 1 or 3
    sqr=conductor(c)>>2
    if sqr==2
      a=c[0];b=c[1]
      if b==c[3] sqr=-2 end
    else
      a=c[8]
      if iseven(length(f)) a=-a end
      b=c[(sqr>>1)+8]
      if b!=-c[3*(sqr>>1)-8] sqr=-sqr
      elseif (sqr>>1)%4==3 b=-b
      end
    end
    d=1
  end
  if d*c!=a+b*root(sqr) return nothing end
  return Quadratic(a,b,sqr,den*d)
end

function Base.show(io::IO,q::Quadratic)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  rq=string(q.a)
  if q.b!=0 
    if iszero(q.a) rq=""
    elseif q.b>0 rq*="+" end
    rq*=q.b==1 ? "" : q.b==-1 ? "-" : string(q.b)
    r=string(q.root)
    rq*=TeX ? "\\sqrt{$r}" : repl ? "√$r" : "root($r)"
    if !iszero(q.a) && q.den!=1 rq="("*rq*")" end
  end
  print(io,rq)
  if q.den!=1 && rq!="0" print(io,(repl||TeX) ? "/" : "//",q.den) end
end

const inforoot=Ref(false)
xprint(x...;p...)=print(IOContext(stdout,:limit=>true,p...),x...)
function proot(x,n,r)
  xprint("root(",x;quadratic=false)
  if n!=2 xprint(",",n) end
  xprint(")=",r,"\n";quadratic=false)
end
const Irootdict=Dict{Tuple{Int,Int},Union{Cyc{Int},Int,Cyc{BigInt}}}()

"""
`root(x,n=2)`

computes  an `n`-th root of `x`  when we know how to  do it. We know how to
compute  `n`-th roots  of roots  of unity,  and `n`-th  or `2n`-th roots of
perfect `n`-th powers of integers or rationals.

```julia-repl
julia> root(-1)
Cyc{Int64}: ζ₄

julia> root(E(4))
Root1: ζ₈

julia> root(27//8,6)
Cyc{Rational{Int64}}: √6/2
```
"""
function root(x::Union{Int,BigInt},n=2) # only defined for Int so no conflict with PuiseuxPols
  if isone(n) || (!lazy && isone(x)) return x end
  if n==2 && x>=0 
    r=isqrt(x);if r^2==x return r end 
  end
  get!(Irootdict,(n,x)) do
    if x==1 || (x==-1 && isodd(n)) return x end
    if x<0 && n==2 return E(4)*root(-x) end
    l=factor(x)
    if any(y->(2y)%n!=0,values(l)) 
      if x==-1 return root(E(2),n) end
      error("root($x,$n) not implemented")
    end
    a=prod(p^div(pow,n) for (p,pow) in l)
    b=[p for (p,pow) in l if pow%n!=0]
    for p in b
      if p%4==1 a*=sum(k->E(p,k^2),1:p)
      elseif p==2 a*=E(8)-E(8,3)
      elseif p%4==3 a*=E(4,3)*sum(k->E(p,k^2),1:p)
      end
    end
    if inforoot[] proot(x,n,a) end
    a
  end
end

root(x::Rational,n=2)=root(numerator(x),n)//root(denominator(x),n)

# find the "canonical" best of the n possible roots
function root(r::Root1,n=2)
  if iszero(n) return one(r) end
  d=order(r)
  j=1
  n1=n
  while true
    k=gcd(n1,d)
    n1=div(n1,k)
    j*=k
    if k==1 break end
  end
  E(j*d,exponent(r)*gcdx(n1,d)[2])
end

const Crootdict=Dict{Tuple{Int,Cyc},Cyc}()
function root(x::Cyc,n=2)
  if isone(n) || isone(x) || iszero(x) return x end
  get!(Crootdict,(n,x)) do
  d=denominator(x)
  if d!=1 return root(numerator(x),n)//root(d,n) end
  d=gcd(coefficients(x))
  if d!=1 return root(div(x,d),n)*root(d,n) end
    r=Root1(x)
    if isnothing(r) 
      if conductor(x)>1 error("cannot compute root($x,$n)") end
      return root(num(x),n)
    end
    res=Cyc(root(r,n))
    if inforoot[] proot(x,n,res) end
    res
  end
end

Base.gcd(v::Vector{<:Cyc})=reduce(gcd,v;init=zero(Cyc))
function Base.gcd(a::Cyc,b::Cyc)
  if isone(denominator(a//b))  b
  elseif isone(denominator(b//a)) a
  else gcd(gcd(collect(values(a.d))),gcd(collect(values(b.d))))
  end
end
Base.gcd(a::Cyc,b::Number)=gcd(gcd(collect(values(a.d))),b)
Base.gcd(b::Number,a::Cyc)=gcd(gcd(collect(values(a.d))),b)

# testmat(12)^2
# 1.5.3 347.534 ms (4367402 allocations: 366.17 MiB)
# 1.5.3 565.431 ms (5861810 allocations: 775.28 MiB) HModuleElts
# 1.8.5 vec 265.328 ms (2111135 allocations: 234.22 MiB)
# 1.8.5 svec 285.653 ms (3568605 allocations: 265.86 MiB)
# 1.9.0 174.725 ms (1856571 allocations: 188.89 MiB)
function testmat(p)
  ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
end

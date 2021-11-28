"""
This  package deals with cyclotomic numbers,  the complex numbers which are
linear  combinations  of  roots  of  unity  with  rational coefficients. It
depends on the packages `ModuleElt` and `Primes`.

The  cyclotomic numbers form a field, the  cyclotomic field. It is also the
maximal  extension of the rationals which  has an abelian Galois group. Its
ring of integers are the sums of unity with integral coefficients.

Cyclotomics are very important for finite groups, since character values of
finite groups are cyclotomic integers.

This package is a port of the GAP implementation of cyclotomics, which uses
a  normal form  given by  writing them  in the  Zumbroich basis.  This form
allows to find the smallest Cyclotomic field which contains a given number,
and  decide in particular  if a cyclotomic  is zero. Let ζₙ=exp(2iπ/n). The
Zumbroich  basis is  a particular  subset of  size φ(n) of 1,ζₙ,ζₙ²,…,ζₙⁿ⁻¹
which forms a basis of ℚ (ζₙ).

I  started  this  file  by  porting  Christian  Stump's Sage code, which is
simpler to understand than GAP's code. The reference for the algorithms is

T. Breuer, Integral bases for subfields of cyclotomic fields AAECC 8 (1997)

As  does  GAP,  I  lower  automatically  numbers  after  each  computation;
currently  the code is  somewhat slower (depending  on the operation it has
the same speed or is slower up to 50%) than the C code in GAP but there are
probably many opportunities to optimize that I missed.

What GAP does which I do not do is convert automatically a Cyclotomic which
is  rational to a Rational,  a Rational which is  integral to an Integer, a
BigInt  which is  small to  an Int,  etc… This  is a tremendously important
optimization  but because of type stability in Julia it needs a new type of
number to be added to Julia, which I am not competent enough to try.

We  define two types in  this package: `Root1` represents  a root of unity,
and `Cyc` a cyclotomic number. The main way to build a Cyclotomic number is
to  use the function `E(n,k=1)` which  constructs the `Root1` `ζₙᵏ`, and to
make   linear  combinations  of  such  numbers  with  integer  or  rational
coefficients.

# Examples
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

The  group of  roots of  unity is  isomorphic to  ℚ /ℤ  , thus  `Root1` are
represented internally by a rational number in `[0,1[`.

```julia-repl
julia> Root1(;r=1//4) # this contructor ensures the fraction is in [0,1[
Root1: ζ₄

julia> c=E(4)*E(3) # faster computation if staying indside roots of unity
Root1: ζ₁₂⁷

julia> c=Complex{Float64}(E(3))  # convert to Complex{float} is sometimes useful
-0.4999999999999999 + 0.8660254037844387im
```

In  presence of a float of type  `T`, a `Cyc` is converted to `Complex{T}`.
In  presence of a `Cyc`,  an integer, a rational,  or a complex of these is
converted to a `Cyc`.

```julia-repl
julia> 0.0+E(3)
-0.4999999999999999 + 0.8660254037844387im

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
be  keys in hashes or  elements of sets. Cyclotomics  which are integers or
rationals  compare  correctly  to  `Real`s  (contrary  to  irrational  real
`Cyc`s):

```julia-repl
julia> -1<Cyc(0)<1
true
```

You  can pick apart a cyclotomic in various ways. The fastest is to use the
iterator  `pairs` which, for a cyclotomic  `a` of conductor `e` iterates on
the  pairs `(i,c)` such that  `a` has a non-zero  coefficient `c` on `ζₑⁱ`.
You  can also get the coefficient `ζₑⁱ` by indexing `a[i]` but it is slower
than  `pairs` to iterate on coefficients this  way. Finally you can get the
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
"""
module Cyclotomics
export coefficients, root, E, Cyc, conductor, galois, Root1, Quadratic

#---- formatting utilities duplicated here to avoid dependency ----------
const ok="([^-+*/]|√-|{-)*"
const par="(\\([^()]*\\))"
const nobf=Regex("^[-+]?$ok$par*$ok(/+)?[0-9]*\$")
const nob=Regex("^[-+]?$ok$par*$ok\$")

function bracket_if_needed(c::String;allow_frac=false)
  e=allow_frac ? nobf : nob
  if match(e,c)!==nothing c
  else "("*c*")" 
  end
end

function format_coefficient(c::String;allow_frac=false)
  if c=="1" ""
  elseif c=="-1" "-"
  else bracket_if_needed(c;allow_frac)
  end
end

const supvec=['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹']
function stringexp(io::IO,n::Integer)
  if isone(n) ""
  elseif get(io,:TeX,false) 
    "^"*(n in 0:9 ? string(n) : "{"*string(n)*"}")
  elseif get(io,:limit,false)
    res=Char[]
    if n<0 push!(res,'⁻'); n=-n end
    for i in reverse(digits(n)) push!(res,supvec[i+1]) end
    String(res)
  else "^"*string(n)
  end
end

const subvec=['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉']
function stringind(io::IO,n::Integer)
  if get(io,:TeX,false) 
    n in 0:9 ? "_"*string(n) : "_{"*string(n)*"}"
  elseif get(io,:limit,false)
    res=Char[]
    if n<0 push!(res,'₋'); n=-n end
    for i in reverse(digits(n)) push!(res,subvec[i+1]) end
    String(res)
  else "_"*string(n)
  end
end

#---- number theory utilities duplicated here to avoid dependency ----------
" the numbers less than n and prime to n "
function prime_residues(n)
  if n==1 return [0] end
  filter(i->gcd(n,i)==1,1:n-1) # inefficient; some sieving would be better
end

import Primes
const dict_factor=Dict(2=>Primes.factor(Dict,2))
"""
`factor(n::Integer)`
make `Primes.factor` fast for small Ints by memoizing it
"""
factor(n::Integer)=get!(()->Primes.factor(Dict,n),dict_factor,n)

#---- utility duplicated here to avoid dependency ----------
constant(a)=isempty(a) || all(==(first(a)),a)

#------------------------ type Root1 ----------------------------------
struct Root1 <: Number # E(c,n)
  r::Rational{Int}
  global Root1_(x)=new(x)
end

" `E(n,p=1)` makes the `Root1` equal to `ζₙᵖ`"
E(c,n=1)=Root1_(mod(Int(n),Int(c))//Int(c))

conductor(a::Root1)=denominator(a.r)
Base.exponent(a::Root1)=numerator(a.r)

Root1(;r=0)=E(denominator(r),numerator(r)) # does a mod1

function Root1(c::Real)
  if c==1 Root1_(0//1)
  elseif c==-1 Root1_(1//2)
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
  c=conductor(r)
  if repl || TeX
    if c==1 print(io,"1")
    elseif c==2 print(io,"-1")
    else r=(TeX ? "\\zeta" : "ζ")*stringind(io,c)
      if d>=1 r*=stringexp(io,d) end
      print(io,r)
    end
  else
    print(io,d==1 ? "E($c)" : "E($c,$d)")
  end
end

function Base.cmp(a::Root1,b::Root1)
  r=cmp(conductor(a),conductor(b))
  if !iszero(r) return r end
  cmp(exponent(a),exponent(b))
end

Base.isless(a::Root1,b::Root1)=cmp(a,b)==-1
Base.isless(c::Root1,d::Real)=Cyc(c)<Cyc(d)
Base.isless(d::Real,c::Root1)=Cyc(d)<Cyc(c)
Base.:(==)(a::Root1,b::Root1)=a.r==b.r
Base.one(a::Root1)=Root1_(0//1)
Base.zero(a::Root1)=zero(Cyc{Int})
Base.isone(a::Root1)=iszero(a.r)
Base.:*(a::Root1,b::Root1)=Root1(;r=a.r+b.r)

Base.:^(a::Root1,n::Integer)=Root1(;r=n*a.r)
Base.:^(a::Root1,r::Rational{<:Integer})=root(a,denominator(r))^numerator(r)
Base.inv(a::Root1)=Root1(;r=-a.r)
Base.conj(a::Root1)=inv(a)
Base.:/(a::Root1,b::Root1)=a*inv(b)

#------------------------ type Cyc ----------------------------------
const impl=:MM # I tried 4 different implementations. For testmat(12)^2
    # ModuleElt is fastest
    # HModuleElt is 50% slower than ModuleElt
    # :svec is 20% slower than ModuleElt
    # :vec is 40% slower than ModuleElt

if impl==:vec
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  d::Vector{T} # the i-th element is the coefficient on ζⁱ⁻¹
  global function Cyc_(d::Vector{T}) where T<:Real 
    new{T}(d)
  end
end
function Cyc(c::Integer,v::AbstractVector)
  if c!=length(v) error("c=$c but length(v)=$(length(v))\n") end
  Cyc_(v)
end
conductor(c::Cyc)=length(c.d)
Base.pairs(c::Cyc)=Iterators.map(x->x[1]-1=>x[2],
                              Iterators.filter(x->x[2]!=0,enumerate(c.d)))
Base.getindex(c::Cyc,i::Integer)=c.d[i+1]
elseif impl==:svec
using SparseArrays
struct Cyc{T}<:Number    # a cyclotomic number
  d::SparseVector{T,Int} # d[i]==coeff on ζⁱ⁻¹ (where i∈ zumbroich_basis(n))
  global function Cyc_(d::SparseVector{T}) where T<:Real 
#   if !issorted(d.nzind) || any(iszero,d.nzval) error(d) end
    new{T}(d)
  end
end
conductor(c::Cyc)=length(c.d)
Base.pairs(c::Cyc)=Iterators.map(x->x[1]-1=>x[2],zip(c.d.nzind,c.d.nzval))
function Cyc(c::Integer,v::SparseVector)
  if c!=length(v) error("c=$c but length(v)=$(length(v))\n") end
  Cyc_(v)
end
Base.getindex(c::Cyc,i::Integer)=c.d[i+1]

elseif impl==:MM
using ModuleElts
const MM=ModuleElt # you can try with HModuleElt
struct Cyc{T <: Real}<: Number   # a cyclotomic number
  n::Int
  d::MM{Int,T} # list of pairs: i=>coeff on ζⁱ (where i∈ zumbroich_basis(n))
end
conductor(c::Cyc)=c.n
Base.pairs(c::Cyc)=c.d
Base.getindex(c::Cyc,i::Integer)=c.d[i]
end

"""
   `conductor(c::Cyc)`
   `conductor(a::AbstractArray)`

returns the smallest positive integer  n such that `c∈ ℚ (ζₙ)` (resp. all
elements of `a` are in `ℚ (ζₙ)`).

```julia-repl
julia> conductor(E(9))
9

julia> conductor([E(3),1//2,E(4)])
12
```
"""
conductor(a::AbstractArray)=lcm(conductor.(a))
conductor(i::Integer)=1 # for convenience
conductor(i::Rational)=1 # for convenience

const zumbroich_basis_dict=Dict(1=>[0]) # The Zumbroich basis is memoized
"""
  Cyclotomics.zumbroich_basis(n::Int) 

  returns  the Zumbroich basis of  ℚ (ζₙ) as the  vector of i in 0:n-1 such
  that `ζₙⁱ` is in the basis
"""
function zumbroich_basis(n::Int)::Vector{Int}
  get!(zumbroich_basis_dict,n) do
  if n==1 return [0] end
  function J(k::Int, p::Int) # see [Breuer] Rem. 1 p. 283
    if k==0 if p==2 return 0:0 else return 1:p-1 end
    elseif p==2 return 0:1
    else return div(1-p,2):div(p-1,2)
    end
  end
  nfact=factor(n)
  res=
  let J=J
    [[div(n*i,p^k) for i in J(k-1,p)] for (p,np) in nfact for k in 1:np]
  end
  v=sum.(vec(collect(Iterators.product(res...))))
  sort(v.%n)
  end
end

"""
`coefficients(c::Cyc)`

for  a cyclotomic `c` of conductor `n`,  returns a vector `v` of length `n`
such that `c==∑ᵢ vᵢ₋₁ ζⁱ`.

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
function coefficients(c::Cyc{T})where T
if impl==:svec return Array(c.d)
elseif impl==:vec return c.d
else
  res=zeros(T,conductor(c))
  for (i,v) in pairs(c) res[i+1]=v end
  res
end
end
  
"""
`denominator(c::Cyc{Rational})`

returns the smallest `d` such that `d*c` has integral coefficients (thus is
an algebraic integer).
"""
Base.denominator(c::Cyc)=lcm(denominator.(values(c.d)))

Base.numerator(c::Cyc{<:Rational{T}}) where T =Cyc{T}(c*denominator(c))

const Elist_dict=Dict{Tuple{Int,Int},Pair{Bool,Vector{Int}}}((1,0)=>(true=>
                       [impl==:MM ? 0 : 1])) # to memoize Elist
"""
  Cyclotomics.Elist(n,i)  
  
  expresses  ζₙⁱ  in  zumbroich_basis(n):  it  is  a  sum  of some ζₙʲ with
  coefficients all 1 or all -1. The result is a Pair sgn=>inds where sgn is
  true  if coefficients are all 1 and false otherwise, and inds is the list
  of i in 0:n-1 such that ζₙⁱ occurs with a non-zero coefficient (the i in
  1:n such that ζₙⁱ⁻¹.. for :vec and :svec)
"""
function Elist(n::Int,i1::Int=1)
  i=mod(i1,n)
  get!(Elist_dict,(n,i)) do
    mp=Int[]
    j=i
    for (p,np) in factor(n)
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
if impl==:svec || impl==:vec
    if isempty(mp) return true=> [i+1] end
elseif impl==:MM
    if isempty(mp) return true=> [i] end
end
    v=vec(sum.(Iterators.product((div(n,p)*(1:p-1) for p in mp)...)))
if impl==:svec || impl==:vec
    iseven(length(mp))=>sort((i.+v).%n).+1
elseif impl==:MM
    iseven(length(mp))=>sort((i.+v).%n)
end
  end
end

if impl==:vec
Cyc(i::Real)=Cyc(1,[i])
elseif impl==:svec
Cyc(i::Real)=Cyc_(i==0 ? spzeros(typeof(i),1) : SparseVector(1,[1],[i]))
else
Cyc(i::Real)=Cyc(1,MM(i==0 ? Pair{Int,typeof(i)}[] : [0=>i];check=false))
end
Cyc{T}(i::Real) where T<:Real=Cyc(T(i))

const E_dict=Dict(0//1=>Cyc(1))

function Cyc(a::Root1)
  get!(E_dict,a.r) do
    c=conductor(a)
    e=exponent(a)
    if c%4==2 return -E(div(c,2),div(e,2)+div(c+2,4)) end
    s,l=Elist(c,e) #::Pair{Bool,Vector{Int}}
if impl==:vec || impl==:svec
    v=zerocyc(Int,c)
    v[l].=ifelse(s,1,-1)
else
    v=MM(l.=>ifelse(s,1,-1);check=false)
end
    Cyc(c,v)
  end
end

Cyc{T}(a::Root1) where T<:Real=Cyc{T}(Cyc(a))

Base.zero(c::Cyc{T}) where T=Cyc{T}(0)
if impl==:svec
Base.iszero(c::Cyc)=nnz(c.d)==0 # faster than the other definition
else
Base.iszero(c::Cyc)=iszero(c.d)
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
if impl==:vec || impl==:svec
  Cyc(conductor(c),T.(c.d))
else
  Cyc(conductor(c),convert(MM{Int,T},c.d))
end
end

# num(c): value of c as a real when conductor(c)==1
if impl==:vec || impl==:svec
  num(c::Cyc)=c.d[1]
else
  num(c::Cyc{T}) where T =iszero(c) ? zero(T) : last(first(c.d))
end

function (::Type{T})(c::Cyc)where T<:Union{Integer,Rational}
  if conductor(c)==1 return T(num(c)) end
  throw(InexactError(:convert,T,c))
end

function Base.convert(::Type{T},c::Cyc;check=true)where T<:AbstractFloat
  if check && !isreal(c) throw(InexactError(:convert,T,c)) end
  real(convert(Complex{T},c))
end

function (::Type{T})(c::Cyc)where T<:AbstractFloat
  real(convert(Complex{T},c;check=true))
end

(::Type{T})(a::Root1) where T<:Number = T(Cyc(a))

function Complex{T}(c::Cyc)where T<:AbstractFloat
  sum(((k,v),)->v*cispi(2*T(k)/conductor(c)),pairs(c);init=Complex{T}(0.0))
end

function Complex{T}(c::Cyc)where T<:Union{Integer,Rational}
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
Base.complex(c::Cyc{T}) where T =(conductor(c)==1||conductor(c)==4) ? Complex{T}(c) : Complex{float(T)}(c)
Base.complex(a::Root1)=complex(Cyc(a))

Base.isinteger(c::Cyc)=conductor(c)==1 && isinteger(num(c))
Base.isinteger(a::Root1)=a.r==0//1 || a.r==1//2

Base.isreal(c::Cyc)=conductor(c)==1 || c==conj(c)

function Base.real(c::Cyc{T}) where T<:Real
 if conductor(c)==1 return num(c) end
  (c+conj(c))/2
end

Base.real(a::Root1)=real(Cyc(a))

function Base.imag(c::Cyc{T}) where T<:Real
 if conductor(c)==1 return 0 end
  (c-conj(c))/2
end

Base.imag(a::Root1)=imag(Cyc(a))

# addroot: add c*E(n,i) to res
function addroot(res,n,i,c)
  (s,v)=Elist(n,i)
  if !s c=-c end
if impl==:MM
  for k in v push!(res,k=>c) end
elseif impl==:vec || impl==:svec
@inbounds  view(res,v).+=c
end
end

# zerocyc: initialise a proper res for addroot for conductor==n
if impl==:vec
zerocyc(::Type{T},n) where T=zeros(T,n)
elseif impl==:svec
zerocyc(::Type{T},n) where T=spzeros(T,n)
elseif impl==:MM
zerocyc(::Type{T},n) where T=T[]
end

function raise(n::Int,c::Cyc) # write c in Q(ζ_n) if conductor(c) divides n
  if n==conductor(c) return c end
  m=div(n,conductor(c))
  res=zerocyc(eltype(c.d),n)
  for (i,v) in pairs(c)
    addroot(res,n,i*m,v)
  end
  Cyc(n,impl==:MM ? MM(res) : res)
end

function promote_conductor(a::Cyc,b::Cyc)
  if conductor(a)==conductor(b) return (a, b) end
  l=lcm(conductor(a),conductor(b))
  (raise(l,a),raise(l,b))
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:AbstractFloat}
  Complex{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:AbstractFloat}
  Complex{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Complex{T2}})where {T1,T2<:Union{Integer,Rational{<:Integer}}}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{T2})where {T1,T2<:Union{Integer,Rational{<:Integer}}}
  Cyc{promote_type(T1,T2)}
end

function Base.promote_rule(a::Type{Cyc{T1}},b::Type{Cyc{T2}})where {T1,T2}
  Cyc{promote_type(T1,T2)}
end

# total order is necessary to put Cycs in a sorted list
# for conductor(c)==1  a<b is as expected
function Base.cmp(a::Cyc,b::Cyc)
  t=cmp(conductor(a),conductor(b))
  if !iszero(t) return t end
  conductor(a)==1 ? cmp(num(a),num(b)) : cmp(a.d,b.d)
end

Base.:(==)(a::Cyc,b::Cyc)=conductor(a)==conductor(b) && a.d==b.d
Base.isless(a::Cyc,b::Cyc)=cmp(a,b)==-1
Base.isless(c::Cyc,d::Real)=c<Cyc(d)
Base.isless(d::Real,c::Cyc)=Cyc(d)<c

# hash is necessary to put Cycs as keys of a Dict or make a Set
Base.hash(a::Cyc, h::UInt)=hash(a.d, hash(conductor(a), h))

function Base.show(io::IO, ::MIME"text/html", a::Cyc)
  print(io, "\$")
  show(IOContext(io,:TeX=>true),a)
  print(io, "\$")
end

function Base.show(io::IO, ::MIME"text/plain", r::Cyc)
  if !haskey(io,:typeinfo) print(io,typeof(r),": ") end
  show(io,r)
end

function normal_show(io::IO,p::Cyc{T})where T
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if T<:Rational{<:Integer}
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
        r=(TeX ? "\\zeta" : "ζ") * stringind(io,conductor(p))
        r*= stringexp(io,deg)
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

function Base.show(io::IO, p::Cyc{T})where T
  quadratic=get(io,:quadratic,true)
  repl=get(io,:limit,false)
  TeX=get(io,:TeX,false)
  if conductor(p)==1
    n=num(p)
    if T<:Integer || (T<:Rational{<:Integer} && denominator(n)==1)
      if repl||TeX||haskey(io,:typeinfo) print(io,numerator(n))
      else print(io,"Cyc{",T,"}(",numerator(n),")")
      end
      return
    end
  end
  rqq=[normal_show(io,p)]
  if quadratic && (T<:Integer || T<:Rational{<:Integer})
    q=Quadratic(p)
    if !isnothing(q) push!(rqq,repr(q;context=io)) end
    for test in [1-E(4),1+E(4),E(3),E(3,2),1-E(3),1-E(3,2),1+E(3),1+E(3,2),root(-3)]
      if !iszero(conductor(p)%conductor(test)) continue end
      q=Quadratic(p/test)
      if isnothing(q) continue end
      rq=repr(q;context=io)
      rq=format_coefficient(rq;allow_frac=true)
      t=format_coefficient(normal_show(io,test))
      if !isempty(rq) && rq[1]=='-' rq="-"*t*rq[2:end] else rq=t*rq end
      push!(rqq,rq)
    end
  end
  print(io,rqq[argmin(length.(rqq))])
end

Base.gcd(v::Vector{<:Cyc})=one(v[1])
Base.gcd(a::Cyc,b::Cyc)=one(a)
Base.gcd(a::Cyc,b::Number)=one(a)
Base.gcd(b::Number,a::Cyc)=one(a)

function Base.:+(x::Cyc,y::Cyc)
  a,b=promote(x,y)
  if iszero(a) return b
  elseif iszero(b) return a
  end
if impl==:vec
  a,b=promote_conductor(a,b)
  lower(Cyc(conductor(a),a.d+b.d))
elseif impl==:svec
  a,b=promote_conductor(a,b)
  lower(Cyc_(dropzeros!(a.d+b.d)))
# n=lcm(conductor(a),conductor(b))
# na=div(n,conductor(a))
# nb=div(n,conductor(b))
# res=zerocyc(eltype(a.d),n)
# for (i,va) in pairs(a) addroot(res,n,na*i,va) end
# for (i,vb) in pairs(b) addroot(res,n,nb*i,vb) end
# lower(Cyc(n,res))
else
  n=lcm(conductor(a),conductor(b))
  res=eltype(a.d)[]
  na=div(n,conductor(a))
  nb=div(n,conductor(b))
  for (i,va) in pairs(a) addroot(res,n,na*i,va) end
  for (i,vb) in pairs(b) addroot(res,n,nb*i,vb) end
  lower(Cyc(n,MM(res)))
end
end

Base.:-(a::Cyc)=Cyc(conductor(a),-a.d)
Base.:-(a::Cyc,b::Cyc)=a+(-b)

if impl==:vec || impl==:svec
Base.://(c::Cyc,a::Real)=Cyc(conductor(c),c.d.//a)
else
Base.://(c::Cyc,a::Real)=Cyc(conductor(c),c.d//a)
end

function Base.div(c::Cyc,a::Real)
if impl==:MM
  n=merge(div,c.d,a)
  Cyc(iszero(n) ? 1 : conductor(c),n)
else
  res=div.(c.d,a)
  iszero(res) ? zerocyc(eltype(res),1) : Cyc(conductor(c),res)
end
end

function Base.:*(c::Cyc,a::Real)
if impl==:MM
  Cyc(iszero(a) ? 1 : conductor(c),c.d*a)
else
  res=c.d*a
  iszero(a) ? zero(Cyc{eltype(res)}) : Cyc(conductor(c),res)
end
end
Base.:*(a::Real,c::Cyc)=c*a

Base.://(a::Cyc,c::Cyc)=a*inv(c)
Base.://(a::Real,c::Cyc)=a*inv(c)
Base.://(c::Root1,a::Real)=Cyc(c)//a
Base.://(a::Real,c::Root1)=a//Cyc(c)
Base.://(c::Root1,a::Cyc)=Cyc(c)//a
Base.://(a::Cyc,c::Root1)=a//Cyc(c)
Base.:/(c::Cyc,a::Real)=c//a
Base.:/(a::Cyc,c::Cyc)=a//c
Base.:/(a::Real,c::Cyc)=a//c

function Base.:*(a::Cyc,b::Cyc)
  a,b=promote(a,b)
  if iszero(a) return a end
  if iszero(b) return b end
  if conductor(a)==1 return b*num(a)
  elseif conductor(b)==1 return a*num(b)
  end
  n=lcm(conductor(a),conductor(b))
  res=zerocyc(eltype(a.d),n)
  na=div(n,conductor(a))
  nb=div(n,conductor(b))
  for (i,ai) in pairs(a), (j,bj) in pairs(b) 
    addroot(res,n,na*i+nb*j,ai*bj)
  end
  lower(Cyc(n,impl==:MM ? MM(res) : impl==:svec ? dropzeros!(res) : res))
end

function lower(c::Cyc{T})where T # write c in smallest Q(ζ_n) where it sits
  n=conductor(c)
 # println("lowering $(conductor(c)):$(c.d)")
  if n==1 return c end
  if iszero(c) return zero(Cyc{T}) end
  for (p,np) in factor(n)
    m=div(n,p)
if impl==:vec
    kk=filter(i->c.d[i]!=0,eachindex(c.d))
    val=c.d[kk]
    if np>1 
      if all(k->(k-1)%p==0,kk) 
        v=zeros(eltype(c.d),m)
        view(v,map(x->1+div(x-1,p),kk)).=val
        return lower(Cyc(m,v))
      end
    elseif iszero(length(kk)%(p-1))
      kk=kk.-1
      cnt=zeros(Int,m)
      for k in kk cnt[1+(k%m)]+=1 end
      if !all(x->iszero(x) || x==p-1,cnt) continue end
      u=findall(!iszero,cnt).-1
      kk=@. div(u+m*mod(-u,p)*invmod(m,p),p)%m
      if !issorted(kk) sort!(kk) end
      v=zeros(eltype(c.d),m)
      if p==2  
        view(v,kk.+1).=[c.d[1+(k*p)%n] for k in kk]
        return lower(Cyc(m,v))
      elseif all(k->constant(map(i->c.d[1+(m*i+k*p)%n],1:p-1)),kk)
        view(v,kk.+1).=[-c.d[1+(m+k*p)%n] for k in kk]
        return lower(Cyc(m,v))
      end
    end
elseif impl==:svec
    kk=c.d.nzind
    val=c.d.nzval
    if np>1 
      if all(k->(k-1)%p==0,kk) 
        return lower(Cyc(m,SparseVector(m,map(x->1+div(x-1,p),kk),val)))
      end
    elseif iszero(length(kk)%(p-1))
      kk=kk.-1
      cnt=zeros(Int,m)
      for k in kk cnt[1+(k%m)]+=1 end
      if !all(x->iszero(x) || x==p-1,cnt) continue end
      u=findall(!iszero,cnt).-1
      kk=@. div(u+m*mod(-u,p)*invmod(m,p),p)%m
      if !issorted(kk) sort!(kk) end
      if p==2  
        return lower(Cyc(m,SparseVector(m,1 .+kk,[c.d[1+(k*p)%n] for k in kk])))
      elseif all(k->constant(map(i->c.d[1+(m*i+k*p)%n],1:p-1)),kk)
        return lower(Cyc(m,SparseVector(m,1 .+kk,[-c.d[1+(m+k*p)%n] for k in kk])))
      end
    end
elseif impl==:MM
    if np>1 
      if all(k->first(k)%p==0,c.d) 
        return lower(Cyc(m,MM(div(k,p)=>v for (k,v) in c.d;check=false)))
      end
    elseif iszero(length(c.d)%(p-1))
      cnt=zeros(Int,m)
      for (k,v) in c.d cnt[1+(k%m)]+=1 end
      if !all(x->iszero(x) || x==p-1,cnt) continue end
      u=findall(!iszero,cnt).-1
      kk=@. div(u+m*mod(-u,p)*invmod(m,p),p)%m
      if p==2  
        return lower(Cyc(m,MM(k=>c.d[(k*p)%n] for k in kk)))
      elseif all(k->constant(map(i->c.d[(m*i+k*p)%n],1:p-1)),kk)
        return lower(Cyc(m,MM(k=>-c.d[(m+k*p)%n] for k in kk)))
      end
    end
end
  end
  c
end

galois(c::Rational,n::Int)=c

galois(c::Integer,n::Int)=c

"""
  galois(c::Cyc,n::Int) applies to c the galois automorphism
  of Q(ζ_conductor(c)) raising all roots of unity to the n-th power.
  n should be prime to conductor(c).
# Examples
```julia-repl
julia> galois(1+E(4),-1) # galois(c,-1) is the same as conj(c)
Cyc{Int64}: 1-ζ₄

julia> galois(root(5),2)==-root(5)
true
```
"""
function galois(c::Cyc,n::Int)
  if gcd(n,conductor(c))!=1 
    error("$n should be prime to conductor($c)=$(conductor(c))")
  end
  if iszero(c) return c end
  res=zerocyc(eltype(c.d),conductor(c))
  for (e,p) in pairs(c)
    addroot(res,conductor(c),(e*n)%conductor(c),p)
  end
  Cyc(conductor(c),impl==:MM ? MM(res) : res)
end

function galois(c::Root1,n::Int)
  d=denominator(c.r)
  if gcd(n,d)!=1 error("$n should be prime to conductor($c)=$d") end
  Root1(;r=n*c.r)
end

Base.conj(c::Cyc)=galois(c,-1)

function Base.inv(c::Cyc)
  if conductor(c)==1
    r=num(c)
    if r==1 || r==-1 return Cyc(r) else return Cyc(1//r) end
  end
  l=setdiff(unique(map(i->galois(c,i),prime_residues(conductor(c)))),[c])
  r=prod(l)
  n=num(c*r)
  n==1 ? r : (n==-1 ? -r : r//n)
end

Base.:^(a::Cyc, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                   Base.power_by_squaring(inv(a),-n)

Base.abs(c::Cyc)=c*conj(c)

"""
`Root1(c)`
    
`c` should be a cyclotomic number (a `Cyc`), or a `Real`. `Root1` returns a
`Root1` object containing the rational `e/n` with `0≤e<n` (that is, `e/n∈ ℚ
/ℤ`) if `c==E(n,e)`, and `nothing` if `c` is not a root of unity.

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
""" 
function Root1(c::Cyc)
  if !(all(x->last(x)==1,pairs(c)) || all(x->last(x)==-1,pairs(c)))
    return nothing
  end
  for i in prime_residues(conductor(c))
    if c==E(conductor(c),i) return Root1_(i//conductor(c)) end
    if -c==E(conductor(c),i)
      if conductor(c)%2==0 return E(conductor(c),div(conductor(c),2)+i)
      else return E(2conductor(c),conductor(c)+2*i)
      end
    end
  end
  # return nothing
end

Base.:(==)(a::Root1,b::Number)=Cyc(a)==b
function Base.:*(a::Cyc,b::Root1)
  n=lcm(conductor(a),conductor(b))
  na=div(n,conductor(a))
  nb=div(n,conductor(b))
  res=zerocyc(eltype(a.d),n)
  for (i,va) in pairs(a) addroot(res,n,na*i+nb*exponent(b),va) end
  lower(Cyc(n,impl==:MM ? MM(res) : impl==:svec ? dropzeros!(res) : res))
end
Base.:*(b::Root1,a::Cyc)=a*b

Base.:+(a::Root1,b::Root1)=Cyc(a)+Cyc(b)
Base.:-(a::Root1,b::Root1)=Cyc(a)-Cyc(b)
Base.:-(r::Root1)=-Cyc(r)
Base.promote_rule(a::Type{Root1},b::Type{Cyc{T}}) where T =b
Base.promote_rule(a::Type{Root1},b::Type{<:Integer})=Cyc{b}
Base.promote_rule(a::Type{Root1},b::Type{<:Rational{<:Integer}})=Cyc{b}
Base.promote_rule(a::Type{Root1},b::Type{<:AbstractFloat})=Complex{b}
Base.promote_rule(a::Type{Root1},b::Type{Complex{T}}) where T =Cyc{promote_type(T,Int)}
#------------------- end of Root1 ----------------------------------------

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
`q.root`,  `q.den` representing `c` as `(q.a + q.b root(q.root))//q.den` if
such a representation is possible or returns `q===nothing` otherwise.

# Examples
```julia-repl
julia> Quadratic(1+E(3))
(1+√-3)/2

julia> Quadratic(1+E(5))

```
"""
function Quadratic(c::Cyc{T})where T
  l1=coefficients(c)
  den=lcm(denominator.(l1))
  c=Cyc{typeof(den)}(c*den)
  if conductor(c)==1 return Quadratic(numerator(num(c)),0,1,den) end
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
    gal=Set(galois.(c,prime_residues(conductor(c))))
    if length(gal)!=2 return nothing end
    a=numerator(convert(T,sum(gal)))      # trace of 'c' over the rationals
    if length(f)%2==0 b=2*c[1]-a
    else b=2*c[1]+a
    end
    if a&1==0 && b&1==0 a>>=1; b>>=1; d=1
    else d=2
    end
  elseif v2==2
    sqr=conductor(c)>>2
    if sqr==1 a=c[0];b=-c[1]
    else
      a=c[4]
      if length(f)%2==0 a=-a end
      b=-c[sqr+4]
    end
    if sqr%4==1 sqr=-sqr; b=-b end
    d=1
  else		# v2 = 3
    sqr=conductor(c)>>2
    if sqr==2
      a=c[0];b=c[1]
      if b==c[3] sqr=-2 end
    else
      a=c[8]
      if length(f)%2==0 a=-a end
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
xprint(x...;p...)=print(IOContext(io,:limit=>true,p...),x...)
function proot(x,n,r)
  xprint("root(",x)
  if n!=2 xprint(",",n) end
  xprint(")=",r,"\n")
end
const Irootdict=Dict{Tuple{Int,Int},Union{Cyc{Int},Int}}()

"""
`root(x,n=2)`

computes  the `n`-th root of `x` when we know  how to do it. We know how to
compute  `n`-th  roots  for  roots  of  unity, square roots of integers and
`n`-th  roots  of  perfect  `n`-th  powers  of  integers or square roots of
integers.

```julia-repl
julia> root(-1)
Cyc{Int64}: ζ₄

julia> root(E(4))
Root1: ζ₈

julia> root(27,6)
Cyc{Int64}: √3
```
"""
function root(x::Integer,n=2)
  if isone(n) || isone(x) return x end
  if !(n isa Int) n=Int(n) end
  get!(Irootdict,(n,x)) do
    if x==1 || (x==-1 && n%2==1) return x end
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
  d=conductor(r)
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
  if isone(n) || isone(x) return x end
  if !(n isa Int) n=Int(n) end
  get!(Crootdict,(n,x)) do
  d=denominator(x)
  if d!=1 return root(numerator(x),n)/root(d,n) end
  d=gcd(coefficients(x))
  if d!=1 return root(div(x,d),n)*root(d,n) end
    r=Root1(x)
    if isnothing(r) 
      if conductor(x)>1 return nothing end
      return root(num(x),n)
    end
    res=Cyc(root(r,n))
    if inforoot[] proot(x,n,res) end
    res
  end
end

# 347.534 ms (4367402 allocations: 366.17 MiB) in 1.5.3
# 565.431 ms (5861810 allocations: 775.28 MiB) in 1.5.3, HModuleElts
function testmat(p)
  ss=[[i,j] for i in 0:p-1 for j in i+1:p-1]
  [(E(p,i'*reverse(j))-E(p,i'*j))//p for i in ss,j in ss]
end
end

# Property tests for CyclotomicNumbers.jl.
#
# These check the algebra itself, not the printing: every claim about a Cyc is
# verified against an independent Complex{BigFloat} evaluation of its stored
# coefficients, so a wrong normal form, a missed lowering or a bad Galois action
# shows up as a failure. They use only the public API plus Elist/zumbroich_basis/
# lower! (the three internals the module documents), and no dependency beyond
# Test, so they can be included from runtests.jl as they are.

# a deterministic generator, so the tests need no `using Random`
mutable struct SmallRNG
  s::UInt64
end
SmallRNG(seed::Integer=20240101)=SmallRNG(UInt64(seed))
function randint!(g::SmallRNG,lo::Int,hi::Int)
  g.s=g.s*6364136223846793005+1442695040888963407
  lo+Int((g.s>>33)%UInt64(hi-lo+1))
end

# value of a Cyc as a complex number, optionally with ζ_n ↦ ζ_n^k applied
function cycvalue(c::Cyc,k::Integer=1)
  n=conductor(c); s=Complex{BigFloat}(0)
  for (i,v) in pairs(c)
    s+=Complex{BigFloat}(BigFloat(v))*cis(2*BigFloat(pi)*BigFloat(mod(i*k,n))/n)
  end
  s
end
cycvalue(x::Number,k::Integer=1)=Complex{BigFloat}(BigFloat(real(x)),BigFloat(imag(x)))
samevalue(a,b;tol=BigFloat(10)^-40)=abs(cycvalue(a)-cycvalue(b))<tol

# ζ_n^i written in Q(ζ_n) *without* lowering (that is what Elist is for).
# raise exists in the versions that split the raw constructor out; otherwise
# build the pair list directly.
if isdefined(CyclotomicNumbers,:raise)
  rawroot(n::Int,i::Int)=CyclotomicNumbers.raise(n,Cyc(E(n,i)))
  embedd(n::Int,y::Cyc)=CyclotomicNumbers.raise(n,y)
else
  rawroot(n::Int,i::Int)=CyclotomicNumbers.rawCyc(n,Int,(mod(i,n)=>1,))
  embedd(n::Int,y::Cyc)=CyclotomicNumbers.rawCyc(n,Int,
    (i*div(n,conductor(y))=>v for (i,v) in pairs(y)))
end

eulerphi(n::Integer)=n==1 ? 1 : prod(p^(e-1)*(p-1) for (p,e) in CyclotomicNumbers.factor(n))
randomcyc(g::SmallRNG,n::Int)=sum(randint!(g,-3,3)*Cyc(E(n,i)) for i in 0:n-1)

@testset "CyclotomicNumbers properties" verbose=true begin
  rng=SmallRNG()

  @testset "Elist expands ζ_n^i correctly" begin
    for n in 1:60, i in 0:n-1
      @test abs(cycvalue(rawroot(n,i))-cis(2*BigFloat(pi)*BigFloat(i)/n))<BigFloat(10)^-40
    end
  end

  @testset "zumbroich_basis is a basis and matches Elist" begin
    for n in 1:120
      @test length(CyclotomicNumbers.zumbroich_basis(n))==eulerphi(n)
    end
    for n in 1:60
      used=sort(unique(vcat([CyclotomicNumbers.Elist(n,i)[2] for i in 0:n-1]...)))
      @test used==CyclotomicNumbers.zumbroich_basis(n)
    end
  end

  @testset "conductor is minimal" begin
    # x ∈ Q(ζ_{n/p}) iff every k ≡ 1 mod n/p, k ≠ 1 mod n, fixes x
    for n in 1:60, _ in 1:4
      x=randomcyc(rng,n); cn=conductor(x); cn==1 && continue
      for p in unique(keys(CyclotomicNumbers.factor(cn)))
        d=div(cn,p)
        ks=[k for k in 1:cn-1 if gcd(k,cn)==1 && mod(k,d)==1 && mod(k,cn)!=1]
        isempty(ks) && continue
        @test !all(k->abs(cycvalue(x,k)-cycvalue(x))<BigFloat(10)^-40,ks)
      end
    end
  end

  @testset "lower! recovers the small field" begin
    for n in 2:40, c in 1:n-1
      n%c==0 || continue
      for _ in 1:2
        y=randomcyc(rng,c); conductor(y)==c || continue
        x=embedd(n,y)
        @test samevalue(x,y)                    # embedding preserves the value
        xl=CyclotomicNumbers.lower!(deepcopy(x))
        @test conductor(xl)==c                  # ... and lowering recovers c
        @test samevalue(xl,y)
      end
    end
  end

  @testset "ring axioms and values" begin
    ns=[1,2,3,4,5,6,8,9,10,12,15,16,18,20,21,24,25,27,30,36,45,60]
    for _ in 1:80
      x,y,z=randomcyc(rng,rand(ns)),randomcyc(rng,rand(ns)),randomcyc(rng,rand(ns))
      @test (x+y)+z==x+(y+z)
      @test (x*y)*z==x*(y*z)
      @test x*(y+z)==x*y+x*z
      @test x*y==y*x
      @test conj(conj(x))==x
      @test samevalue(x*y,cycvalue(x)*cycvalue(y))
      @test x-x==zero(Cyc{Int})
      @test conductor(x-x)==1
      @test iszero(x-x)
    end
  end

  @testset "galois action" begin
    for n in 1:30, k in 1:n-1
      gcd(k,n)==1 || continue
      x=randomcyc(rng,n); cn=conductor(x)
      s=sum((v*Cyc(E(cn,mod(i*k,cn))) for (i,v) in pairs(x));init=zero(Cyc{Int}))
      @test galois(x,k)==s
      @test conductor(galois(x,k))==cn
    end
  end

  @testset "roots of unity" begin
    for n in 1:18, i in 0:n-1, k in 1:6
      @test root(E(n,i),k)^k==E(n,i)
    end
    for n in 1:40, i in 0:n-1
      @test Root1(Cyc(E(n,i)))==E(n,i)
    end
  end

  @testset "Quadratic is complete on quadratic subfields" begin
    for (D,_) in [(-3,3),(5,5),(-7,7),(-4,4),(8,8),(12,12),(13,13),(-11,11),(17,17),
                  (-8,8),(24,24),(-15,15),(21,21),(-20,20),(28,28),(33,33),(-24,24),
                  (-3,9),(5,15),(13,39),(-7,21),(-3,27),(5,45),(17,51)]
      sq=root(D)
      for _ in 1:5
        c=randint!(rng,-4,4)+randint!(rng,-4,4)*sq
        iszero(c) && continue
        q=Quadratic(c)
        @test !isnothing(q)
        isnothing(q) || @test q.den*c==q.a+q.b*root(q.root)
      end
    end
  end
end

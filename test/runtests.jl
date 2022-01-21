# auto-generated tests from julia-repl docstrings
using Test, CyclotomicNumbers
function mytest(f::String,a::String,b::String)
  println(f," ",a)
  omit=a[end]==';'
  a=replace(a,"\\\\"=>"\\")
  a=repr(MIME("text/plain"),eval(Meta.parse(a)),context=:limit=>true)
  if omit a="nothing" end
  a=replace(a,r" *(\n|$)"s=>s"\1")
  a=replace(a,r"\n$"s=>"")
  b=replace(b,r" *(\n|$)"s=>s"\1")
  b=replace(b,r"\n$"s=>"")
  i=1
  while i<=lastindex(a) && i<=lastindex(b) && a[i]==b[i]
    i=nextind(a,i)
  end
  if a!=b print("exec=$(repr(a[i:end]))\nmanl=$(repr(b[i:end]))\n") end
  a==b
end
@testset verbose = true "Gapjm" begin
@testset "CyclotomicNumbers.jl" begin
@test mytest("CyclotomicNumbers.jl","E(3,2)","Root1: ζ₃²")
@test mytest("CyclotomicNumbers.jl","E(3)+E(4)","Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹")
@test mytest("CyclotomicNumbers.jl","E(12,11)-E(12,7)","Cyc{Int64}: √3")
@test mytest("CyclotomicNumbers.jl","repr(E(12,11)-E(12,7),context=(:limit=>true,:quadratic=>false))","\"-ζ₁₂⁷+ζ₁₂¹¹\"")
@test mytest("CyclotomicNumbers.jl","a=E(3)+E(3,2)","Cyc{Int64}: -1")
@test mytest("CyclotomicNumbers.jl","conductor(a)","1")
@test mytest("CyclotomicNumbers.jl","typeof(Int(a))","Int64")
@test mytest("CyclotomicNumbers.jl","inv(1+E(4))","Cyc{Float64}: 0.5-0.5ζ₄")
@test mytest("CyclotomicNumbers.jl","1//(1+E(4))","Cyc{Rational{Int64}}: (1-ζ₄)/2")
@test mytest("CyclotomicNumbers.jl","Cyc(1//2+im)","Cyc{Rational{Int64}}: (1+2ζ₄)/2")
@test mytest("CyclotomicNumbers.jl","conj(1+E(4))","Cyc{Int64}: 1-ζ₄")
@test mytest("CyclotomicNumbers.jl","real(E(3))","Cyc{Rational{Int64}}: -1/2")
@test mytest("CyclotomicNumbers.jl","Rational{Int}(real(E(3)))","-1//2")
@test mytest("CyclotomicNumbers.jl","imag(E(3))","Cyc{Rational{Int64}}: √3/2")
@test mytest("CyclotomicNumbers.jl","c=Cyc(E(9))","Cyc{Int64}: -ζ₉⁴-ζ₉⁷")
@test mytest("CyclotomicNumbers.jl","Root1(c)","Root1: ζ₉")
@test mytest("CyclotomicNumbers.jl","Root1(1+E(4))","nothing")
@test mytest("CyclotomicNumbers.jl","Root1(;r=1//4)","Root1: ζ₄")
@test mytest("CyclotomicNumbers.jl","c=E(4)*E(3)","Root1: ζ₁₂⁷")
@test mytest("CyclotomicNumbers.jl","c=Complex{Float64}(E(3))","-0.4999999999999999 + 0.8660254037844387im")
@test mytest("CyclotomicNumbers.jl","0.0+E(3)","Cyc{Float64}: 1.0ζ₃")
@test mytest("CyclotomicNumbers.jl","E(3)+1//2","Cyc{Rational{Int64}}: √-3/2")
@test mytest("CyclotomicNumbers.jl","E(3)+im","Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹")
@test mytest("CyclotomicNumbers.jl","complex(E(4))","0 + 1im")
@test mytest("CyclotomicNumbers.jl","complex(E(3))","-0.4999999999999999 + 0.8660254037844387im")
@test mytest("CyclotomicNumbers.jl","-1<Cyc(0)<1","true")
@test mytest("CyclotomicNumbers.jl","a=E(3)+E(4)","Cyc{Int64}: ζ₁₂⁴-ζ₁₂⁷-ζ₁₂¹¹")
@test mytest("CyclotomicNumbers.jl","collect(pairs(a))","3-element Vector{Pair{Int64, Int64}}:\n  4 => 1\n  7 => -1\n 11 => -1")
@test mytest("CyclotomicNumbers.jl","a[6],a[7]","(0, -1)")
@test mytest("CyclotomicNumbers.jl","coefficients(a)","12-element Vector{Int64}:\n  0\n  0\n  0\n  0\n  1\n  0\n  0\n -1\n  0\n  0\n  0\n -1")
@test mytest("CyclotomicNumbers.jl","valtype(a)","Int64")
@test mytest("CyclotomicNumbers.jl","conductor(E(9))","9")
@test mytest("CyclotomicNumbers.jl","conductor([E(3),1//2,E(4)])","12")
@test mytest("CyclotomicNumbers.jl","coefficients(Cyc(E(9)))","9-element Vector{Int64}:\n  0\n  0\n  0\n  0\n -1\n  0\n  0\n -1\n  0")
@test mytest("CyclotomicNumbers.jl","galois(1+E(4),-1)","Cyc{Int64}: 1-ζ₄")
@test mytest("CyclotomicNumbers.jl","galois(root(5),2)==-root(5)","true")
@test mytest("CyclotomicNumbers.jl","conjugates(1+root(5))","2-element Vector{Cyc{Int64}}:\n 1+√5\n 1-√5")
@test mytest("CyclotomicNumbers.jl","r=Root1(-E(9,2)-E(9,5))","Root1: ζ₉⁸")
@test mytest("CyclotomicNumbers.jl","order(r)","9")
@test mytest("CyclotomicNumbers.jl","exponent(r)","8")
@test mytest("CyclotomicNumbers.jl","Cyc(r)","Cyc{Int64}: -ζ₉²-ζ₉⁵")
@test mytest("CyclotomicNumbers.jl","Root1(-E(9,4)-E(9,5))","nothing")
@test mytest("CyclotomicNumbers.jl","Quadratic(E(3,2)-2E(3))","(1-3√-3)/2")
@test mytest("CyclotomicNumbers.jl","Quadratic(1+E(5))","nothing")
@test mytest("CyclotomicNumbers.jl","root(-1)","Cyc{Int64}: ζ₄")
@test mytest("CyclotomicNumbers.jl","root(E(4))","Root1: ζ₈")
@test mytest("CyclotomicNumbers.jl","root(27,6)","Cyc{Int64}: √3")
end
end

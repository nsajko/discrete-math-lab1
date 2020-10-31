# Copyright 2020 Neven Sajko. All rights reserved.

module Lab1

using Printf
using LinearAlgebra
using SpecialMatrices

export main

const I = BigInt
const F = BigFloat
const V = Vector{F}

# Mutates neither x nor a.
lab1_method1(x::V, a::V, n::I)::F = dot(Vandermonde(x)' \ a, map(z::F -> z^n, x))

# Mutates a, but not x.
function lab1_method2(x::V, a::V, n::I)::F
	n == 0 && return a[1]
	n == 1 && return a[2]
	n == 2 && return a[3]
	local c1 = x[1] + x[2] + x[3]
	local c2 = -(x[1]*x[2] + x[2]*x[3] + x[3]*x[1])
	local c3 = x[1]*x[2]*x[3]
	for i in 3:n
		a = (a[2], a[3], c1*a[3] + c2*a[2] + c3*a[1])
	end
	a[3]
end

# Reads a line from stdin and parses the floating point number
# that is assumed to be contained in the line.
rd_lf()::F = BigFloat(readline(), precision=8192)

function main()
	local x::V = F[0, 0, 0]
	local a::V = F[0, 0, 0]

	print("Unesite prvo rješenje x_0 karakteristične jednadžbe: ")
	x[1] = rd_lf()
	print("Unesite drugo rješenje x_1 karakteristične jednadžbe: ")
	x[2] = rd_lf()
	print("Unesite treće rješenje x_2 karakteristične jednadžbe: ")
	x[3] = rd_lf()

	print("Unesite vrijednost nultog člana niza a_0: ")
	a[1] = rd_lf()
	print("Unesite vrijednost prvog člana niza a_1: ")
	a[2] = rd_lf()
	print("Unesite vrijednost drugog člana niza a_2: ")
	a[3] = rd_lf()

	print("Unesite redni broj n traženog člana niza: ")
	local n::I = parse(BigInt, readline())

	local out1 = lab1_method1(x, a, n)
	local out2 = lab1_method2(x, a, n)
	@printf("Vrijednost n-tog člana niza pomoću formule: %.15g\n", out1)
	@printf("Vrijednost n-tog člana niza iz rekurzije: %.15g\n", out2)
end

end

Lab1.main()

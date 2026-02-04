using SpecialFunctions

function double_factorial(m::Int)
	if m <= 1
		return 1.0
	end
	if isodd(m)
		n = (m + 1) ÷ 2
		return 2.0^n * gamma(n + 0.5) / sqrt(pi)
	else
		k = m ÷ 2
		return 2.0^k * factorial(k)
	end
end

function boys(m::Int, x::Float64)
	if x == 0.0
		return 1.0 / (2m + 1)
	end
	if x < 1e-8
		return 1.0 / (2m + 1) - x / (2m + 3)
	end
	s = m + 0.5
	return 0.5 * x^(-s) * (gamma(s) - gamma(s, x))
end

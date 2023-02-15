using Plots

x = Float64[]
y = Float64[]
for i in 1:10.:100, j in 1:10.:100
    a, b = randn(2) .* 1
    # a, b = 0, 0 
    push!(x, i+a)
    push!(y, j+b)
end

scatter(x, y, legen=:none)
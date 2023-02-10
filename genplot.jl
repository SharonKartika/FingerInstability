using Plots
using DelimitedFiles
using Dates
gr()
theme(:juno)

function generateAnimation(W, H)
    data = readdlm("intermediateResults/positionData.csv", ',', Float64)
    date = Dates.format(Dates.now(), "dd-mm-yy-HHMMSS")
    begin
        anim = @animate for i = 1:size(data)[1]
            x, y = data[i, 1:2:end], data[i, 2:2:end]
            # x = x .- W.*round.(x./W);
            # y = y .- H.*round.(y./H);
            plot(x, y, seriestype=:scatter,
                xlims=[-W / 2, W / 2],
                ylims=[-H / 2, H / 2],
                legend=nothing,
                # axis=nothing,
                markerstrokewidth=0,
                c=:royalblue,
                aspectratio=:equal)
            # plot(x ,y ,seriestype=:scatter,
            # legend=nothing, axis=nothing,
            # markerstrokewidth=0,c=:royalblue,
            # aspectratio=:equal)
        end
        gif(anim, "results/Result$(date).gif", fps=30)
    end
end

function simulateAndPlot(; N=200, nt=100, W=1200, H=1200)
    rt = 70.0
    t1 = @elapsed run(`g++ main.cpp -o main`)
    t2 = @elapsed run(`./main $(N) $(nt) $(rt)`)
    println("Time taken to simulate: $(t1+t2) seconds\n")
    t3 = @elapsed plt = generateAnimation(W, H)
    println("Time taken to plot: $(t3) seconds\n")
    plt
end

function plotFirstFrame(W=1200, H=1200)
    data = readdlm("intermediateResults/boundaryData.csv", ',', Float64)
    x, y = data[1, 1:2:end], data[1, 2:2:end]
    plot(x, y, legend=:none, color=:blue,
        seriestype=:scatter,
        xlims=[-W, W], ylims=[-H, H])
end



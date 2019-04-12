using Distributions
using Statistics
using DelimitedFiles
using Random
using LightGraphs
using LightGraphsExtras
using GraphPlot
using Compose

function mutational_neightbors(seq1, base = ['A','U','C','G'])
    mut_neighbors::Array{String,1} = ["" for i in 1:((length(base)-1)*length(seq1))]
    #mut_neighbors::Array{String,1} = []
    #println(mut_neighbors)
    #println(base)
    index = 1
    for i in 1:length(seq1)
        if seq1[i] != 'N'
            for ch in base
                if ch != seq1[i]
                    mut_neighbors[index] = seq1[1:i-1]*ch*seq1[i+1:end]
                    index += 1
                    #push!(mut_neighbors, seq1[1:i-1]*ch*seq1[i+1:end])
                end
            end
        end
    end
    return mut_neighbors
end


function load_data(path, file, shufffle = false)
    arti_rep1 = readdlm(path * file, ',')
    arti_seq = arti_rep1[2:end,1]
    arti_seq = convert(Array{String,1}, arti_seq)

    arti_exp = arti_rep1[2:end, end]
    for i in 1:length(arti_exp)
        if arti_exp[i] == "NA"
            arti_exp[i] = missing
        end
    end
    arti_exp = convert(Array{Union{Float64,Missing},1}, arti_exp)
    rng = MersenneTwister(1234);
    if shufffle == true
        shuffle(rng, arti_exp)
    end
    return arti_seq, arti_exp
end

variants, arti_array = load_data("/home/rljahn/Dropbox/Julia repo/mutational distance/","arti_SDR_union_count25.csv")
#fepb_array = load_data("/home/rljahn/Dropbox/Julia repo/mutational distance/","fepB_SDR_union_count25.csv")
arti_index = [i for i in 1:length(arti_array)]

# yourdic = Dict(zip(keys, vals))
fitness_dict = Dict(zip(variants, arti_array))
name_dict = Dict(zip(arti_index,variants))
index_dict = Dict(zip(variants,arti_index))

# Create Graph
using LightGraphs
using LightGraphsExtras
using SimpleWeightedGraphs
using TikzGraphs

n_nodes = length(arti_array)

g = SimpleWeightedDiGraph(n_nodes)
for i = 1:n_nodes
    mut_neigh = mutational_neightbors(name_dict[i])
    f1 = fitness_dict[name_dict[i]]
    for seq in mut_neigh
        f2 = fitness_dict[seq]
        if (f2 > f1) & (!ismissing(f1)) & (!ismissing(f2))
            add_edge!(g, i, index_dict[seq],f2-f1)
        end
    end
end


nvertices = nv(g) # number of vertices
nedges = ne(g)    # number of edges

g2 = WheelGraph(262144) # ~1000 is fast for display
display(gplot(g2))
gplot(g2, layout=random_layout)
import Cairo, Fontconfig
#using ColorTypes
@time draw(PNG("Wheel_262144.png", 25cm, 25cm), gplot(g2, layout=random_layout))
#enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
#enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
current_start = index_dict["UAAGGAGGU"]
diff_result = diffusion(g, 0.5, 1000, normalize = true, initial_infections = [current_start] ) # sample(vertices(g), 1)
NBrandom_walk_result = non_backtracking_randomwalk(g, 1, 100)
random_walk_result = randomwalk(g, 1, 100)
#nodefillc = [RGBA(0.0,0.8,0.8,i) for i in alphas]

#a_star(g,1,3)
#a_star(g, s, t[, distmx][, heuristic])

#index_dict["GGGGGGGGG"] => 174763
#index_dict["CCCCCCCCC"] => 87382

# plot

#display(gplot(g))
"""
g = WheelGraph(10); am = Matrix(adjacency_matrix(g))
loc_x, loc_y = layout_spring_adj(am)
draw_layout_adj(am, loc_x, loc_y, filename="wheel10.svg")

TikzGraphs.plot(g)

t = TikzGraphs.plot(g)
using TikzPictures # this is required for saving
TikzPictures.save(PDF("graph"), t)
TikzPictures.save(SVG("graph"), t)
TikzPictures.save(TEX("graph"), t)

TikzGraphs.plot(g, ["A", "B", "C", "D"])
"""

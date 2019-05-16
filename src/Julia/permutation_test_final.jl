using Distributions
using Statistics
using DelimitedFiles
using Random
using Plots


path = "/home/rljahn/FitnessLandscape_SD/data/Read Count Table/UnionTable/"

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

function mutational_neightbors(seq1, base = ['A','U','C','G'])::Array{String,1}
    mut_neighbors::Array{String,1} = ["" for i in 1:((length(base)-1)*length(seq1))]
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

function pearson_correlation(x,y)
    ind_x = [!ismissing(i) for i in x]
    ind_y = [!ismissing(i) for i in y]
    ind = ind_x .& ind_y
    x2::Array{Float64,1} = x[ind]
    y2::Array{Float64,1} = y[ind]
    return cor(x2,y2), length(x2),x2,y2
end

arti_seq, arti_array = load_data(path, "arti_SDR_union_count25.csv")
arti_dict = Dict(zip(arti_seq,arti_array))
fepb_seq, fepb_array = load_data(path, "fepB_SDR_union_count25.csv")
fepb_dict = Dict(zip(fepb_seq,fepb_array))
dmsc_seq, dmsc_array = load_data(path, "dmsC_SDR_union_count25.csv")
yfp_seq, yfp_array = load_data("/home/rljahn/FitnessLandscape_SD/data/Read Count Table/", "yfp_SDR_rep1_count25.csv")

function neighbor_center_correlation(seq_list::Array{String,1},exp_list::Array{Union{Missing, Float64},1},show::Bool = true)::Float64
    lib_dict = Dict{String,Union{Missing, Float64}}(zip(seq_list,exp_list))
    center::Array{Float64,1} = []
    neighbor_mean_value_array::Array{Float64,1} = []
    for seq in seq_list
        neighbors = mutational_neightbors(seq)
        if (!ismissing(lib_dict[seq])) & all([!isequal(lib_dict[ne],missing) for ne in neighbors])
            push!(center, lib_dict[seq])
            push!(neighbor_mean_value_array, mean([lib_dict[ne] for ne in neighbors]))
        end
    end
    # Calculate Pearson Correlation
    r = cor(center, neighbor_mean_value_array)
    if show == true
        println("Pearson correlation is ", r)
    end
    return r
end

#@time r_arti = neighbor_center_correlation(arti_seq,arti_array)

function MonteCarlo_permutation_test(lib_seq,lib_exp, sample_size = 1000, filename="")
    r_from_data::Float64 = neighbor_center_correlation(lib_seq,lib_exp)
    println("r from sample is $r_from_data")
    greater_count::Int64 = 0
    r_distribution::Array{Float64,1} = []
    rng = MersenneTwister(1234);
    for i in 1:sample_size
        lib_exp_new = shuffle(rng, lib_exp)
        r = neighbor_center_correlation(lib_seq,lib_exp_new, false)
        push!(r_distribution,r)
        #println(r)
        if r > r_from_data
            greater_count += 1
        end
    end
    println("Write!")
    writedlm("/Users/imac/Desktop/Julia_repo/" * filename * "_" * "r_distribution.csv", r_distribution, ",")
    println("P-value under $sample_size sampling = ", greater_count/sample_size)
    #histogram(r_distribution)
    #savefig("Users/imac/Desktop/Julia_repo/" * filename * "_" * "histogramOFr.svg")
end

function plot_r_distribution(filename = "")
    r_distribution = readdlm("/Users/imac/Desktop/Julia_repo/" * filename * "_" * "r_distribution.csv")
    histogram(r_distribution, legend=false)
    xlabel!("r")
    ylabel!("Counts")
    savefig(filename * "_" * "histogramOFr.svg")
end
# sum(skipmissing([1, missing]))

#plot_neighbor_vs_center(X,Y)


#println(">>>>>> arti >>>>>>")
#@time MonteCarlo_permutation_test(arti_seq,arti_array, 10000, "arti")
#println(">>>>>> dmsC >>>>>>")
#@time MonteCarlo_permutation_test(dmsc_seq,dmsc_array, 10000, "dmsC")
#println(">>>>>> fepB >>>>>>")
#@time MonteCarlo_permutation_test(fepb_seq,fepb_array, 10000, "fepB")
println(">>>>>> yfp >>>>>>")
@time MonteCarlo_permutation_test(yfp_seq,yfp_array, 10000, "yfp")
#plot_r_distribution("arti")
#plot_r_distribution("dmsC")
#plot_r_distribution("fepB")
plot_r_distribution("yfp")
#end

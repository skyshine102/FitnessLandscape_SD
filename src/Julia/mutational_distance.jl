using DelimitedFiles
using Plots
using Random

#using DataFrames
#using CSV

# read files and
# create two arrays to store the sequences and expression values
#arti_rep1 = CSV.read("/home/rljahn/Dropbox/Julia repo/mutational distance/arti_SDR_rep1.csv")


path = "/home/rljahn/FitnessLandscape_SD/data/Read Count Table/UnionTable/"
function theoretical_value(n::Int64)
    result::Float64 = 0.0
    for i::Int64 = 1:(n+1)
        result += (i-1) * binomial(n,i-1) * 3^(n-i+1)
        println(result)
    end
    result = result/(4^n)
    return result
end


function mutational_distance(file, shuffule::Bool = false)
    arti_rep1 = readdlm(path * file, ',')
    arti_seq = arti_rep1[2:end,1]
    arti_seq = convert(Array{String,1}, arti_seq)

    arti_exp = arti_rep1[2:end, end]
    for i in 1:length(arti_exp)
        if arti_exp[i] == "NA"
            arti_exp[i] = NaN
        end
    end
    arti_exp = convert(Array{Float64,1}, arti_exp)
    rng = MersenneTwister(1234);
    if shuffule == true
        println("Shuffled!!!")
        shuffle!(rng, arti_exp)
    end

    #gpmapping = Dict(zip(arti_seq,arti_exp))

    # a function to compute hamming distance
    # e.g. hamming_distance(arti_seq[1],arti_seq[end])
    function hamming_similarity(seq1::String,seq2::String, len::Int64 = 9)::Float64
        return sum([1 for m = 1:len if seq1[m] == seq2[m]])
        #if length(seq1) == length(seq1)
        #    return length(seq1) - sum(collect(seq1) .== collect(seq2))
        #end
    end

    # a function to classify sequences into bins/buckets
    # based on theirs expression values
    function bucket_classifier(array::Array{Float64,1}, mini::Float64 = 0.0, maxi::Float64 = 3.2, number_bins::Int64 = 20)::Array{Union{Nothing, Int64},1}
        function comparator(input_value::Float64, bins = mini:(maxi-mini)/number_bins:maxi )
            for i = 1:length(bins)
                if input_value > bins[i] && input_value <= bins[i+1]
                    return i
                end
            end
        end
        values = comparator.(array)
        return values
    end

    buckets = [[] for i = 1:20]
    bucket_indeces = bucket_classifier(arti_exp)
    for (index, i) in enumerate(bucket_indeces)
        if i != nothing
            push!(buckets[i],arti_seq[index])
        end
    end


    # check the variants in each bins
    println([length(buckets[i]) for i in 1:20])
    println("The first bin has: ", buckets[1])
    println("The last bin has: ", buckets[end])
    #dictionary = Dict(zip(arti_seq,bucket_classifier(arti_exp)))

    n::Int64 = 20
    distance_matrix = zeros(n,n)

    for i::Int64 = 1:n
        for j::Int64 = 1:i
            println("Current progress is ($i,$j)!!!")
            t_start = time()
            dist::Float64 = 0.0
            num::Float64 = 0.0
            if i == j
                for k::Int64 = 1:length(buckets[i])
                    for l::Int64 = 1:(k-1)
                        dist += (hamming_similarity(buckets[i][k],buckets[j][l])*2)
                        num += (1.0*2)
                    end
                    dist += hamming_similarity(buckets[i][k],buckets[j][k])
                    num += 1.0
                end
            else
                for seq1::String in buckets[i]
                    for seq2::String in buckets[j]
                        dist += hamming_similarity(seq1,seq2)
                        num += 1.0
                    end
                end
            end

            #for seq1::String in buckets[i]
            #    for seq2::String in buckets[j]
            #        dist += hamming_similarity(seq1,seq2)
            #        num += 1.0
            #    end
            #end

            println("similarity score = $dist")
            println("number = $num")
            if num != 0
                distance_matrix[i,j] = dist/num
            end
            println("result = $(distance_matrix[i,j])")
            tmp = time()-t_start
            println("Time used: $tmp seconds")
        end
    end

    distance_matrix = distance_matrix + transpose(distance_matrix)
    for i = 1: size(distance_matrix)[1]
        distance_matrix[i,i] /= 2
    end

    heatmap(distance_matrix, c=:RdBu, aspect_ratio =1)
    writedlm("/Users/imac/Desktop/Julia_repo/" * "matrix_"* file * "_20190311.csv", distance_matrix, ',')
    #savefig("/Users/imac/Desktop/Julia_repo/" * "heatmap_" * file * ".svg")
    #println("Maximum time used in each loop: $max_time")
end

println(">>>>>> arti !!!")
@time mutational_distance("arti_SDR_union_count25.csv",true)

println(">>>>>> dmsC !!!")
@time mutational_distance("dmsC_SDR_union_count25.csv",true)

println(">>>>>> fepB !!!")
@time mutational_distance("fepB_SDR_union_count25.csv",true)


function plot_part(file)
    table1 = readdlm("/Users/imac/Desktop/Julia_repo/" * "matrix_"* file * "_20190311.csv",  ',')
    heatmap(table1, c=cgrad([:blue, :green, :yellow, :orange, :red]), clim = (1.5,4), aspect_ratio =1)
    savefig("/Users/imac/Desktop/Julia_repo/" * "heatmap_part" * file * "_20190311.svg")
end

plot_part("arti_SDR_union_count25.csv")
plot_part("fepB_SDR_union_count25.csv")
plot_part("dmsC_SDR_union_count25.csv")

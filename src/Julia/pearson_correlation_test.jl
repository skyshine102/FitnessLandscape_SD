using Distributions
using Statistics
using DelimitedFiles
using Random


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

    return arti_exp

end

arti_array = load_data("/home/rljahn/Dropbox/Julia repo/mutational distance/","arti_SDR_union_count25.csv")
fepb_array = load_data("/home/rljahn/Dropbox/Julia repo/mutational distance/","fepB_SDR_union_count25.csv")

function pearson_correlation(x,y)
    ind_x = [!ismissing(i) for i in x]
    ind_y = [!ismissing(i) for i in y]
    ind = ind_x .& ind_y
    x2::Array{Float64,1} = x[ind]
    y2::Array{Float64,1} = y[ind]
    return cor(x2,y2), length(x2),x2,y2
end


r, sample_size, x2, y2 = pearson_correlation(arti_array,fepb_array)
alpha = 0.05
df = sample_size - 2


if r != 1
    t_star = r * sqrt(df) / sqrt(1 - r^2)
end

d = TDist(df)
cdf(d,t_star)

#a = r * sigma_y / sigma_x
print(0.979788 /(std(y2)/std(x2)))

#quantile.(Normal(), [0.5, 0.95])

#a = Vector{Float64}(1,1,1,1)
#AbstractVector{T}

#0.9999999999999993

using DataFrames, GLM
data = DataFrame(X=x2, Y=y2)
ols = lm(@formula(Y ~1 + X), data)

using StatsBase
r_spearman = corspearman(x2,y2)
t_star_spearman = r_spearman * sqrt(df) / sqrt(1 - r_spearman^2)


function permutation_test(array1::Array{Float64,1},array2::Array{Float64,1}, n_samples::Int64)
    # Compute observed correlation: r_obs
    r_obs::Float64 = cor(array1,array2)
    # Initialize permutation replicates: perm_replicates
    perm_replicates::Array{Float64,1} = zeros(n_samples)
    # Draw replicates
    for i in 1:n_samples
        # Permute illiteracy measurments: illiteracy_permuted
        array1_permuted = shuffle(array1)
        # Compute Pearson correlation
        perm_replicates[i] = cor(array1_permuted, array2)
    end
    # Compute p-value: p
    num_obs_greater_r = sum(perm_replicates.>r_obs)
    println("# of (r_obs > r)", num_obs_greater_r)
    p = num_obs_greater_r/n_samples
    println("p-val =", p)
    return p
end

@time permutation_test(x2,y2,10000) # 25 secs

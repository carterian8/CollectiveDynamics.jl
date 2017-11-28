### ============== ### ============== ### ============== ###
## Compute system expansion through
## mean distance between particles
## (export data)
## Martin Zumaya Hernandez
## 21 / 02 / 2017
### ============== ### ============== ### ============== ###

# addprocs(4)

### ================================== ###

# using Plots
# @everywhere using DistributedArrays

### ================================== ###

# function calc_Rij_2D(pos, Rij)
#
#     N = size(Rij, 1)
#
#     @parallel for i in 1:2:2N, j in (i+2):2:2N
#
#         k = div(i, 2) + 1
#         l = div(j, 2) + 1
#
#         Rij[l, k] = norm([pos[i], pos[i+1]] - [pos[j], pos[j+1]])
#     end
#
# end

### ================================== ###

function make_dir_from_path(path)

    try
        mkdir(path)
    catch error
        println("Main data folder already exists")
    end

end

### ================================== ###

function get_times(Ti, Tf)

    times = [convert(Int, exp10(i)) for i in Ti:Tf]
    tau = Int64[]

    push!(tau, 1)

    for i in 1:(length(times) - 1)

        for t in (times[i]+1):times[i+1]

            if t % times[i] == 0 || t % div(times[i], exp10(1)) == 0
                push!(tau, t)
                # println("//////// ", t)
            end
        end

    end

    return tau
end

### ================================== ###

function calc_Rij_3D(pos::SharedArray, Rij::SharedArray)

    N = size(Rij, 1)

    for i in 1:3:3N

        for j in (i+3):3:3N

            k = div(i, 3) + 1
            l = div(j, 3) + 1

            Rij[l, k] = norm([pos[i], pos[i+1], pos[i+2]] - [pos[j], pos[j+1], pos[j+2]])
        end
    end

end

### ================================== ###

@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,2),nchunks+1)]
    1:size(q,1), splits[idx]+1:splits[idx+1]
end

### ================================== ###

@everywhere function calc_means(pos, vel, mean, nn_mean, vel_mean, N)

    Rij = zeros(N, N)

    j_range = myrange(pos)[2]

    for j in j_range
        # @show j
        calc_Rij_3D(pos[:,j], Rij)

        mean[j] = mean(Symmetric(Rij, :L))
        nn_mean[j] = mean(sort(Symmetric(Rij, :L), 1)[2,:])
        vel_mean[j] = norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]))

    end
end

### ================================== ###

folder = ARGS[1]
N = parse(Int, ARGS[2])
k = ARGS[3]
# w = ARGS[4]


# folder = "NLOC_MET_3D"
# folder = "NLOC_TOP_3D"
# folder = "NLOC_DATA"
# folder = "SVM_GRID_3D"
# folder = "SVM_GRID_FN_2D"

# folder = "NLOC_MET_3D_EXT"
# folder = "new/NLOC_MET_3D_EXT"
# folder = "new/NLOC_TOP_3D_EXT"

# N = 4096
# N = 4000
# N = 1024
# N = 512
# N = 256
# N = 128
# N = 100

# k = "9.0"
# k = "0.5"
# k = "0.001"
# k = "0.0125"
# w = "0.5"

w = "0.5"

Ti = 0
Tf = 6

times = get_times(Ti, Tf)

### ================================== ###

data_folder_path = joinpath(homedir(),"art_DATA",folder,"DATA","data_N_$(N)","data_N_$(N)_k_$(k)_w_$(w)")

output_folder_path = joinpath(homedir(),"art_DATA",folder,"EXP_N")

make_dir_from_path(output_folder_path)
make_dir_from_path(output_folder_path * "/exp_data_N_$(N)")

params = "_k_$(k)_w_$(w)"
println(params)

### ================================== ###

exp_file        = open(output_folder_path * "/exp_data_N_$(N)" * "/exp" * params * ".dat", "w+")
order_file      = open(output_folder_path * "/exp_data_N_$(N)" * "/order" * params * ".dat", "w+")
nn_mean_file    = open(output_folder_path * "/exp_data_N_$(N)" * "/nn_mean" * params * ".dat", "w+")
# vel_mean_file = open(output_folder_path * "/exp_data_N_$(N)" * "/vel_mean" * params * ".dat", "w+")

### ================================== ###

sh_mean     = SharedArray{Float64}(length(times)-1)
sh_nn_mean  = SharedArray{Float64}(length(times)-1)
sh_vel_mean = SharedArray{Float64}(length(times)-1)

reps = [match(r"\w+_(\d+).\w+", x).captures[1] for x in filter(x -> ismatch(r"pos_\d+.\w+", x), readdir(data_folder_path))]

# r = 2
for r in reps

    println("rep=",r)

    raw_data = reinterpret(Float64,read(data_folder_path * "/pos_$(r).dat"))

    if length(raw_data) == 3N*(length(times)-1)
        pos_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    else
        pos_data = reshape(raw_data[3N+1:end], 3N, div(length(raw_data[3N+1:end]), 3N))
    end

    raw_data = reinterpret(Float64,read(data_folder_path * "/vel_$(r).dat"))

    if length(raw_data) == 3N*(length(times)-1)
        vel_data = reshape(raw_data, 3N, div(length(raw_data), 3N))
    else
        vel_data = reshape(raw_data[3N+1:end], 3N, div(length(raw_data[3N+1:end]), 3N))
    end

    sh_pos = SharedArray{Float64,2}(size(pos_data))
    sh_vel = SharedArray{Float64,2}(size(vel_data))

    for i in 1:length(pos_data)
        sh_pos[i] = pos_data[i]
        sh_vel[i] = vel_data[i]
    end

    # d_pos = distribute(reshape(raw_data, 3N, div(length(raw_data), 3N)), procs = workers(), dist = [1,length(workers())])

    @sync begin
        for p in workers()
            @async remotecall_wait(calc_means, p, d_pos, d_vel, d_mean, d_nn_mean, d_vel_mean, N)
        end
    end

    println("finished calc_means")

    write(exp_file, fetch(sh_mean))
    write(nn_mean_file, fetch(sh_nn_mean))
    write(order_file, fetch(sh_vel_mean))

end

close(exp_file)
close(nn_mean_file)

close(order_file)
# close(vel_mean_file)

println("Done")

### ================================== ###

# gui()
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), means, m = :o, xscale = :log10, yscale = :log10, leg = false)
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), nn_means, m = :o, xscale = :log10, yscale = :log10, leg = false)
#
# plot(union(times[1:369],collect(2exp10(5):exp10(5):exp10(6))), [norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N])) for j in 1:size(vel_data, 2)], xscale = :log10, leg = false)

# pos_data
# sh_pos
# sh_vel
#
# @everywhere using DistributedArrays
#
# addprocs(4)
#
# N = 16
#
# t = rand(5, N)
#
# d_t = distribute(t, procs = workers(), dist = [1, length(workers())])
#
# d_s = dzeros(size(d_t,2))
#
# for p in workers()
#     # println(remotecall_fetch(myrange, p, sh_pos))
#     # println(remotecall_fetch(calc_means, p, sh_pos, sh_vel, sh_mean, sh_nn_mean, sh_vel, N))
#     println(p)
#     println(remotecall_fetch(myrange, p, sh_pos))
#     println(remotecall_fetch(myrange, p, sh_vel))
#     println(remotecall_fetch(localindexes, p, sh_mean))
#     println(remotecall_fetch(localindexes, p, sh_nn_mean))
#     println(remotecall_fetch(localindexes, p, sh_vel_mean))
#     println()
# end
#
# @everywhere function d_test(d)
#     println(localindexes(d))
#     println(d[:l])
# end
#
# @everywhere function d_test1(d)
#     for i in localindexes(d)[1]
#         d[i] = 10.
#     end
# end
#
# remotecall_fetch(d_test, 2, d_t)
# remotecall_fetch(d_test1, 2, d_s)
#
# DArray(size(d_s), procs(d_s)) do I
#     println(I)
# end
#
# size(sh_pos)
#
# @everywhere function myrange(q::SharedArray)
#     idx = indexpids(q)
#     if idx == 0 # This worker is not assigned a piece
#         return 1:0, 1:0
#     end
#     nchunks = length(procs(q))
#     splits = [round(Int, s) for s in linspace(0,size(q,2),nchunks+1)]
#     1:size(q,1), splits[idx]+1:splits[idx+1]
# end
#
#
# @everywhere function calc_means(pos, vel, mean, nn_mean, vel_mean, N)
#
#     Rij = zeros(N, N)
#
#     j_range = myrange(pos)[2]
#
#     for j in j_range
#         # @show j
#         calc_Rij_3D(pos[:,j], Rij)
#
#         mean[j] = mean(Symmetric(Rij, :L))
#         nn_mean[j] = mean(sort(Symmetric(Rij, :L), 1)[2,:])
#         vel_mean[j] = norm(mean([[vel_data[i, j], vel_data[i+1, j], vel_data[i+2, j]] for i in 1:3:3N]))
#
#     end
# end
#
# sh_mean
#
# println(remotecall_fetch(calc_means, 7, sh_pos, sh_vel, sh_mean, sh_nn_mean, sh_vel, N))
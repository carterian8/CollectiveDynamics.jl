### ============== ### ============== ###
##    Simple Vicksek Model  2D         ##
##    Martin Zumaya Hernandez          ##
##    EXAMPLE SIMULATION SCRIPT        ##
### ============== ### ============== ###

### ============ INCLUDE PACKAGES ============ ###

using CollectiveDynamics.SVM2D
using Plots

### ============ SYSTEM'S PARAMETERS ============ ###

N   = parse(Int, ARGS[1]) # number of particles
ρ   = parse(Float64, ARGS[2]) # density
η   = parse(Float64, ARGS[3]) # noise intensity
ω   = parse(Float64, ARGS[4]) # maximal angular velocity
φ   = parse(Float64, ARGS[5]) # particle field of view
T   = parse(Int, ARGS[6]) # integration time steps
rep = parse(Int, ARGS[7]) # ensemble index
viz = true # visualize steps

### =============== ### =============== ###

v0 = 0.2 # 1.0 # particle's speed
dt = 1.0 # integration timestep

### =============== ### =============== ###

M = 10 # number of boxes per dimension

# l = 0.5 # interaction range is double the distance each particle moves in one time step
#l = 0.1 # interaction range is ten times the distance each particle moves in one time step

# r0 = (v0 * dt) / l # interaction range
ratio = 1.03
# r0 = v0 / (ratio * ω) # interaction range
r0 = 1.0

L = sqrt(N / ρ) # size of box

### =============== ### =============== ###

cell_size = step(range(0., stop=L , length=M))

### ============ SYSTEM'S INITIALIZATION ============ ###

box = Box(L, M) # suqare box of size L and M partitions per dimension
flock = Flock(N, L, dt, v0, r0, η, ω, φ) # flock initialization

### ============== ### ============== ###

output_path = set_output_data_structure_vsk("SVM_2D", N, ρ)

pos_file = open(output_path * "/pos_$(rep).dat", "w+")
vel_file = open(output_path * "/vel_$(rep).dat", "w+")

times = [convert(Int, exp10(i)) for i in 0:T]

plt = !viz ? nothing : quiver(
    [ flock.pos[i][1] for i in 1:N ],
    [ flock.pos[i][2] for i in 1:N ],
    quiver=(
        [ flock.vel[i][1] for i in 1:N ],
        [ flock.vel[i][2] for i in 1:N ]
    ),
    legend=false
)
savefig(plt, output_path * "/step0.png")

### ============ TIME EVOLUTION ============ ###

for i in 1:(length(times) - 1)

    if i > 1

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, box, cell_size)

            if t % times[i] == 0 || t % times[i-1] == 0
                println("//////// ", t)
                if viz
                    plt = quiver(
                        [ flock.pos[i][1] for i in 1:N ],
                        [ flock.pos[i][2] for i in 1:N ],
                        quiver=(
                            [ flock.vel[i][1] for i in 1:N ],
                            [ flock.vel[i][2] for i in 1:N ]
                        ),
                        legend=false
                    )
                    savefig(plt, output_path * "/step$(t).png")
                end
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    else

        for t in (times[i]+1):times[i+1]

            evolve_system(flock, box, cell_size)

            if t % times[i] == 0
                println("//////// ", t)
                if viz
                    plt = quiver(
                        [ flock.pos[i][1] for i in 1:N ],
                        [ flock.pos[i][2] for i in 1:N ],
                        quiver=(
                            [ flock.vel[i][1] for i in 1:N ],
                            [ flock.vel[i][2] for i in 1:N ]
                        ),
                        legend=false
                    )
                    savefig(plt, output_path * "/step$(t).png")
                end
                write(pos_file, vcat(flock.pos...))
                write(vel_file, vcat(flock.vel...))
            end
        end

    end

end

close(pos_file)
close(vel_file)

println("Done all")

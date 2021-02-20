### ============== ### ============== ###
### SIMPLE VICSEK MODEL 2D
### PERIODIC BOUNDARY CONDITION
### ============== ### ============== ###

function set_output_data_structure_vsk(path, N, ρ)

    parent_folder_path = "$(homedir())/art_DATA/$(path)"
    folder_path        = parent_folder_path * "/DATA/data_N_$(N)"
    reps_path          = folder_path * "/data_N_$(N)_rho_$(ρ)"

    try
        mkdir("$(homedir())/art_DATA")
    catch error
        println("Main data folder already exists")
    end

    try
        mkdir(parent_folder_path)
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(parent_folder_path * "/DATA")
    catch error
        println("Parent folder already exists")
    end

    try
        mkdir(folder_path)
    catch error
        println("Folder already exists")
    end

    try
        mkdir(reps_path)
    catch error
        println("Parameter folder already exists")
    end

    return reps_path
end

### ============== ### ============== ###
"""
    Box(L, M)
Simulation Box
# Constructor Arguments
* L -> size of box
* M -> number of cells per dimension

# Fields
* L -> linear size
* M -> number of cells per dimension

* bulk_cells   -> discrete coordinates
* bottom_cells -> discrete coordinates
* top_cells    -> discrete coordinates
* left_wall    -> discrete coordinates
* right_wall   -> discrete coordinates
* front_wall   -> discrete coordinates
* back_wall   -> discrete coordinates
* corners -> discrete coordinates

* center_cell -> real coordinates of cell center
* vel_cell -> average velocity of particles within cell

* p_per_cell -> Dict: cell_id => [particles id within cell]
"""
mutable struct Box
    L :: Float64
    M :: Float64

    bulk_cells   :: Array{Array{Int64,1},2}
    bottom_cells :: Array{Array{Int64,1},1}
    top_cells    :: Array{Array{Int64,1},1}
    left_cells    :: Array{Array{Int64,1},1}
    right_cells   :: Array{Array{Int64,1},1}
    corners   :: Array{Array{Int64,1},2}

    p_per_cell :: Dict{Array{Int,1}, Array{Int, 1}}

    vel_cell :: Array{Array{Float64,1},2}
    center_cell :: Array{Array{Float64,1},2}

    function Box(L, M)

        ### ============== ### ============== ###
        ## indices of cells

        bulk_cells = [[i,j] for i in 2:M-1, j in 2:M-1]

        bottom_cells = [[i,M] for i in 2:M-1]
        top_cells = [[i,1] for i in 2:M-1]

        left_cells = [[1,j] for j in 2:M-1]
        right_cells = [[M,j] for j in 2:M-1]

        corners = [[i,j] for i in [1,M], j in [1,M]]

        ### ============== ### ============== ###

        # particles id's in each cell
        p_per_cell = Dict{Array{Int,1}, Array{Int, 1}}()

        # average velocity per cell
        vel_cell = Array{Array{Float64,1}}(undef,M,M)

        # Initialize velocity vector per cell and compute center of each cell
        for i in 1:length(vel_cell)
            vel_cell[i] = zeros(Float64, 2)
        end

        # position of center of each box
        center_cell = Array{Array{Float64,1}}(undef,M,M)

        r = range(0., stop=L, length=M+1)
        mp = collect(Float64, r[1:length(r) - 1] .+ 0.5 * step(r))

        for i in 1:M, j in 1:M
            # println(i,"\t", j, "\t", k)
            center_cell[i,j] = [mp[i], mp[j]]
        end

        new(L, float(M), bulk_cells, bottom_cells, top_cells, left_cells, right_cells, corners, p_per_cell, vel_cell, center_cell)
    end
end

### ============== ### ============== ###
"""
    Flock(N, L, v0, p, d)
Inertial Flock type
# Constructor Arguments
* N -> number of particles
* L -> size of box
* dt -> integration time step
* v0 -> particles speed
* r0 -> local interaction range
* η -> noise intensity
* ω -> maximal angular velocity
* φ -> particle field of view

# Fields
* N -> number of particles
* v0 -> particles speed
* r0 -> local interaction range
* dt -> integration time step
* η -> noise intensity
* ω -> maximal angular velocity
* φ -> particle field of view

* pos -> particles' positions
* vel -> particles' velocities
* r -> particles' interaction ranges
* v_r -> interaction vector
* k_r -> local connectivity

* p_cell_id -> particle's cell id
"""
mutable struct Flock
    N  :: Int64
    v0 :: Float64
    r0 :: Float64
    dt :: Float64
    η :: Float64
    ω :: Float64
    φ :: Float64

    pos :: Array{Array{Float64, 1}}
    vel :: Array{Array{Float64, 1}}
    r :: Array{Array{Float64, 1}}
    v_r :: Array{Array{Float64, 1}}
    k_r :: Array{Float64, 1}

    p_cell_id :: Array{Array{Int64,1}}

    function Flock(N, L, dt, v0, r0, η, ω, φ)

        # particles' positions and velocities
        pos = [ [L*rand(), L*rand()] for i in 1:N ]
        vel = v0 * [ normalize([2*rand() - 1, 2*rand() - 1]) for i in 1:N ]
        # vel = v0 * [ones(2) for i in 1:N]

        # interaction ranges
        r = r0 * [ones(2) for i in 1:N]

        # interaction vector
        v_r = [zeros(2) for i in 1:N]

        # connectivity per particle
        k_r = zeros(Float64, N)

        # box index for each particle in each dimension
        p_cell_id = [zeros(Int,2) for i in 1:N]

        new(N, v0, r0, dt, η, deg2rad(ω), deg2rad(φ), pos, vel, r, v_r, k_r, p_cell_id)
    end
end

## ============== ### ============== ###
function av_vel_cell(box, vel)
    for k in keys(box.p_per_cell)
        box.vel_cell[k[1], k[2]] = mean([vel[i] for i in box.p_per_cell[k]])
    end

end

### ============== ### ============== ###
function assign_cell(flock, box, cell_size)

    box.p_per_cell = Dict{Array{Int,1}, Array{Int, 1}}()

    for i in 1:flock.N
        flock.p_cell_id[i] = convert(Array{Int}, div.(floor.(flock.pos[i]), cell_size)) .+ 1

        haskey(box.p_per_cell, flock.p_cell_id[i]) ? push!(box.p_per_cell[flock.p_cell_id[i]], i) : box.p_per_cell[flock.p_cell_id[i]] = [i]
    end

end

### ============== ### ============== ###
function check_bulk_cells(flock, p_id, cell_id, p_per_cell)

    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1

    # println("center cell: ", cell_id)

    for i in -1:1, j in -1:1
        key = [cell_id[1]+i, cell_id[2]+j]
        k_t[c], v_t[c] = part_bulk_cell_interaction(p_id, flock, get!(p_per_cell, key, []))
        c += 1
        # println(key, "\t", get!(p_per_cell, key, []))
    end

    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]
    
    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###
function check_corners(flock, p_id, cell_id, p_per_cell, L, M)

    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1
    # println("corner cell: ", cell_id)

    if cell_id[1] == 1 && cell_id[2] == 1
        for i in [M, 1, 2], j in [M, 1, 2]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == 1 && cell_id[2] == M
        for i in [M, 1, 2], j in [M-1, M, 1]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == 1
        for i in [M-1, M, 1], j in [M, 1, 2]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    elseif cell_id[1] == M && cell_id[2] == M
        for i in [M-1, M, 1], j in [M-1, M, 1]
            key = [i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    end

    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]

    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###
function check_bottom_cells(flock, p_id, cell_id, p_per_cell, L, M)

    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1
    # println("center cell: ", cell_id)

    # if in(cell_id[1], 2:M-1) && cell_id[2] == M
        for i in -1:1, j in [M-1, M, 1]
            key = [cell_id[1] + i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end
    
    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]

    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###
function check_top_cells(flock, p_id, cell_id, p_per_cell, L, M)
    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1
    # println("center cell: ", cell_id)

    # if in(cell_id[1], 2:M-1) && cell_id[2] == 1
        for i in -1:1, j in [M, 1, 2]
            key = [cell_id[1] + i, j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]

    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###
function check_left_cells(flock, p_id, cell_id, p_per_cell, L, M)

    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1
    # println("center cell: ", cell_id)

    # if cell_id[1] == 1 && in(cell_id[2], 2:M-1)
        for i in [M, 1, 2], j in -1:1
            key = [i, cell_id[2] + j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]

    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###
function check_right_cells(flock, p_id, cell_id, p_per_cell, L, M)

    k_t = Vector{Float64}(undef,10)
    v_t = Vector{Vector{Float64}}(undef,10)

    c = 1
    # println("center cell: ", cell_id)

    # if cell_id[1] == M && in(cell_id[2], 2:M-1)
        for i in [M-1, M, 1], j in -1:1
            key = [i, cell_id[2] + j]
            k_t[c], v_t[c] = part_bound_cell_interaction(p_id, flock, get!(p_per_cell, key, []), L)
            c += 1
            # println(key, "\t", get!(p_per_cell, key, []))
        end
    # end

    # Include current particle in interactions
    k_t[c] = 1
    v_t[c] = flock.vel[p_id]

    return calc_vt(flock, p_id, k_t, v_t)
end

### ============== ### ============== ###

function get_neighbors(p_id, flock, box)
    if flock.p_cell_id[p_id] in box.bulk_cells
        # println("bulk_cells")
        flock.v_r[p_id] = check_bulk_cells(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell)
    elseif flock.p_cell_id[p_id] in box.bottom_cells
        # println("bottom_cells")
        flock.v_r[p_id] = check_bottom_cells(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L, box.M)
    elseif flock.p_cell_id[p_id] in box.top_cells
        # println("top_cells")
        flock.v_r[p_id] = check_top_cells(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L, box.M)
    elseif flock.p_cell_id[p_id] in box.left_cells
        # println("left_cells")
        flock.v_r[p_id] = check_left_cells(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L, box.M)
    elseif flock.p_cell_id[p_id] in box.right_cells
        # println("right_cells")
        flock.v_r[p_id] = check_right_cells(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L, box.M)
    elseif flock.p_cell_id[p_id] in box.corners
        # println("corners")
        flock.v_r[p_id] = check_corners(flock, p_id, flock.p_cell_id[p_id], box.p_per_cell, box.L, box.M)
    end
end

### ============== ### ============== ###
function part_bulk_cell_interaction(p_id, flock, cell_parts)

    k = 0.
    v_r = zeros(Float64, 2)
    half_φ = flock.φ * 0.5

    for j in cell_parts
        # Make jth point relative to point at p_id
        rel_pos_j = flock.pos[j] - flock.pos[p_id]

        # Do a distance and field of view check
        if norm(rel_pos_j) > 0. && norm(rel_pos_j) <= flock.r0
            # Calculate angle between relative jth point and velocity vector at point p_id
            angle_i_j = angle_btwn_vecs(flock.vel[p_id],rel_pos_j)
            # println("part_bulk_cell: Within r0")
            if angle_i_j <= half_φ
                # println("part_bulk_cell: Within FOV")
                k += 1.
                v_r = v_r + flock.vel[j]
            end
        end
    end

    return k, v_r
end

### ============== ### ============== ###
function part_bound_cell_interaction(p_id, flock, cell_parts, L)

    k = 0.
    v_r = zeros(Float64, 2)
    d_v = zeros(Float64, 2)
    half_φ = flock.φ * 0.5

    for j in cell_parts

        # compute relative vector  taking into account PBC
        for i in eachindex(flock.pos[p_id])
            δ = flock.pos[j][i] - flock.pos[p_id][i]
            δ > 0.5 * L ? d_v[i] = δ - L : d_v[i] = δ
        end

        # check if particles are within interaction range and field of view
        if norm(d_v) > 0. && norm(d_v) <= flock.r0
            # Calculate angle between relative jth point and velocity vector at point p_id
            angle_i_j = angle_btwn_vecs(flock.vel[p_id], d_v)
            # println("part_bound_cell: Within r0")
            if angle_i_j <= half_φ
                # println("part_bound_cell: Within FOV")
                k += 1.
                v_r = v_r + flock.vel[j]
            end
        end
    end

    return k, v_r
end

### ============== ### ============== ###

function calc_vt(flock, p_id, k_t, v_t)

    if sum(k_t) <= 0.
        println("kt was 0")
        return zeros(2)
    end
    
    # Get polar coords of current velocity...
    Θi = atan(flock.vel[p_id][2], flock.vel[p_id][1])
    if Θi < 0
        Θi = Θi + (2*pi)
    end
    ri = sqrt(
        (flock.vel[p_id][1] * flock.vel[p_id][1]) +
        (flock.vel[p_id][2] * flock.vel[p_id][2])
    )

    # Calc difference in orientation between the current orientation and the 
    # average orientation of the interacting particles...
    avg_vt = sum(v_t) ./ sum(k_t)
    avg_Θ = atan(avg_vt[2], avg_vt[1])
    if avg_Θ < 0
        avg_Θ = avg_Θ + (2*pi)
    end
    big_Θ = neg_pi_pos_pi_rads(avg_Θ - Θi)

    # Calc the new angle...
    new_Θi = 0
    if abs(big_Θ) < flock.ω
        new_Θi = avg_Θ
    elseif big_Θ >= flock.ω
        new_Θi = Θi + flock.ω
    elseif big_Θ <= -flock.ω
        new_Θi = Θi - flock.ω
    end

    # Ensure its [0, 360)
    if new_Θi < 0 
        new_Θi = new_Θi + 360
    end

    # Convert back to cartesian coords...
    new_vt = zeros(2)
    new_vt[1] = ri .* cos.(new_Θi)
    new_vt[2] = ri .* sin.(new_Θi)

    return new_vt
    
end

### ============== ### ============== ###

function neg_pi_pos_pi_rads(angle_rads)
    angle_degs = rad2deg(angle_rads)
    int_degs = floor(angle_degs)
    fractional_degs = angle_degs - int_degs

    if (int_degs > 360)
        # Reduce the angle
        int_degs = int_degs % 360
    end

    # Force 0 <= angle < 360
    int_degs = (int_degs + 360) % 360
    
    norm_angle = int_degs + fractional_degs;
    if norm_angle > 180
        norm_angle -= 360
    end

    return deg2rad(norm_angle)

end

### ============== ### ============== ###

function angle_btwn_vecs(a, b)
    tmp = dot(a,b)/(norm(a)*norm(b))
    if abs(tmp) > 1.0
        # this can happen with floating point error
        tmp = sign(tmp)
    end
    return acos(tmp)
end

### ============== ### ============== ###
function update_part(pos, vel, v_r, η, L)

    prop_angle = atan(vel[2], vel[1])
    signal_angle = 0.0
    u_vel = zeros(Float64, 2)

    v_r[1] != 0.0 || v_r[2] != 0.0 ? signal_angle = atan(v_r[2], v_r[1]) - prop_angle : signal_angle = 0.0

    total_angle = signal_angle + η * (2.0 * rand() * pi - pi);

    c = cos(total_angle)
    s = sin(total_angle)

    u_vel[1] = vel[1]*c - vel[2]*s;
    u_vel[2] = vel[1]*s + vel[2]*c;

    for i in eachindex(vel)

        vel[i] = u_vel[i]
        pos[i] += u_vel[i]

        # if(pos[i] > L * 0.5) pos[i] -= L end
        # if(pos[i] <= -L * 0.5) pos[i] += L end

        if(pos[i] > L ) pos[i] -= L end
        if(pos[i] <= 0.0) pos[i] += L end
    end

end

### ============== ### ============== ###
###          SYSTEM EVOLUTION         ###
### ============== ### ============== ###

function evolve_system(flock, box, cell_size)

    # compute each particle cell_id and assgin particles id to each box
    assign_cell(flock, box, cell_size)

    # for i in eachindex(flock.v_r)
    #     println(i, "\t",flock.vel[i])
    # end

    # Search interactions in adjacent cells
    for i in eachindex(flock.pos)
        get_neighbors(i, flock, box)
        # println("pass")
    end

    # for i in eachindex(flock.v_r)
    #     println(i, "\t",flock.v_r[i])
    # end

    # update particles' position
    broadcast(update_part, flock.pos, flock.vel, flock.v_r, flock.η, box.L)

end

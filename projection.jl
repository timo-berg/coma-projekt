using Random, Distributions, Images, Sixel

function abbild(x, y, z)
    p_vec = collect((x,y,z))
    pointOnPlane = p_vec * 250/p_vec[3]
    if -250 < pointOnPlane[1] < 250 && -250 < pointOnPlane[2] < 250
        (floor(Int, pointOnPlane[1]), floor(Int, pointOnPlane[2]))
    else
        nothing
    end
end


function is_visible(p, m, r)
    α, β, γ = p
    x, y, z = m

    # p does not lie in the image
    if isnothing(abbild(α, β, γ))
        false
    end

    a = α^2+β^2+γ^2
    b = -2*(α*x+β*y+γ*z)
    c = x^2+y^2+z^2-r^2

    # Root of negative number not defined
    try
        sqrt(b^2-4*a*c)
    catch exception
        if isa(exception, DomainError)
            false
        end
    end

    # If root is zero, it lies on the visibel part of the sphere (only one 
    # intersect) if not, it has two intersects and is therefore on the backside
    if sqrt(b^2-4*a*c) == 0
        true
    else
        false
    end
end

function samples(x, y, b, h, m, r, dichte)
    δ = rand(Uniform(0,1), dichte, 2)
    sample_array = Array{Tuple{Float64,Float64,Float64},1}(undef,dichte)

    for i in 1:dichte
        θ = (x+δ[i,1])*π/h
        ϕ = (y+δ[i,2])*2*π/b
        sample_array[i] = sphere_projection(θ, ϕ, m, r)
    end

    sample_array
end

function sphere_projection(θ, ϕ, m, r)
    x, y, z = m
    x+r*sin(θ)*cos(ϕ), y+r*sin(θ)*sin(ϕ), z+r*cos(θ)
end

function snapshot_sphere(b, h, daten, m, r, dichte)
    # Projected plane: 500x500 array of RGBA-tuple
    bildebene = Array{Tuple{Float64,Float64,Float64,Float64}}(undef, 500, 500)

    # Iterate over pixel
    for x in 1:b
        for y in 1:h
            # Generate samples and select only visible ones
            sample_array = samples(x, y, b, h, m, r, dichte)
            selection_array = map((x) -> is_visible(x, m, r), sample_array)
            sample_array = sample_array[selection_array]
            # Write visible samples in bildebene. Average multiple entries
            for p in 1:size(sample_array)[1]
                # Get sample
                x_sphere, y_sphere, z_sphere = sample_array[p]
                # Project sample onto the plane
                x_plane, y_plane = abbild(x_sphere, y_sphere, z_sphere)
                # Write pixel value at the projected place
                bildebene[x_plane+250, y_plane+250] = daten[x*y]
            end
        end
    end

    bildebene
end


##################-Test-Section-##########################
# using ImageInTerminal
# include("projection.jl")
# b = 320
# h = 320
# m = (0, 0, 100)
# r = 50
# dichte = 1
# img = load("test.png")
# img_rgba = map((x) -> convert(RGBA, x), img)
# daten = map((color) -> (color.r, color.g, color.b, color.alpha), img_rgba)
# projected_image_tuple = snapshot_sphere(b, h, daten, m, r, dichte)
# projected_image = map((x) -> RGBA{N0f8}(x[1],x[2],x[3],x[4]), projected_image_tuple)
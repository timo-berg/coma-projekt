using Base: Float64
using Random, Distributions, Images, ImageView, Plots

function abbild(p)
    if p[3] < 250
        return nothing
    end

    pointOnPlane = collect(p) * 250/p[3]
    if -250 <= pointOnPlane[1]* 250/p[3] < 250 && -250 <= pointOnPlane[2]* 250/p[3] < 250
        return (floor(Int, pointOnPlane[1]), floor(Int, pointOnPlane[2]))
    else
        return nothing
    end
end

function is_visible(p, m, r)
    α, β, γ = p
    x, y, z = m

    # p does not lie in the image
    if isnothing(abbild(p))
        return false
    end

    a = α^2+β^2+γ^2
    b = -2*(α*x+β*y+γ*z)
    c = x^2+y^2+z^2-r^2

    # Root of negative number not defined
    try
        sqrt(b^2-4*a*c)
    catch exception
        if isa(exception, DomainError)
            return false
        end
    end

    t = 250/p[3]

    s_1 = (-b+sqrt(b^2-4*a*c))/(2*a)
    s_2 = (-b-sqrt(b^2-4*a*c))/(2*a)

    intersect1 = 1 >= round(s_1, digits=3) > t
    intersect2 = 1 >= round(s_2, digits=3) > t

    # If root is zero, it lies on the visibel part of the sphere (only one
    # intersect) if not, it has two intersects and is therefore on the backside
    if intersect1 && intersect2
        return false
    elseif intersect1
        return true
    elseif intersect2
        return true
    end

    false
end

function samples(x, y, b, h, m, r, dichte)
    δ = rand(Uniform(-1,1), dichte, 2)
    # δ = zeros(dichte, 2)
    sample_array = Array{Tuple{Float64,Float64,Float64},1}(undef,dichte)

    for i in 1:dichte
        θ = (x+δ[i,1])*(π/h)
        ϕ = (y+δ[i,2])*2*(π/b)
        sample_array[i] = sphere_projection(θ, ϕ, m, r)
    end

    return sample_array
end

function sphere_projection(θ, ϕ, m, r)
    x, y, z = m
    x-r*cos(θ), y-r*sin(θ)*sin(ϕ), z+r*sin(θ)*cos(ϕ)
end

function snapshot_sphere(b, h, daten, m, r, dichte)
    # Projected plane: 500x500 array of RGBA-tuple
    bild_ebene = Array{Tuple{Float64,Float64,Float64,Float64}}(undef, 500, 500)
    bild_counter = zeros(Int, 500, 500)
    # Iterate over pixel
    for l in 1:(b*h)
        y = l % b
        x = floor(Int, l // b)

        # Generate samples and select only visible ones
        sample_array = samples(x, y, b, h, m, r, dichte)
        selection_array = map((x) -> is_visible(x, m, r), sample_array)
        sample_array = sample_array[selection_array]

        # Write visible samples in bildebene. Average multiple entries
        for p in 1:size(sample_array)[1]
            # Project sample onto the plane
            x_plane, y_plane = abbild(sample_array[p])
            # Write pixel value at the projected place
            bild_counter[x_plane+249, y_plane+249] += 1
            old_value = collect(bild_ebene[x_plane+249, y_plane+249])
            new_value = old_value + (collect(daten[l])-old_value)/bild_counter[x_plane+249, y_plane+249]
            bild_ebene[x_plane+249, y_plane+249] = tuple(new_value...)
        end
    end

    return bild_ebene
end


function sample_z(b,h,m,r,dichte)
    z_values = zeros(Float64, h, b)
    for x in 1:h
        for y in 1:b
            z_values[x,y] = samples(x,y,b,h,m,r,1)[1][3]
        end
    end

    z_values
end

function sample_visible(b,h,m,r,dichte)
    visible = zeros(Bool, h, b)
    for x in 1:h
        for y in 1:b
            visible[x,y] = is_visible(tuple(samples(x,y,b,h,m,r,1)[1]...),m,r)
        end
    end

    visible
end


##################-Test-Section-##########################

# img = transpose(load("textures/ice.jpg"));
# b, h = size(img)
# r = floor(b/1.5)
# m = (0, 0, 260+2*r)
# dichte = 10

# img_rgba = map((x) -> convert(RGBA, x), img);
# daten = map((color) -> (color.r, color.g, color.b, color.alpha), img_rgba);
# projected_image_tuple = snapshot_sphere(b, h, daten, m, r, dichte);
# projected_image = map((x) -> RGBA{N0f8}(x[1],x[2],x[3],x[4]), projected_image_tuple)

# ImageView.imshow(projected_image)

# save("output/proj_checker.png", sum_image)

# z_values = sample_z(b,h,m,r,dichte)
# heatmap(z_values)

# visible = sample_visible(b,h,m,r,dichte)
# heatmap(visible)



# is_visible(samples(Int(b/2),Int(h/2),b,h,m,r,dichte)[1], m, r)
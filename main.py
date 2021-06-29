from julia import Main
from PIL import Image
from PIL import ImageShow
import numpy as np
Main.include("projection.jl")

# Load image and convert it to RGBA
imgobj = Image.open('textures/ice.jpg')
pixels = imgobj.convert('RGBA')
data = np.array(pixels)

# Convert RGBA-arrays to RGBA-tuple
tuple_data = [tuple(data[x,y,:]) for x in range(data.shape[0]) for y in range(data.shape[1])]

# Set image and projection parameter
h = data.shape[0]
b = data.shape[1]
r = b//1.5
m = (0, 0, 260+2*r)
dichte = 10

# Project image and convert float RGBA to int RGBA
projected_array = Main.snapshot_sphere(b, h, tuple_data, m, r, dichte)
projected_array_int = [tuple(map(int, item)) for sublist in projected_array for item in sublist] 

# Plot image
projected_image = Image.new('RGBA', [500,500])
projected_image.putdata(projected_array_int)
ImageShow.show(projected_image)
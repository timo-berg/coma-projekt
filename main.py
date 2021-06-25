from julia import Main
from PIL import Image
from PIL import ImageShow
import numpy as np

imgobj = Image.open('textures/test.png')
pixels = imgobj.convert('RGBA')
data = np.array(pixels)

tuple_data = [tuple(data[x,y,:]) for x in range(data.shape[0]) for y in range(data.shape[1])]


Main.include("projection.jl")


b = 320
h = 320
r = 260
m = (0, 0, 260+r)
dichte = 50

projected_array = Main.snapshot_sphere(b, h, tuple_data, m, r, dichte)
projected_array_int = [tuple(map(int, item)) for sublist in projected_array for item in sublist] 

projected_image = Image.new('RGBA', [500,500])
projected_image.putdata(projected_array_int)
ImageShow.show(projected_image)
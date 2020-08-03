from cmath import exp, nan, isnan, pi, isclose
from math import sqrt, sin, cos
from itertools import count, product
from colour import Color
from PIL import Image
import numpy as np
from types import SimpleNamespace
import tqdm
from newton_iterate import do_newton_iterate


CRITICAL_DIST = 1e-7
ACCEPTABLE_ZERO = 1e-3
MAXITER = 30

POINTS = (tuple(x+y*1j for x, y in product([20.0,35.0],[0.0,37.0])) +
          tuple(x+y*1j for x, y in product([56.0,79.0],[23.0,37.0])) +
          tuple(x+y*1j for x, y in product([110.0,129.0],[0.0,33.0])) +
          tuple(x+y*1j for x, y in product([18.0,32.0],[60.0,75.0])) +
          tuple(x+y*1j for x, y in product([72.0,105.0],[56.0,75.0])) +
          tuple(x+y*1j for x, y in product([124.0,142.0],[42.0,75.0])))

def average(data):
    if (length := len(data)) > 0 :
        return sum(data)/length

def get_centres(pointset, n=4):
    while pointset:
        yield average(pointset[:n])
        pointset = pointset[n:]


#XRES = 153
#YRES = 81

# XRES = 1366
# YRES = 768
NANGLES = 1000
XRES = 820
YRES = 820

#XLIMITS = (0.0, 153.0)
#YLIMITS = (0.0, 81.0)


XLIMITS = (-2.4, 2.4)
YLIMITS = (-2.4, 2.4)


def get_geometry():
    _xdelta = abs(XLIMITS[1] - XLIMITS[0])
    _ydelta = abs(YLIMITS[1] - YLIMITS[0])

    _factor = 1.0/max(_xdelta, _ydelta)

    xmax = _xdelta * _factor
    ymax = _ydelta * _factor

    return xmax, XRES, ymax, YRES

def gen_data(angle):
    _xdelta = abs(XLIMITS[1] - XLIMITS[0])
    _ydelta = abs(YLIMITS[1] - YLIMITS[0])

    _factor = 1.0/max(_xdelta, _ydelta)

    N = 7
    amp1 = 1.0 + 0.5*sin(angle)
    amp2 = 1.0 + 0.5*cos(angle)
    poles = [exp(k*pi/N*1j) for k in range(2*N)]
    roots = ([exp(12*k*pi/N*1j+angle*1j)*amp1 for k in range(N)] +
             [exp(12*k*pi/N*1j-angle*1j)*amp2 for k in range(N)])


    LL_CORNER = XLIMITS[0] + YLIMITS[0]*1j
    poles = [(z - LL_CORNER) * _factor for z in poles]
    roots = [(z - LL_CORNER) * _factor for z in roots]

    colors = {
        roots[0]: Color("DarkRed"),
        roots[1]: Color("RoyalBlue"),
        roots[2]: Color("Tomato"),
        roots[3]: Color("Turquoise"),
        roots[4]: Color("DarkGreen"),
        roots[5]: Color("BlueViolet"),
        roots[6]: Color("green"),
        roots[7]: Color("LawnGreen"),
        roots[8]: Color("blue"),
        roots[9]: Color("RoyalBlue"),
        roots[10]: Color("Turquoise"),
        roots[11]: Color("SteelBlue"),
        roots[12]: Color("Olive"),
        roots[13]: Color("PaleGreen"),
        None: Color("black")
    }
    return roots, poles, colors


def affine(P1, P2, x):
    x1, y1 = P1
    x2, y2 = P2
    return y1 + (y2-y1)/(x2-x1)*(x-x1)

def colorize(colors, info):
    if info['z'] not in colors:
        candidate = min([k for k in colors.keys() if k is not None],
                        key=lambda x: abs(x-info['z']))
        assert isclose(info['z'], candidate), f"{info['z']}, {candidate}"
    the_color = colors[info['z']]
    lum_formula = lambda iterations: min(affine((MAXITER, 0),
                                                (1, 0.5),
                                                iterations), 0.5)
    the_color.luminance = lum_formula(info['niter'])
    return the_color.get_rgb()

def main():

    xmax, xres, ymax, yres = get_geometry()
    xlinspace = np.linspace(0, xmax, xres)
    ylinspace = np.linspace(0, ymax, yres)
    base_array = [x+1j*y for y in ylinspace for x in xlinspace]

    for i, angle in tqdm.tqdm(list(enumerate(np.linspace(0, 2*pi, NANGLES)))[:-1]):
        roots, poles, colors = gen_data(angle)


        output = f"out/newton-anim-{i:06d}.png"
        tqdm.tqdm.write(f"Computing iterations:, angle={angle*180/pi:04.1f}:")
        array = [do_newton_iterate(MAXITER, z, roots, poles) for z in tqdm.tqdm(base_array)]
        tqdm.tqdm.write("Coloring...")
        array = [colorize(colors, z) for z in tqdm.tqdm(array)]
        array = np.array(array)
        array = (array * 255).astype(np.uint8)
        array = array.reshape(yres, xres, 3)
        image = Image.fromarray(array)
        image.save(output)


if __name__=="__main__":
    main()

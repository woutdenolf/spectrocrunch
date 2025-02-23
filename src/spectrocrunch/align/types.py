from ..utils.Enum import Enum

dataType = Enum(["h5", "h5ext", "singlefile", "nparray"])
alignType = Enum(["full", "usetransfo", "calctransfo"])
# B-spline, moving least-squares
transformationType = Enum(
    ["translation", "rigid", "similarity", "affine", "projective"]
)

# Affine transformation:
#   Coordinate transformation: X' = R.X + T
#   Change of reference frame: X = R^-1.X' - R^-1.T
#
# Homogeneous coordinates: M = [[R,T],[P,1]]
#   [[Y'],[1]] = M . [[X],[1]]
#   X' = Y'/Z
#
# rotation: R = [[cos,-sin],[sin,cos]]
# reflection: R = [[-1,0],[0,-1]]
# scaling: R = [[s,0],[0,s]]
# aspect: R = [[a,0],[1/a,0]]
# shear: R = [[1,s],[s,1]]
#
# rigid(eucledian):                   reflection + rotation + translation   dof = 3      R = [[a,sqrt(1-a^2),tx],[sqrt(1-a^2),a,ty],[0,0,1]]
# similarity:                         rigid + scaling                       dof = 4      R = [[a,-b,tx],[b,a,ty],[0,0,1]]
# affine (parallel projection):       similarity + aspect + shear           dof = 6      R = [[a,b,tx],[c,d,ty],[0,0,1]]
# homography(perspective projection): affine + ...                          dof = 8      R = [[a,b,tx],[c,d,ty],[px,py,1]]
#
# http://www.robots.ox.ac.uk/~az/tutorials/cvpr03_part1.pdf
# http://morpheo.inrialpes.fr/people/Boyer/Teaching/M2R/geoProj.pdf
# https://ags.cs.uni-kl.de/fileadmin/inf_ags/3dcv-ws11-12/3DCV_WS11-12_lec04.pdf

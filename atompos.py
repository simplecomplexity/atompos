#! /usr/bin/env python
# Silicon crystal lattice position calculator

import math

## functions ##############################################################
## measure the distance between 2 points ########
def dist(x1, y1, z1, x2, y2, z2):
  d = math.sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2))
  return d

## main program ###########################################################
# number of cells for each axis
cell_x = 2;
cell_y = 2;
cell_z = 2;

# number of atoms (8 atoms for unit Silicon cell)
atoms = 8 * cell_x * cell_y * cell_z

# calculation of cell size
unit_cell = 5.43
cell_x_length = unit_cell * cell_x
cell_y_length = unit_cell * cell_y
cell_z_length = unit_cell * cell_z

# basic unit cell co-ordinate : 8 silicons
sub_b = [ [ 0.000,  0.000,  0.000], \
          [ 0.500,  0.500,  0.000], \
          [ 0.500,  0.000,  0.500], \
          [ 0.000,  0.500,  0.500], \
          [ 0.250,  0.250,  0.250], \
          [ 0.750,  0.750,  0.250], \
          [ 0.750,  0.250,  0.750], \
          [ 0.250,  0.750,  0.750] ]
# tetrahedral interstitial site vector
tet_b =   [ 0.500,  0.500,  0.500 ]
# hexagonal interstitial site vector
hex_b = [ [ 0.125,  0.125,  0.375], \
          [-0.125, -0.125,  0.375], \
          [-0.375, -0.125,  0.125], \
          [-0.125, -0.375,  0.125] ]
# bond-centerd site vector
bc_b  = [ [ 0.125,  0.125,  0.125], \
          [ 0.125, -0.125, -0.125], \
          [-0.125,  0.125, -0.125], \
          [-0.125, -0.125,  0.125] ]
# anti-bond interstitial site vector
ab_b  = [ [ 0.375,  0.375,  0.375], \
          [ 0.375, -0.375, -0.375], \
          [-0.375,  0.375, -0.375], \
          [-0.375, -0.375,  0.375] ]

## substitutional site data ######################
sub = []
for k in xrange(cell_z):
  for j in xrange(cell_y):
    for i in xrange(cell_x):
      for s in xrange(8):
        sub.append([(sub_b[s][0]+i)/cell_x, \
                    (sub_b[s][1]+j)/cell_y, \
                    (sub_b[s][2]+k)/cell_z])

print 'substitutional site = ' + str(len(sub)) + ' sites'
for i in xrange(len(sub)):
  print "[%.14f,%.14f,%.14f]" % (sub[i][0], sub[i][1], sub[i][2])

## tetrahedral site data #########################
tet = []
for k in xrange(cell_z):
  for j in xrange(cell_y):
    for i in xrange(cell_x):
      for s in xrange(4):
        x = (sub_b[s][0]+tet_b[0]+i)/cell_x
        y = (sub_b[s][1]+tet_b[1]+j)/cell_y
        z = (sub_b[s][2]+tet_b[2]+k)/cell_z
        if x >= 0 and x <= 1 and y >= 0 and y <= 1 and z >= 0 and z <= 1 :
          tet.append([x,y,z])

print 'tetrahedral site = ' + str(len(tet)) + ' sites'
for i in xrange(len(tet)):
  d = dist(0.5, 0.5, 0.5, tet[i][0], tet[i][1], tet[i][2])
  print "[%.14f,%.14f,%.14f] %.4f" % (tet[i][0], tet[i][1], tet[i][2], d)

## hexagnal site data ############################
hex = []
for k in xrange(cell_z):
  for j in xrange(cell_y):
    for i in xrange(cell_x):
      for s in xrange(4):
        for h in xrange(len(hex_b)):
          x = (sub_b[s+4][0]+hex_b[h][0]+i)/cell_x
          y = (sub_b[s+4][1]+hex_b[h][1]+j)/cell_y
          z = (sub_b[s+4][2]+hex_b[h][2]+k)/cell_z
          if x >= 0 and x <= 1 and y >= 0 and y <= 1 and z >= 0 and z <= 1 :
            hex.append([x,y,z])

print 'hexagonal site = ' + str(len(hex)) + ' sites'
for i in xrange(len(hex)):
  d = dist(0.5, 0.5, 0.5, hex[i][0], hex[i][1], hex[i][2])
  print "[%.14f,%.14f,%.14f] %.4f" % (hex[i][0], hex[i][1], hex[i][2], d)

## bond-centerd site data ########################
bc  = []
for k in xrange(cell_z):
  for j in xrange(cell_y):
    for i in xrange(cell_x):
      for s in xrange(4):
        for b in xrange(len(bc_b)):
          x = (sub_b[s][0]+bc_b[b][0]+i)/cell_x
          y = (sub_b[s][1]+bc_b[b][1]+j)/cell_y
          z = (sub_b[s][2]+bc_b[b][2]+k)/cell_z
          if x >= 0 and x <= 1 and y >= 0 and y <= 1 and z >= 0 and z <= 1 :
            bc.append([x,y,z])

print 'bond-centered site = ' + str(len(bc)) + ' sites'
for i in xrange(len(bc)):
  d = dist(0.5, 0.5, 0.5, bc[i][0], bc[i][1], bc[i][2])
  print "[%.14f,%.14f,%.14f] %.4f" % (bc[i][0], bc[i][1], bc[i][2], d)

## anti-bond interstitial site data ##############
ab  = []
for k in xrange(cell_z):
  for j in xrange(cell_y):
    for i in xrange(cell_x):
      for s in xrange(4):
        for a in xrange(len(ab_b)):
          x = (sub_b[s][0]+ab_b[a][0]+i)/cell_x
          y = (sub_b[s][1]+ab_b[a][1]+j)/cell_y
          z = (sub_b[s][2]+ab_b[a][2]+k)/cell_z
          if x >= 0 and x <= 1 and y >= 0 and y <= 1 and z >= 0 and z <= 1 :
            ab.append([x,y,z])

print 'anti-bond intestitial site = ' + str(len(ab)) + ' sites'
for i in xrange(len(ab)):
  d = dist(0.5, 0.5, 0.5, ab[i][0], ab[i][1], ab[i][2])
  print "[%.14f,%.14f,%.14f] %.4f" % (ab[i][0], ab[i][1], ab[i][2], d)


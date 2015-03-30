from __future__ import print_function, absolute_import, division 
from numpy import *


def GenerateSquare(nx, dx, ny, dy):
    nodelist = {}
    id_counter = 0
    for j in range(0, ny):
        for i in range(0, nx):
            x = i * dx
            y = j * dy
            nodelist.update({id_counter: array([x, y])})
            id_counter += 1

    # generate triangle
    connectivities = {}
    el_counter = 0
    for j in range(0, ny - 1):
        base_i = j * nx
        base_i1 = (j + 1) * nx
        for i in range(0, nx - 1):
            connectivities.update(
                {el_counter: [0, [base_i + i, base_i1 + 1 + i, base_i1 + i]]})
            el_counter += 1
            connectivities.update(
                {el_counter: [0, [base_i + i, base_i + 1 + i, base_i1 + 1 + i]]})
            el_counter += 1
            
    #generate face connectivities
    face_connectivities = {}
    cond_counter = 0
    
    #create faces of the bottom side
    for i in range(0, nx - 1):
        face_connectivities.update( {cond_counter: [0, [i, i+1] ] } )
        cond_counter += 1
        
    #create faces of the right side
    for i in range(nx-1, ny*nx-1, nx ):
        face_connectivities.update( {cond_counter: [0, [i, i+nx] ]} )
        cond_counter += 1

    #create faces of the top side
    for i in range(ny*nx-1, (ny-1)*nx, -1 ):
        face_connectivities.update( {cond_counter: [0, [i, i-1] ]} )    
        cond_counter += 1

    #create faces of the left side
    for i in range((ny-1)*nx, 0, -nx ):
        face_connectivities.update( {cond_counter: [0, [i, i-nx] ]} )
        cond_counter += 1
        
        
    return nodelist, connectivities, face_connectivities



#nodelist,connectivities,face_connectivities = GenerateSquare(2,0.1,3,0.1)

# for aaa in connectivities.values():
    # print aaa

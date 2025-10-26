
"""
    Some experiments with calling Fortran from Python

    You have to build the shared library `libfat.so` first using:
        FoBiS.py build -f fobis -mode shared-intel

    author: Jacob Williams, 5/29/2017
"""

from ctypes import *

# load the shared library:
FAT = CDLL('../../lib/libfat.so')

# The routines we want are in the `c_interface_module`:
#
#  * initialize_geopotential_model(itype,gravfile,n,m)
#  * destroy_geopotential_model(cp)
#  * get_acceleration(cp,n,m,rvec,acc)

# define the procedure interfaces:
initialize_geopotential_model = FAT.initialize_geopotential_model
destroy_geopotential_model    = FAT.destroy_geopotential_model
get_acceleration              = FAT.get_acceleration
initialize_geopotential_model.argtypes = [c_int,c_char_p,c_int,c_int]
initialize_geopotential_model.restype  = POINTER(c_int)
destroy_geopotential_model.argtypes    = [POINTER(c_int)]
destroy_geopotential_model.restype     = None
get_acceleration.argtypes              = [POINTER(c_int),c_int,c_int,POINTER(c_double),POINTER(c_double)]
get_acceleration.restype               = None

VecType = c_double*3  # to define 3x1 vectors
gravfile = c_char_p(b'../../grav/GGM03C.GEO')  # gravity coefficient file name
itype = c_int(1)  # model type (1 is only one supported currently)
n = c_int(8)  # degree
m = c_int(8)  # order

#----------------------------
print('')
print('... string test...')
print('')

return_a_string = FAT.return_a_string
return_a_string.restype = None
return_a_string.argtypes = [c_int,POINTER(c_char_p)]

buff_len = 9 # length of the buffer to use
buffer = c_char_p(b' '.ljust(buff_len))  # allocate the buffer
for i in range(13):
    return_a_string(c_int(i),byref(buffer))
    print('->'+buffer.value+'<-')
print('')
#----------------------------

print('calling initialize_geopotential_model...')

cp = initialize_geopotential_model(itype,gravfile,n,m)

print('calling get_acceleration...')

acc = VecType(0.0,0.0,0.0)  # acceleration vector

for r in range(10000, 20000, 1000):
    rvec = VecType(r,r,r)  # position vector
    get_acceleration(cp,n,m,rvec,acc)
    print(acc[0],acc[1],acc[2])

print('calling destroy_geopotential_model...')

destroy_geopotential_model(cp)


from astropy.io import ascii

data = ascii.read("vikaeus.dat", delimiter=' ')  
print(data)
print(data.columns)


print(data['Z'].data)
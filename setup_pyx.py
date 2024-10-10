import pyximport
pyximport.install()
from c_msdl import c_msdl
print (c_msdl(1000000))
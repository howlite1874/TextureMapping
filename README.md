# TextureMapping
To compile:

module add legacy-eng
qmake -project QT+=opengl
qmake
make

Function:convert model to a properly textured surface and to compute a normal map, with 
the xyz of the normals mapped to rgb

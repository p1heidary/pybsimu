
# cdef extern from "src/vec3d.hpp":
# 	cdef cppclass Vec3D:
# 		double p[3];
# 		Vec3D( double x, double y, double z ) { p[0] = x; p[1] = y; p[2] = z; }
	
# cdef extern from "src/solid.hpp":
# 	cdef cppclass Solid:
# 		Solid() except +
# 		bool inside( Vec3D &x )
from libc.stdint cimport uint32_t, int64_t

cdef extern from "<iostream>" namespace "std":
    cdef cppclass string:
    	string(const char*) except +
    # cdef cppclass ostream:
    #     ostream& write(const char*, int) except +

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+

cdef extern from "src/types":
	cdef enum field_extrpl_e

cdef extern from "src/meshvectorfield.hpp":
	cdef cppclass MeshVectorField:
		MeshVectorField() except +

cdef extern from "src/meshscalarfield.hpp":
	cdef cppclass MeshScalarField:
		MeshScalarField() except +

cdef extern from "src/ibsimu.hpp":
	cdef cppclass ibsimu:
		ibsimu() except +
		~ ibsimu()
		void set_message_threshold(int)
		void set_thread_count(int)
		# std::ostream message(int) fixme

cdef extern from "src/error.hpp":
	cdef cppclass Error:
		Error() except +
		void print_error_message( std::ostream, bool );

cdef extern from "src/geometry.hpp":
	cdef cppclass Geometry:
		Geometry( geom_mode_e, Int3D, Vec3D, double ) except + # fixus custom type
		void set_solid( int, const Solid * )
		void set_boundary( int, const Bound )
		void build_mesh()
		Vec3D origo()
		Vec3D max()

cdef extern from "src/epot_field.hpp":
	cdef cppclass EpotField:
		EpotField( const Geometry ) except + # fixme const

cdef extern from "src/epot_efield.hpp":
	cdef cppclass EpotEfield:
		EpotEfield( const EpotField ) except + # fixme const
		void set_extrapolation( field_extrpl_e extrpl[6] ) # fixme array passed
		void recalculate()

cdef extern from "src/epot_bicgstabsolver.hpp":
	cdef cppclass EpotBiCGSTABSolver:
		EpotBiCGSTABSolver( Geometry &geom, double eps = 1.0e-4, 
			uint32_t imax = 10000, double newton_eps = 1.0e-4, 
			uint32_t newton_imax = 10, bool gnewton = true ) except +
		void solve( MeshScalarField &epot, const ScalarField &scharge )
		uint32_t get_iter() const

cdef extern from "mydxffile.hpp":
	cdef cppclass MyDXFFile:
		MyDXFFile() except +
		~MyDXFFile()
		void set_warning_level( int )
		void read( const std::string &filename ); # fixme std &

cdef extern from "dxf_solid.hpp":
	cdef cppclass DXFSolid:
		DXFSolid( MyDXFFile *dxffile, const std::string &layername ) except +
		void scale( double )

# cdef extern from "src/gtkplotter.hpp":
# cdef extern from "src/trajectorydiagnostics.hpp":
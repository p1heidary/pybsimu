from libc.stdint cimport uint32_t, int64_t, int32_t
# from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "<iostream>" namespace "std":
	cdef cppclass string:
		string(const char*) except +
	cdef cppclass ostream:
		ostream& write(const char*, int) except +

cdef extern from "<fstream>" namespace "std":
	cdef cppclass ofstream:
		ofstream(const char*) except +

cdef extern from "src/types.hpp":
	cdef enum field_extrpl_e:
		pass
	cdef enum geom_mode_e:
		pass

cdef extern from "src/error.hpp":
	cdef cppclass Error:
		Error() except +
		void print_error_message( ostream, bool )

cdef extern from "src/vec3d.hpp":

	cdef cppclass Int3D:
		Int3D() except +
		Int3D( int32_t , int32_t , int32_t ) except +

	cdef cppclass Vec3D:
		Vec3D() except +
		Vec3D( double , double , double ) except +

cdef extern from "src/meshvectorfield.hpp":
	cdef cppclass MeshVectorField:
		MeshVectorField() except +

cdef extern from "src/meshscalarfield.hpp":
	cdef cppclass MeshScalarField:
		MeshScalarField() except +

cdef extern from "src/scalarfield.hpp":
	cdef cppclass ScalarField:
		pass

cdef extern from "src/ibsimu.hpp":
	cdef cppclass ibsimu:
		ibsimu() except +
		# void set_message_threshold(int)
		# void set_thread_count(int)
		ostream message(int)

cdef extern from "src/solid.hpp":
	cdef cppclass Solid:
		Solid() except +

cdef extern from "src/geometry.hpp":
	
	cdef cppclass Bound:
		pass

	cdef cppclass Geometry:
		Geometry( geom_mode_e, Int3D, Vec3D, double ) except + 
		void set_solid( int, const Solid * )
		void set_boundary( int, const Bound )
		void build_mesh()
		Vec3D origo()
		Vec3D max()

cdef extern from "src/epot_field.hpp":
	cdef cppclass EpotField:
		EpotField( const Geometry ) except + 

cdef extern from "src/epot_efield.hpp":
	cdef cppclass EpotEfield:
		EpotEfield( const EpotField ) except + 
		# void set_extrapolation( field_extrpl_e* extrpl ) # fixme array passed
		void recalculate()

cdef extern from "src/particledatabase.hpp":
	cdef cppclass ParticleDataBase3D:
		ParticleDataBase3D( const Geometry& ) except +

cdef extern from "src/trajectorydiagnostics.hpp":
	cdef cppclass TrajectoryDiagnosticData:
		TrajectoryDiagnosticData() except +

cdef extern from "src/epot_bicgstabsolver.hpp":
	cdef cppclass EpotBiCGSTABSolver:
		EpotBiCGSTABSolver( Geometry&, double eps = 1.0e-4, 
			uint32_t imax = 10000, double newton_eps = 1.0e-4, 
			uint32_t newton_imax = 10, bool gnewton = True ) except +
		void solve( MeshScalarField&, const ScalarField& ) # fixme refrence
		uint32_t get_iter() const

cdef extern from "mydxffile.hpp":
	cdef cppclass MyDXFFile:
		MyDXFFile() except +
		# void set_warning_level( int )
		void read( const string& ) # fixme std &
		void set_warning_level ( int )

cdef extern from "dxf_solid.hpp":
	cdef cppclass DXFSolid:
		DXFSolid( MyDXFFile*, const string& ) except +
		void scale( double )
		void define_2x3_mapping()

# cdef extern from "src/gtkplotter.hpp":
# cdef extern from "src/trajectorydiagnostics.hpp":
#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "multimeshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "readascii.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;


VectorField *build_bfield( void )
{
    bool fout_all[3] = {true, true, true};

    bool sol1_fout[3] = {true, true, false};
    MeshVectorField sol1( MODE_CYL, sol1_fout, 1.0e-3, 1.0, "solfield_100A.txt" );
    field_extrpl_e sol1_extrpl[6] = { FIELD_ZERO, FIELD_ZERO, FIELD_ZERO, FIELD_ZERO, FIELD_ZERO, FIELD_ZERO };
    sol1.set_extrapolation( sol1_extrpl );

    Int3D size( 101, 101, 1001 );
    Vec3D origo( -50e-3, -50e-3, -500e-3 );
    double h = 1.0e-3;
    MeshVectorField *sol1_xyz = new MeshVectorField( MODE_3D, fout_all, size, origo, h, sol1 );

    return( sol1_xyz );
}


void simu( int argc, char **argv )
{
    double h = 2e-3;
    double sizereq[3] = { 100.0e-3, 100.0e-3, 1000.0e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( -50.0e-3, -50e-3, -500.0e-3 );
    Geometry geom( MODE_3D, meshsize, origo, h );

    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    VectorField *bfield = build_bfield();

    EpotEfield efield( epot );
    ParticleDataBase3D pdb( geom );

    uint32_t N = 10000;
    double I = 150.0e-6;
    double q = 8.0;
    double m = 40.0;
    double E0 = 10e3*q;
    double e = 100e-6;
    double xmax = 10e-3;
    double xpmax = 10e-3;
    double g = xpmax*xpmax/e;
    double b = xmax*xmax/e;
    double a = -sqrt(b*g-1.0);
    Vec3D c( 0, 0, -500e-3 );
    Vec3D dir1( 1, 0, 0 );
    Vec3D dir2( 0, 1, 0 );
    pdb.add_3d_KV_beam_with_emittance( N, I, q, m, E0, a, b, e, a, b, e, c, dir1, dir2 );
    pdb.iterate_trajectories( scharge, efield, *bfield );

    if( true ) {
	MeshScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
	
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_bfield( bfield );
	plotter.set_scharge( &scharge );
	plotter.set_trajdens( &tdens );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }

    // Write (x,x',y,y') at input and output
    ofstream of( "data.txt" );
    for( uint32_t a = 0; a < pdb.size(); a++ ) {
	const ParticleP3D &p0 = pdb.trajectory_point(a,0);
	const ParticleP3D &p1 = pdb.particle(a).x();
	of << p0(1) << " "
	   << p0(2)/p0(6) << " "
	   << p0(3) << " "
	   << p0(4)/p0(6) << " "
	   << p1(1) << " "
	   << p1(2)/p1(6) << " "
	   << p1(3) << " "
	   << p1(4)/p1(6) << "\n";
    }
}


int main( int argc, char **argv )
{
    try {
	//ibsimu.set_message_output( "ibsimu" + stamp + ".txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}

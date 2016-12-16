#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"


using namespace std;


double Vpuller = -50e3;
double Vgnd = -40e3;

double Te = 5.0;
double Up = 5.0;
double E0 = 0.5*Te;

double q = 1.0;
double m = 1.0;
double Tt = 0.5;
double r0 = 8e-3;
double J = 300.0;

double sc_alpha = 0.9;
double h = 0.25e-3;
double nperh = 100;


void simu( int *argc, char ***argv )
{
    Vec3D origo( -2e-3, 
		 0, 
		 0 );
    double sizereq[3] = { 200e-3,
                           50e-3, 
                            0 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
		    1 );
    Geometry geom( MODE_2D, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 2 );
    dxffile->read( "slit.dxf" );
    DXFSolid *s1 = new DXFSolid( dxffile, "plasma" );
    s1->scale( 1e-3 );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "puller" );
    s2->scale( 1e-3 );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "gnd" );
    s3->scale( 1e-3 );
    geom.set_solid( 9, s3 );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) );

    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Vpuller) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, Vgnd) );

    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma init_plasma( AXIS_X, 2e-4 );
    solver.set_initial_plasma( Up, &init_plasma );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    MeshVectorField bfield;

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );

    for( size_t a = 0; a < 20; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	if( a == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    solver.set_pexp_plasma( rhoe, Te, Up );
	}

	solver.solve( epot, scharge_ave );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

        pdb.clear(); 
	int Npart = floor(r0/h*nperh+0.5);
	pdb.add_2d_beam_with_energy( Npart, J, q, m,
				     E0, 0, Tt,
				     geom.origo(0), 0,
				     geom.origo(0), r0 );
	
        pdb.iterate_trajectories( scharge, efield, bfield );

	// 90 % space charge compensation at x>100 mm
	/*
	for( uint32_t i = 0; i < geom.size(0); i++ ) {
	    double x = geom.origo(0) + geom.h()*x;
	    if( x > 100e-3 ) {
		for( uint32_t j = 0; j < geom.size(1); j++ )
		    scharge(i,j) *= 0.1;
	    }
	}
	*/

	if( a == 0 ) {
	    scharge_ave = scharge;
	} else {
	    double sc_beta = 1.0-sc_alpha;
	    uint32_t nodecount = scharge.nodecount();
	    for( uint32_t b = 0; b < nodecount; b++ ) {
		scharge_ave(b) = sc_alpha*scharge(b) + sc_beta*scharge_ave(b);
	    }
	}

	TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_Y );
        diagnostics.push_back( DIAG_YP );
        pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0)-geom.h(), diagnostics );
        Emittance emit( tdata(0).data(), tdata(1).data() );

        // Output
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit.alpha() << " "
             << emit.beta() << " "
             << emit.epsilon() << "\n";
        dout.close();

    }
    
    if( false ) {
	MeshScalarField tdens(geom);
	pdb.build_trajectory_density_field(tdens);
	
	GTKPlotter plotter( argc, argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_efield( &efield );
	plotter.set_bfield( &bfield );
	plotter.set_trajdens( &tdens );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }

    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_Y );
    diagnostics.push_back( DIAG_YP );
    pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0)-geom.h(), diagnostics );
    Emittance emit( tdata(0).data(), tdata(1).data() );

    ofstream dout( "diag.txt", ios_base::app );
    dout << h << " " << emit.epsilon() << "\n";
}


int main( int argc, char **argv )
{
    h = atof( argv[1] );

    remove( "emit.txt" );

    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( &argc, &argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}

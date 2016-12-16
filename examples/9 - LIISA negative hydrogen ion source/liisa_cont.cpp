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
#include "random.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "transformation.hpp"
#include "readascii.hpp"



using namespace std;


double sc_alpha = 1.0;
double Veinzel2 = 15.00e3;


void simu( int argc, char **argv )
{
    double h = 1e-3;
    Vec3D origo( -30e-3, -30e-3, 150e-3 );
    double sizereq[3] = {  60e-3,
                           60e-3, 
                          140e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
		    (int)floor(sizereq[2]/h)+1 );
    Geometry geom( MODE_3D, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 1 );
    dxffile->read( "LIISA_geom.dxf" );
    MyDXFEntities *e = dxffile->get_entities();
    MyDXFEntitySelection *sel = e->selection_all();
    e->scale( sel, dxffile, 1.0e-3 );

    DXFSolid *s1 = new DXFSolid( dxffile, "gnd" );
    s1->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "einzel_2" );
    s2->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 8, s2 );

    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.00) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   0.00) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.00) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   Veinzel2) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    ReadAscii ra( "part.txt" );
    for( size_t a = 0; a < 5; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	solver.solve( epot, scharge_ave );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

        pdb.clear(); 
	for( uint32_t b = 0; b < ra.rows(); b++ ) {
	    pdb.add_particle( ra[0][b], -1.0, 1.0, ParticleP3D( 0.0,
							       ra[1][b], ra[2][b],
							       ra[3][b], ra[4][b],
							       ra[5][b], ra[6][b] ) );
	}
        pdb.iterate_trajectories( scharge, efield, bfield );

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
        diagnostics.push_back( DIAG_X );
        diagnostics.push_back( DIAG_XP );
        diagnostics.push_back( DIAG_Y );
        diagnostics.push_back( DIAG_YP );
        pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diagnostics );
        Emittance emit_xxp( tdata(0).data(), tdata(1).data() );
        Emittance emit_yyp( tdata(2).data(), tdata(3).data() );

        // Output
        ofstream dout( "emit_cont.txt", ios_base::app );
        dout << emit_xxp.alpha() << " "
             << emit_xxp.beta() << " "
             << emit_xxp.epsilon() << " "
	     << emit_yyp.alpha() << " "
             << emit_yyp.beta() << " "
             << emit_yyp.epsilon() << "\n";
        dout.close();
    }

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );
    //tdens.debug_print( cout );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_bfield( &bfield );
    plotter.set_efield( &efield );
    plotter.set_scharge( &scharge );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 0 );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_font_size( 16 );
    geomplotter.set_trajdens( &tdens );
    geomplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    geomplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );

    geomplotter.set_view( VIEW_ZX, -1 );
    geomplotter.plot_png( "geom_cont_zx.png" );
    geomplotter.set_view( VIEW_ZY, -1 );
    geomplotter.plot_png( "geom_cont_zy.png" );
}


int main( int argc, char **argv )
{
    remove( "emit.txt" );

    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( cout );
	exit( 1 );
    }

    return( 0 );
}

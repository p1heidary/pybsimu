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
#include "convergence.hpp"



using namespace std;


double Vsource = 5.90e3;
double Vpuller = 0.90e3;
double Vext_einzel = 3.00e3;
double Vaccel = 3.04e3;
double Veinzel1 = 12.00e3;

double sc_alpha = 0.5;
double J = 20.0;
double eir = 12;
double Tt = 0.5;
double E0 = 2.0;
double r0 = 8.0e-3;
double nperh = 100;

double Tp = 0.2;
double Rf = 0.25;
double Up = 5.0;

double Usup = 10.0;


class ForcedPot : public CallbackFunctorB_V {

public:

    ForcedPot() {}
    ~ForcedPot() {}

    virtual bool operator()( const Vec3D &x ) const {
        return( x[2] < 0.1e-3 && x[0]*x[0]+x[1]*x[1] > 7e-3*7e-3 );
    }
};


void simu( int argc, char **argv )
{
    double h = 0.5e-3;
    Vec3D origo( -30e-3, -30e-3, -4e-3 );
    double sizereq[3] = {  60e-3,
                           60e-3, 
                          154e-3 };
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

    DXFSolid *s1 = new DXFSolid( dxffile, "plasma" );
    s1->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "puller" );
    s2->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "ext_einzel" );
    s3->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "accel" );
    s4->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 10, s4 );
    DXFSolid *s5 = new DXFSolid( dxffile, "gnd" );
    s5->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 11, s5 );
    DXFSolid *s6 = new DXFSolid( dxffile, "einzel_1" );
    s6->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 12, s6 );

    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.00) );
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.00) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   Vsource) ); // Dirichlet for cont

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.00) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   Vpuller) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   Vext_einzel) );
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,   Vsource+Vaccel) );
    geom.set_boundary( 11,  Bound(BOUND_DIRICHLET,   Vsource) );
    geom.set_boundary( 12,  Bound(BOUND_DIRICHLET,   Vsource+Veinzel1) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_Z, 0.1e-3 );
    solver.set_nsimp_initial_plasma( &initp );
    ForcedPot force;
    solver.set_forced_potential_volume( 0.0, &force );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    // Define magnetic field
    bool fout[3] = {true, true, true};
    MeshVectorField bfield( MODE_3D, fout, 1.0e-3, 1.0, "LIISA_bfield.txt" );
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
				     FIELD_ZERO, FIELD_ZERO, 
				     FIELD_ZERO, FIELD_ZERO };
    bfield.set_extrapolation( bfldextrpl );
    //bfield.debug_print();

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

    NPlasmaBfieldSuppression psup( epot, Usup );
    pdb.set_bfield_suppression( &psup );

    double rho_tot = 0.0;
    for( size_t a = 0; a < 5; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	if( a == 1 ) {
	    std::vector<double> Ei, rhoi;
	    Ei.push_back( Tp );
	    rhoi.push_back( (1-Rf)*rho_tot );
	    double rhop = Rf*rho_tot;
            solver.set_nsimp_plasma( rhop, Up, rhoi, Ei );
        }

	solver.solve( epot, scharge_ave );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

        pdb.clear(); 
	int Npart = floor((M_PI*r0*r0)/(h*h)*nperh+0.5);
	pdb.add_cylindrical_beam_with_energy( Npart, -J, -1.0, 1.0, 
					      E0, 0.0, Tt, 
					      Vec3D(0,0,geom.origo(2)), // center
					      Vec3D(1,0,0), // dir1
					      Vec3D(0,1,0), // dir2
					      r0 );
	pdb.add_cylindrical_beam_with_energy( Npart, -J*eir, -1.0, 1.0/1836.15, 
					      E0, 0.0, Tt, 
					      Vec3D(0,0,geom.origo(2)), // center
					      Vec3D(1,0,0), // dir1
					      Vec3D(0,1,0), // dir2
					      r0 );
        pdb.iterate_trajectories( scharge, efield, bfield );
	rho_tot = pdb.get_rhosum();

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
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit_xxp.alpha() << " "
             << emit_xxp.beta() << " "
             << emit_xxp.epsilon() << " "
	     << emit_yyp.alpha() << " "
             << emit_yyp.beta() << " "
             << emit_yyp.epsilon() << "\n";
        dout.close();
    }

    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

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

    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_CURR  );
    diagnostics.push_back( DIAG_X  );
    diagnostics.push_back( DIAG_VX );
    diagnostics.push_back( DIAG_Y  );
    diagnostics.push_back( DIAG_VY );
    diagnostics.push_back( DIAG_Z  );
    diagnostics.push_back( DIAG_VZ );
    TrajectoryDiagnosticData tdata;
    pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2), diagnostics );

    ofstream pout( "part.txt" );
    for( uint32_t a = 0; a < tdata.traj_size(); a++ ) {
	for( uint32_t b = 0; b < 7; b++ )
	    pout << tdata(a,b) << " ";
	pout << "\n";
    }
    pout.close();

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 0 );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_font_size( 16 );
    geomplotter.set_trajdens( &tdens );
    geomplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    geomplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );
    std::vector<double> eqlines;
    eqlines.push_back( 1.0 );
    eqlines.push_back( 2.0 );
    eqlines.push_back( 4.0 );
    geomplotter.set_eqlines_manual( eqlines );

    geomplotter.set_view( VIEW_ZX, -1 );
    geomplotter.plot_png( "geom_zx.png" );
    geomplotter.set_view( VIEW_ZY, -1 );
    geomplotter.plot_png( "geom_zy.png" );
    geomplotter.set_view_si( VIEW_XY, 5.0e-3 );
    geomplotter.plot_png( "geom_xy.png" );
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

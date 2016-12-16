#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"


using namespace std;


double sc_alpha = 0.5;
double r0 = 8e-3;

double h = 2e-4;
double nperh = 200.0;
double Npart = nperh*r0/h;

// Beam parameters
double J = 35.0;
double eir = 50;
double Tt = 1.0;
double E0 = 2.0;

// Plasma parameters
double Tp = 0.2;
double Rf = 0.5;
double Up = 5.0;

string stamp;



class ForcedPot : public CallbackFunctorB_V {

public:

    ForcedPot() {}
    ~ForcedPot() {}

    virtual bool operator()( const Vec3D &x ) const {
        return( x[0] < 0.2e-3 && x[1] > 6e-3 );
    }
};


class THCallback : public TrajectoryHandlerCallback {
public:

    THCallback() {}

    virtual ~THCallback() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticlePCyl *pcur = (ParticlePCyl *)( xcur );
        ParticlePCyl *pend = (ParticlePCyl *)( xend );
        // Kill particles with mass less than 0.5 atomic mass and x-coordinate more than 40 mm.
        if( particle->m() < 0.5*MASS_U && (*pcur)[PARTICLE_X] >= 40e-3 ) {
            *pend = *pcur;
            particle->set_status( PARTICLE_COLL );
        }
    }
};


void simu( int argc, char **argv )
{
    double start = -2.0e-3;
    double sizereq[3] = { 82.0e-3,
                           30.0e-3, 
                            0.0e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( start, 0.0, 0.0 );
    Geometry geom( MODE_CYL, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile( "geom.dxf" );
    dxffile->set_warning_level( 2 );
    MyDXFEntities *e = dxffile->get_entities();
    MyDXFEntitySelection *sel = e->selection_all();
    e->scale( sel, dxffile, 1.0e-3 );

    DXFSolid *s1 = new DXFSolid( dxffile, "plasma" );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "puller" );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "einzel" );
    geom.set_solid( 9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "gnd" );
    geom.set_solid( 10, s4 );
    
    geom.set_boundary(  1,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  2,  Bound(BOUND_DIRICHLET,  20.0e3) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   8.0e3) ); // Puller
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   6.0e3) ); // Einzel/Dump
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  20.0e3) ); // GND
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 0.0 );
    solver.set_nsimp_initial_plasma( &initp );
    ForcedPot force;
    solver.set_forced_potential_volume( 0.0, &force );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    MeshVectorField bfield;

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( false );

    THCallback thc;
    pdb.set_trajectory_handler_callback( &thc );

    double rho_tot;
    for( size_t a = 0; a < 25; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	if( a == 1 ) {
	    std::vector<double> Ei, rhoi;
	    Ei.push_back( Tp );
	    rhoi.push_back( (1.0-Rf)*rho_tot );
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
	pdb.add_2d_beam_with_energy( Npart, -J, -1.0, 1.0, 
				     E0, 0.0, Tt, 
				     start, 0,
				     start, r0 );
	pdb.add_2d_beam_with_energy( Npart, -J*eir, -1.0, 1.0/1836.15, 
				     E0, 0.0, Tt, 
				     start, 0,
				     start, r0 );
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
        diagnostics.push_back( DIAG_R );
        diagnostics.push_back( DIAG_RP );
        pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0), diagnostics );
        Emittance emit_rrp( tdata(0).data(), tdata(1).data() );

        // Output
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit_rrp.alpha() << " "
             << emit_rrp.beta() << " "
             << emit_rrp.epsilon() << "\n";
        dout.close();
    }
    
    //geom.save( "geom.dat" );
    //epot.save( "epot.dat" );
    //pdb.save( "pdb.dat" );

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 10 );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_font_size( 20 );
    //MeshScalarField tdens( geom );
    //pdb.build_trajectory_density_field( tdens );
    //geomplotter.set_trajdens( &tdens );
    //geomplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    //geomplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );
    std::vector<double> eqlines;
    eqlines.push_back( 1.0 );
    eqlines.push_back( 2.0 );
    eqlines.push_back( 4.0 );
    geomplotter.set_eqlines_manual( eqlines );
    geomplotter.plot_png( "geom" + stamp + ".png" );
    geomplotter.set_ranges( geom.origo(0), 0, 0.01, 0.01 );
    geomplotter.plot_png( "zoom" + stamp + ".png" );

    ParticleDiagPlotter pdplotter( geom, pdb, AXIS_X, geom.max(0),
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Z, DIAG_ZP );
    pdplotter.plot_png( "emit" + stamp + ".png" );

    if( true ) {
	MeshScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
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
    }
}


int main( int argc, char **argv )
{
    try {
	//ibsimu.set_message_output( "ibsimu_out.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}

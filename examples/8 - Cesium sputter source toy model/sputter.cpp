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
#include "particlediagplotter.hpp"
#include "axisymmetricvectorfield.hpp"


using namespace std;

double Vcathode = 0.0;
double Vionizer = 2e3;
double Vextractor = 15e3;


double mn = 1.0;
double m = 133.0;
double Tt = 0.1;
double Ttn = 0.1;
double E0 = 3.0;
double E0n = 3.0;
double r0 = 25e-3;
double I = 100e-6;
double J = I/(r0*r0*M_PI);


double sc_alpha = 0.4;
uint32_t Npart = 1000;


// Note, mapping used after scale!
double ionizer_xoff = 35e-3;
double ionizer_r = 35e-3;


Vec3D ionizershape( const Vec3D &x )
{
    if( (x[0]-1e3*ionizer_xoff)*(x[0]-1e3*ionizer_xoff) + x[1]*x[1] < 1e3*ionizer_r*1e3*ionizer_r )
	return( Vec3D(numeric_limits<double>::infinity(),
		      numeric_limits<double>::infinity(),
		      numeric_limits<double>::infinity() ) );
    return( x );
}


class SECallback : public TrajectoryEndCallback {

    Geometry &_geom;
    MTRandom  _rand;
    
public:

    SECallback( Geometry &geom ) 
        : _geom(geom), _rand(2) {

	_rand.set_transformation( 0, Gaussian_Transformation() );
	_rand.set_transformation( 1, Gaussian_Transformation() );
    }

    virtual ~SECallback() {}

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) {

        ParticleCyl *pcyl = (ParticleCyl *)( particle );
        Vec3D loc = pcyl->location();
        Vec3D vel = pcyl->velocity();
        double E = 0.5*pcyl->m()*vel.ssqr()/CHARGE_E;

       // Make secondaries based on location and energy
        if( fabs(loc[0]-4e-3) < 1e-6 && loc[1] <= 4e-3 && E > 1e3 ) {

            // Adjust location off the surface
            loc[0] += 0.01*_geom.h();

	    double vt[2];
	    _rand.get( vt );
	    double dvt = sqrt(Ttn*CHARGE_E/(m*MASS_U));
            vel[0] = sqrt( 2.0*E0n*CHARGE_E/(mn*MASS_U) );
	    vel[1] = dvt*vt[0];
	    vel[2] = dvt*vt[1];

            ParticleDataBaseCyl *pdbcyl = (ParticleDataBaseCyl *)( pdb );
            pdbcyl->add_particle( 0.1*pcyl->IQ(), -1.0, mn, ParticlePCyl( 0.0, 
									  loc[0], vel[0], 
									  loc[1], vel[1], 
									  vel[2]/loc[1] ) );
        }
    }
};


void simu( int *argc, char ***argv )
{
    double h = 0.5e-3;
    Vec3D origo( 0e-3, 0, 0 );
    double sizereq[3] = { 180e-3,
                          60e-3, 
                          0 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
		    1);
    Geometry geom( MODE_CYL, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 2 );
    dxffile->read( "sputter.dxf" );
    DXFSolid *s1 = new DXFSolid( dxffile, "cathode" );
    s1->scale( 1e-3 );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "cathode shroud" );
    s2->scale( 1e-3 );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "ionizer" );
    s3->scale( 1e-3 );
    s3->define_2x3_mapping( ionizershape );
    geom.set_solid( 9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "ionizer shroud" );
    s4->scale( 1e-3 );
    geom.set_solid( 10, s4 );
    DXFSolid *s5 = new DXFSolid( dxffile, "extractor" );
    s5->scale( 1e-3 );
    geom.set_solid( 11, s5 );

    geom.set_boundary(  1, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary(  2, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary(  3, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary(  4, Bound(BOUND_NEUMANN, 0.0) );

    geom.set_boundary(  7, Bound(BOUND_DIRICHLET, Vcathode) );
    geom.set_boundary(  8, Bound(BOUND_DIRICHLET, Vcathode) );
    geom.set_boundary(  9, Bound(BOUND_DIRICHLET, Vionizer) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, Vionizer) );
    geom.set_boundary( 11, Bound(BOUND_DIRICHLET, Vextractor) );

    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );

    SECallback secb( geom );
    pdb.set_trajectory_end_callback( &secb );

    MTRandom rng( 2 );
    rng.set_transformation( 0, Gaussian_Transformation() );
    rng.set_transformation( 1, Gaussian_Transformation() );

    for( size_t a = 0; a < 40; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	solver.solve( epot, scharge_ave );
	if( a > 0 && solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

	// Angular coordinates for start and stop point
	double astart = asin(5e-3/ionizer_r);
	double astop = acos((53e-3-ionizer_xoff)/ionizer_r);
	ibsimu.message(1) << "astart = " << 180.0*astart/M_PI << "\n";
	ibsimu.message(1) << "astop = " << 180.0*astop/M_PI << "\n";

        pdb.clear(); 
	for( uint32_t b = 0; b < Npart; b++ ) {
	    double da = (astop-astart)/Npart;
	    double a = astart + (b+0.5)*da;
	    double R = ionizer_r*sin(a);
	    double x = ionizer_xoff+(ionizer_r-geom.h()/10.0)*cos(a);
	    double r = (ionizer_r-geom.h()/10.0)*sin(a);
	    double v = sqrt(2.0*E0*CHARGE_E/(m*MASS_U));

	    double vt[2];
	    rng.get( vt );
	    double dvt = sqrt(Tt*CHARGE_E/(m*MASS_U));

	    double vx = -v*cos(a) - dvt*vt[0]*sin(a);
	    double vr = -v*sin(a) + dvt*vt[0]*cos(a);
	    double vtheta = dvt*vt[1];
	    pdb.add_particle( J*da*ionizer_r*2.0*M_PI*R, 1.0, m, ParticlePCyl(0,x,vx,r,vr,vtheta/r) );
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
    }
    
    MeshScalarField tdens(geom);
    pdb.build_trajectory_density_field(tdens);

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_bfield( &bfield );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
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

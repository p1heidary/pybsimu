
cdef class geometry:
  
  cdef Geometry *geom
  cdef MyDXFFile *dxffile 

  cdef EpotBiCGSTABSolver solver
  cdef EpotField epot
  cdef MeshScalarField scharge
  cdef MeshVectorField bfield
  cdef EpotEfield efield
  cdef ParticleDataBase3D pdb

  cdef TrajectoryDiagnosticData tdata

  def __cinit__(self, mode, mesh, origin, h):
    
    # Int3D my_mesh(mesh[0], mesh[1], mesh[2])
    cdef Int3D my_mesh(mesh[0], mesh[1], mesh[2])
    cdef Vec3D my_origin(origin[0], origin[1], origin[2])
    # my_modes = {'MODE_3D': MODE_3D}
    self.geom  = new Geometry(my_modes(mode), my_mesh, my_origin, h)

  def __dealloc__(self):
    del self.geom
  
  def set_dxf_solid(self, dxf_file):
    
    dxffile.set_warning_level( 2 )
    dxffile.read( dxf_file )
    # my_mappings = {'rotx': DXFSolid::rotz}

  def set_solid(self, name, scale, position, mapping):
    cdef s1 = DXFSolid( dxffile, "gnd" );
    s1.scale( scales[n] );
    s1.define_2x3_mapping( my_mappings(mappings[n]) );
    self.geom.set_solid( position, s1 )

  def set_boundary(self, position, type, value):
    geom.set_boundary( position, Bound(type, value) )

  def prepare():

    geom.build_surface()
    geom.build_mesh()

    self.solver  = EpotBiCGSTABSolver( geom )
    self.epot    = EpotField( geom )
    self.scharge = MeshScalarField( geom )
    self.efield  = EpotEfield( epot )
    self.pdb     = ParticleDataBase3D( geom )

  def set_extrapolation(field_extrapolation):
    efield.set_extrapolation(field_extrapolation)

  def set_mirror(mirror):
    pdb.set_mirror(mirror)

  def run(iteration):
    for i in range(iteration):
      self.solver.solve( self.epot, self.scharge )
      if ( solver.get_iter() == 0 ):
        print "No iterations, breaking major cycle\n";
        raise
        self.efield.recalculate()
        self.pdb.clear()
        self.pdb.add_cylindrical_beam_with_energy( 1000, J, q, m,
                E0, 0, Tt,
                Vec3D(0,0,0),
                Vec3D(1,0,0),
                Vec3D(0,1,0),
                r0 );
  
        self.pdb.iterate_trajectories( self.scharge, self.efield, self.bfield );

        # std::vector<trajectory_diagnostic_e> diagnostics;
        # diagnostics.push_back( DIAG_X )
        # diagnostics.push_back( DIAG_XP )
        # self.pdb.trajectories_at_plane( self.tdata, AXIS_Z, self.geom.max(2)-self.geom.h(), self.diagnostics )
        # Emittance emit( tdata(0).data(), tdata(1).data() )

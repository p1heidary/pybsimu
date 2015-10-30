
cdef class geometry:
	
	cdef Geometry *geom
	
	def __cinit__(self, int a):
		# self.geom  = new Geometry()
		print 'kir'
	
	# def __dealloc__(self):
		# del self.ibs
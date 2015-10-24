
cdef class pybsimu:
	cdef ibsimu *ibs
	def __cinit__(self):
		ibs  = new kernel.ibsimu()
		geom = 
	def __dealloc__(self):
		del ibs
	def run(self, int iteration):
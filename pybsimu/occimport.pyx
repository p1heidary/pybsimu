from OCC.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.TopAbs import TopAbs_IN

class cysolid:

	classifier = BRepClass3d_SolidClassifier()

	def __cinit__(my_shape):
		self.classifier = BRepClass3d_SolidClassifier(my_shape)

	cdef inside(self, pnt):
		
		self.classifier.Perform(pnt, .001)
		
		return self.classifier.State() == TopAbs_IN
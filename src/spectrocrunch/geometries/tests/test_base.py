import unittest

from .. import base
from ...patch import jsonpickle


class test_base(unittest.TestCase):
    def test_serialize(self):
        g1 = base.Base()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.SolidAngle()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.SolidAngle(solidangle=0.8)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.FlatSample()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.FlatSample(anglein=40.1, angleout=-50.0, detector_azimuth=3.0)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.Centric()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.Centric(
            distance=0.1, anglein=40.1, angleout=-50.0, detector_azimuth=3.0
        )
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

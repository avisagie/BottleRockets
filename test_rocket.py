import unittest

from numpy import sin, cos, sqrt, pi
from rocket import Ballistic


class TestBallistic(unittest.TestCase):
    def test_ballistic1(self):
        """
        Simple test in a vaccuum
        """
        radius = 0
        speed = 10
        theta = 65 * 2*pi/360
        ballistic = Ballistic([0,0], [speed*cos(theta), speed*sin(theta)], 0.0,
                              dry_mass=1, 
                              C_drag=0.3,
                              A_cross_sectional_area=pi * (radius**2),
                              timestep=0.001)

        for ii in range(100000):
            time, position, velocity = ballistic.step()
            # print(f'{time:0.03f}: {position}, {velocity}')

            # did it hit the ground?
            if position[1] < 0:
                break

        # analytical solution to zero drag flight aboce a flat plane
        expected_end_time = speed*sin(theta) / (0.5*9.81)
        # print(f"expected end time: {expected_end_time}s")

        if abs(expected_end_time - time) > 0.01:
            raise Exception("Not close enough to end time")

if __name__ == '__main__':
    unittest.main()
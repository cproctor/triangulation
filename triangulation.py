from math import sin, cos, pi, sqrt, acos
from collections import OrderedDict
from random import random, choice

# Abstract geometric concepts

class Point():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __sub__(self, other):
        "The difference between two points is a vector"
        return Vector(self.x - other.x, self.y - other.y)

    def __add__(self, vector):
        assert(isinstance(vector, Vector))
        return Point(self.x + vector.x, self.y + vector.y)

class Vector():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def length(self):
        return sqrt(self.dot(self))

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def scale(self, scalar):
        return Vector(self.x * scalar, self.y * scalar)

    def normalize(self):
        return self.scale(1.0/self.length())

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y)

    def orthogonal(self):
        return Vector(self.y, self.x * -1)

# Concrete classes represent real-world measurements, which may include errors
class MapPoint():
    def __init__(self, x=None, y=None, index=0, name=""):
        self.x = x
        self.y = y
        self.index = 0
        self.name = name
        self.changed = False
        self.edges = []

    def fixed(self):
        "Check whether this point's position has been determined"
        return self.x is not None and self.y is not None

    def compute_position(self):
        """
        Try to determine this point's posiiton. This is only possible if it
        is part of a triangle, where the other two points are fixed.
        First, we project ac onto ab, to get the component vector which is parallel
        to ab. Then we can use the pythagorean formula to compute the length
        of the orthogonal component. 
        """
        assert not self.fixed()

        fixedNeighbors = [n for n in self.neighbors() if n.fixed()]
        for neighbor in fixedNeighbors:
            for other in [n for n in neighbor.neighbors() if n in fixedNeighbors and n is not self]:
                MapTriangle(self, neighbor, other).fix_points()
                return True

    def get_edge_to(self, destination):
        for edge in self.edges:
            if destination in edge.points:
                return edge

    def neighbors(self):
        return [edge.traverse(self) for edge in self.edges]

    def __sub__(self, other):
        "The difference between two points is a vector"
        assert self.fixed() and other.fixed(), "Cannot create a vector between unfixed points"
        return Vector(self.x - other.x, self.y - other.y)

    def __add__(self, vector):
        assert isinstance(vector, Vector)
        return Point(self.x + vector.x, self.y + vector.y)

    def __repr__(self):
        return "<Point {} ({}, {})>".format(self.name,  round(self.x, 2) if self.x is not None else "?", round(self.y, 2) if self.y is not None else "?")

class MapEdge(object):
    def __init__(self, point1, point2, length, index=0):
        if point1 is point2: 
            raise ValueError("Edges from a point to itself are not allowed.")
        self.points = [point1, point2]
        point1.edges.append(self)
        point2.edges.append(self)
        self.length = length # The measured length
        self.index = index

    def error(self):
        if self.fixed():
            computedLength = (self.points[1] - self.points[0]).length()
            return self.length - computedLength
        else: 
            return 0

    def fixed(self):
        return self.points[0].fixed() and self.points[1].fixed()

    def traverse(self, origin):
        if not origin in self.points:
            raise ValueError("Edge {} does not connect to MapPoint {}".format(self, origin))
        return self.points[0] if self.points[0] is not origin else self.points[1]

    def __repr__(self):
        return "<Edge from {} to {} with length {}>".format(self.points[0], self.points[1], self.length)

class MapTriangle():
    def __init__(self, p1, p2, p3):
        if not p1.get_edge_to(p2) and p2.get_edge_to(p3) and p3.get_edge_to(p1): 
            raise ValueError("Triangles must be composed of points linked by edges")

        # Sort elements by their index, so that the most recently added edge will be last.
        self.edges = sorted([p1.get_edge_to(p2), p2.get_edge_to(p3), p3.get_edge_to(p1)], key=lambda e:e.index)
        final_edge = self.edges[2]

        # Define the order of points so that they progress in the order of the last-added edge.
        c = final_edge.points[0]
        a = final_edge.points[1]
        for point in (p1, p2, p3):
            if point is not a and point is not c:
                b = point
        self.points = [a, b, c]

    def follow(self, point):
        "Gets the next point, in counter-clockwise order"
        for i in range(3):
            if point is self.points[i]:
                return self.points[(i + 1) % 3]

    def get_angle(self, point):
        assert point in self.points
        # Let us assume we are looking for angle C
        # It does not matter which point we call A and B, because
        # acos always returns a positive value. 
        others = [p for p in self.points if p is not point]
        a = point.get_edge_to(others[0]).length
        b = point.get_edge_to(others[1]).length
        c = others[0].get_edge_to(others[1]).length

        # law of cosines
        return acos((a*a + b*b - c*c) / (2.0 * a * b))

    def fix_points(self):
        """
        If two of the points are fixed, we can fix the third by computing the 
        angle at one of the fixed points, and using this to compute the length
        of vectors parallel and perpendicular to vector ab. 
        """
        assert len(self.fixed_points()) == 2, "Can only fix points for a triangle with two fixed points. {}".format(fixed)

        # A new attempt, respecting the predecence of edges.   
        # The final edge of a triangle always defines the directionality, 
        # triangles must go counterclockwise. This means that as the final 
        # edge goes from C to A, B is on the left as we go.
        
        # C is always the unfixed point. 
        c = self.unfixed_points()[0]
        a = self.follow(c)
        b = self.follow(a)

        # The final edge determines the directionality of the triangle
        print("FIXING POINTS (B should be to the left of C->A)")
        print("A: {}, B: {}, C: {}".format(a, b, c))
        angle_a = self.get_angle(a)
        ac = a.get_edge_to(c).length

        parallel = (b - a).normalize().scale(ac * cos(angle_a))
        perp = parallel.orthogonal().normalize().scale(-1 * ac * sin(angle_a))
        position = a + parallel + perp

        c.x = position.x
        c.y = position.y
        c.changed = True

    def fixed_points(self):
        return [p for p in self.points if p.fixed()]

    def unfixed_points(self):
        return [p for p in self.points if not p.fixed()]

    def fixed(self):
        return not any(self.unfixed_points())

# A triangulation is constructed from a CSV file containing rows of measurements.
# Each measurement has two point names (point1, point2) and a distance. 
# Measurements should be provided in order so that triangles are constructed in 
# a counter-clockwise direction.
# We define the first edge to go from the origin in the direction of orientation
# (in degrees clockwise from North)
class Triangulation():
    def __init__(self, measurements, orientation=0):
        self.points = OrderedDict()
        self.edges = []
        self.indexCounter = 0
        for i, m in enumerate(measurements):
            p1 = self.get_or_create_point(m['point1'])
            p2 = self.get_or_create_point(m['point2'])
            edge = MapEdge(p1, p2, m['distance'], i)
            self.edges.append(edge)
        origin, endpoint = self.edges[0].points
        origin.x = 0
        origin.y = 0
        origin.changed = True
        endpoint.x = self.edges[0].length * cos(pi/2 - to_radians(orientation))
        endpoint.y = self.edges[0].length * sin(pi/2 - to_radians(orientation))
        endpoint.changed = True

    def compute_positions(self):
        while any(filter(lambda p: p.changed, self.points.values())):
            self.clear_change_flags()
            for point in [p for p in self.points.values() if not p.fixed()]:
                point.compute_position()

    def improve(self, dist=12, iterations = 1000):
        """
        Iteratively try moving points around to see whether we can get 
        MSE down.
        """
        for i in range(iterations):
            point = choice(list(self.points.values()))
            err = sum(pow(e.error(), 2) for e in point.edges)
            direction = random() * 2 * pi
            point.x += dist * sin(direction)
            point.y += dist * cos(direction)
            if err > sum(pow(e.error(), 2) for e in point.edges):
                dist = dist * 0.9
            else: 
                point.x -= dist * sin(direction)
                point.y -= dist * cos(direction)

    def plt_data(self):
        return [[[e.points[0].x, e.points[1].x], [e.points[0].y, e.points[1].y]] for e in self.edges]
    
    def mean_squared_error(self):
        return sum(pow(p.error(), 2) for p in self.edges)/len(self.edges)
            
    def clear_change_flags(self):
        for point in self.points.values():
            point.changed = False
            
    def get_or_create_point(self, pointName):
        if self.points.get(pointName) is None:
            self.points[pointName] = MapPoint(index = self.indexCounter, name=pointName)
            self.indexCounter += 1
        return self.points.get(pointName)
            
# Helpers
def to_radians(degrees):
    return degrees * (2*pi)/360.0

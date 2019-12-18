using System;

namespace protein.Geometry
{
    /** Represents a line in 3D given by a single point on it and its direction vector. */
    public class Line
    {
        public Vector FixedPoint { get; set; }
        public Vector Direction { get; set; }
        public Line(Vector fixedPoint, Vector direction)
        {
            FixedPoint = fixedPoint;
            Direction = direction;
        }
        public Line(Point fixedPoint, Vector direction) : this(fixedPoint.Vector, direction) { }
        public Line(Point point1, Point point2) : this(point1.Vector, (point2.Vector - point1.Vector).Normalize()) { }
    }
}


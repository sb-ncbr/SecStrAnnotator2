using System;

namespace protein.Geometry
{
    /** Represents a line in 3D given by a single point on it and its direction vector. */
    public struct Line
    {
        public readonly Point FixedPoint;
        public readonly Vector Direction;

        public Line(Point fixedPoint, Vector direction)
        {
            FixedPoint = fixedPoint;
            Direction = direction;
        }

        // public Line(Point fixedPoint, Vector direction) : this(fixedPoint.Vector, direction) { }
        
        public Line(Point point1, Point point2) : this(point1, (point2 - point1).Normalize()) { }
    }
}


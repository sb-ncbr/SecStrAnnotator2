using System;

namespace protein.Geometry
{
    /** Represents a plane in 3D given by a single point on it and its normal (perpendicular!) vector. */
    public class Plane
    {
        public Vector FixedPoint { get; set; }
        public Vector Normal { get; set; }
        public Plane(Vector fixedPoint, Vector normal)
        {
            FixedPoint = fixedPoint;
            Normal = normal;
        }
        public Plane(Point fixedPoint, Vector normal) : this(fixedPoint.Vector, normal) { }
    }

}


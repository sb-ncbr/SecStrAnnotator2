using System;

namespace protein.Geometry
{
    /** Represents a plane in 3D given by a single point on it and its normal (perpendicular!) vector. */
    public struct Plane
    {
        public readonly Point FixedPoint;
        public readonly Vector Normal;

        public Plane(Point fixedPoint, Vector normal)
        {
            FixedPoint = fixedPoint;
            Normal = normal;
        }
    }
}


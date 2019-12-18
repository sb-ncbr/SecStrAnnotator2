using System;

namespace protein.Geometry
{
    /** Represents a point in 3D with position Vector. */
    public class Point
    {
        public Vector Vector { get; set; }
        public Point(Vector vector)
        {
            Vector = vector;
        }
    }
}


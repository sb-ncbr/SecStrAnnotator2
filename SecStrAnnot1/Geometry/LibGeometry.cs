using System;

namespace protein.Geometry
{
    public static class LibGeometry
    {
        /** Calculates the angle between two vectors (in [0,pi] radians). */
        public static double AngleInRadians(Vector u, Vector v)
        {
            double dot = u * v;
            double absDot = Math.Abs(dot);
            int sgnDot = Math.Sign(dot);
            double absCross = (u % v).Size;
            double result = (absCross < absDot) ? Math.Asin(absCross / u.Size / v.Size) : Math.Acos(absDot / u.Size / v.Size);
            if (sgnDot < 0)
                result = Math.PI - result;
            return result;
        }

        /** Calculates the angle between two vectors (in [0,180] degrees). */
        public static double AngleInDegrees(Vector u, Vector v)
        {
            return Radians2Degrees(AngleInRadians(u, v));
        }

        /** Converts an angle in radians to degrees. */
        public static double Radians2Degrees(double radians)
        {
            return 180 * radians / Math.PI;
        }

        public static bool ArePerpendicular(Vector u, Vector v)
        {
            return (u.Normalize() * v.Normalize()) < Vector.EPSILON;
        }

        public static bool AreParallel(Vector u, Vector v)
        {
            return (u.Normalize() % v.Normalize()).IsZero();
        }
        
        public static bool AreParallel(Line line, Plane plane)
        {
            return ArePerpendicular(line.Direction, plane.Normal);
        }

        /** Returns a point that is common to two given geometrical objects, or null if such point does not exist. */
        public static Point? Intersection(Line line, Plane plane)
        {
            if (AreParallel(line, plane))
                return null;
            double t = ((plane.FixedPoint - line.FixedPoint) * plane.Normal) / (line.Direction * plane.Normal);
            return line.FixedPoint + t * line.Direction;
        }

        /* Returns the point that is in the middle between two points.*/
        public static Point Middle(Point a, Point b)
        {
            return new Point(0.5 * (a.Vector + b.Vector));
        }

        public static Point ProjectPointOnLine(Point r, Point linePoint, Vector lineDirection)
        {
            lineDirection = lineDirection.Normalize();
            return linePoint + ((r - linePoint) * lineDirection) * lineDirection;
        }

        /** Returns the distance between two objects. */
        public static double Distance(Point a, Point b)
        {
            return (a- b).Size;
        }
    }
}


using System;

namespace protein.Geometry
{
    /** Represents a point in 3D with position Vector. */
    public struct Point
    {
        public readonly Vector Vector;  

        public Point(Vector vector)
        {
            Vector = vector;
        }

        public Point(double x, double y, double z)
        {
            Vector = new Vector(x, y, z);
        }

        public static Point operator +(Point a, Vector v) => new Point(a.Vector + v);

        public static Vector operator -(Point a, Point b) => a.Vector - b.Vector;

        public double X => Vector.X;
        public double Y => Vector.Y;
        public double Z => Vector.Z;

        public Point Transform(Matrix rotation, Matrix translation){
            Matrix oldMatrix = Matrix.FromRows(new Point[]{this});
            Matrix newMatrix = oldMatrix * rotation + translation;
            Point newPoint = new Point(newMatrix.ToRowVectors()[0]);
            return newPoint;
        }

        public override String ToString()
        {
            return String.Format("({0},{1},{2})", X, Y, Z);
        }
        public String ToString(String format)
        {
            return String.Format("({0},{1},{2})", X.ToString(format), Y.ToString(format), Z.ToString(format));
        }
    }
}


using System;
using System.Collections.Generic;

namespace protein
{
	public struct Vector
	{
		public double X;
		public double Y;
		public double Z;

		public const double EPSILON = 10e-9;
		public bool IsZero (){
			return this.Size < EPSILON;
		}

		public Vector (double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		public override String ToString(){
			return String.Format ("[{0},{1},{2}]", X, Y, Z);
		}
		public String ToString(String format){
			return String.Format ("[{0},{1},{2}]", X.ToString (format), Y.ToString (format), Z.ToString (format));
		}

		public static Vector ZERO { get { return new Vector (0, 0, 0); } }

		public static Vector operator +(Vector u, Vector v){
			return new Vector (u.X + v.X, u.Y + v.Y, u.Z + v.Z);
		}

		public static Vector operator -(Vector v){
			return new Vector (-v.X, -v.Y, -v.Z);
		}

		public static Vector operator -(Vector u, Vector v){
			return u + (-v);
		}

		public static Vector operator *(double a, Vector v){
			return new Vector (a * v.X, a * v.Y, a * v.Z);
		}

		public static Vector operator /(Vector v, double a){
			return new Vector (v.X/a, v.Y/a, v.Z/a);
		}

		/*Dot product*/
		public static double operator *(Vector u, Vector v){
			return u.X * v.X + u.Y * v.Y + u.Z * v.Z;
		}


		/*Cross product*/
		public static Vector operator %(Vector u, Vector v){
			return new Vector (u.Y * v.Z - u.Z * v.Y, u.Z * v.X - u.X * v.Z, u.X * v.Y - u.Y * v.X);
		}

		public double Size{ get { return Math.Sqrt (this * this); } }

		public double SqSize{ get { return this * this; } }

		public Vector Normalize(){
			double s = this.Size;
			try {
				return 1/s*this;
			} catch (DivideByZeroException e) {
				throw new DivideByZeroException ("Trying to normalize zero vector.", e);
			}
		}

		public static Vector Sum(IEnumerable<Vector> vectors){
			Vector result = new Vector (0, 0, 0);
			foreach (Vector v in vectors)
				result += v;
			return result;
		}

		public static Vector Average(IEnumerable<Vector> vectors){
			Vector sum = new Vector (0, 0, 0);
			int num = 0;
			foreach (Vector v in vectors) {
				sum += v;
				num++;
			}
			return (1.0 / num) * sum;
		}

		public List<object> AsList(){
			return new List<object>{ X, Y, Z };
		}

		public double[] AsArray(){
			return new double[]{ X, Y, Z };
		}

		public Vector Round(int digits){
			return new Vector (Math.Round (X, digits), Math.Round (Y, digits), Math.Round (Z, digits));
		}

	}
}


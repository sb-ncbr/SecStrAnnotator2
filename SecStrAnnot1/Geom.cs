using System;

namespace protein
{
	public static class Geom
	{
		/** Represents a point in 3D with position Vector. */
		public class Point{
			public Vector Vector {get;set;}
			public Point (Vector vector){
				Vector=vector;
			}
		}

		/** Represents a line in 3D given by a single point on it and its direction vector. */
		public class Line{
			public Vector FixedPoint {get;set;}
			public Vector Direction{get;set;}
			public Line (Vector fixedPoint, Vector direction){
				FixedPoint=fixedPoint;
				Direction=direction;
			}
			public Line (Point fixedPoint, Vector direction):this(fixedPoint.Vector,direction){}
			public Line (Point point1, Point point2):this(point1.Vector,(point2.Vector-point1.Vector).Normalize()){}
		}

		/** Represents a line segment in 3D given by its start and end point. */
		public class LineSegment{
			public Vector Start {get;set;}
			public Vector End{get;set;}
			public LineSegment (Vector start, Vector end){
				Start=start;
				End=end;
			}
			public LineSegment(Point start, Point end):this(start.Vector,end.Vector){}
			public Vector Direction{get{ return (End-Start).Normalize(); }}
			/** Returns a line that is created by enlongation of this line segment. */
			public Line ToLine(){
				return new Line (Start, Direction);
			}
		}

		/** Represents a plane in 3D given by a single point on it and its normal (perpendicular!) vector. */
		public class Plane{
			public Vector FixedPoint {get;set;}
			public Vector Normal{get;set;}
			public Plane (Vector fixedPoint, Vector normal){
				FixedPoint=fixedPoint;
				Normal=normal;
			}
			public Plane (Point fixedPoint, Vector normal):this(fixedPoint.Vector,normal){}
		}

		/** Calculates the angle between two vectors (in [0,pi] radians). */
		public static double AngleInRadians(Vector u, Vector v){
			double dot = u * v;
			double absDot = Math.Abs (dot);
			int sgnDot = Math.Sign (dot);
			double absCross = (u % v).Size;
			double result = (absCross < absDot) ? Math.Asin (absCross / u.Size / v.Size) : Math.Acos (absDot / u.Size / v.Size);
			if (sgnDot < 0)
				result = Math.PI - result;
			return result;
		}

		/** Calculates the angle between two vectors (in [0,180] degrees). */
		public static double AngleInDegrees(Vector u, Vector v){
			return Radians2Degrees (AngleInRadians (u, v));
		}

		/** Converts an angle in radians to degrees. */
		public static double Radians2Degrees(double radians){
			return 180 * radians / Math.PI;
		}

		public static bool ArePerpendicular(Vector u, Vector v){
			return (u.Normalize() * v.Normalize()) < Vector.EPSILON;
		}
		public static bool AreParallel(Vector u, Vector v){
			return (u.Normalize() % v.Normalize()).IsZero();
		}
		public static bool AreParallel(Line line, Plane plane){
			return ArePerpendicular(line.Direction,plane.Normal);
		}

		/** Returns a point that is common to two given geometrical objects, or null if such point does not exist. */
		public static Point Intersection (Line line, Plane plane){
			if (AreParallel(line,plane))
				return null;
			double t = ((plane.FixedPoint - line.FixedPoint) * plane.Normal) / (line.Direction * plane.Normal);
			return new Point (line.FixedPoint + t * line.Direction);
		}

		/* Returns the point that is in the middle between two points.*/
		public static Point Middle(Point a, Point b){
			return new Point (0.5*(a.Vector + b.Vector));
		}

		public static Vector ProjectPointOnLine (Vector r, Vector linePoint, Vector lineDirection){
			lineDirection = lineDirection.Normalize ();
			return linePoint + ((r - linePoint) * lineDirection) * lineDirection;
		}

		/** Returns the distance between two objects. */
		public static double Distance(Point a, Point b){
			return (a.Vector - b.Vector).Size;
		}
	}
}


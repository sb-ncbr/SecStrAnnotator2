using System;
using System.Collections.Generic;
using System.Linq;

namespace protein
{
	public class Superimposer
	{
		public Superimposer ()
		{
		}

		public static List<Atom> RotateAll (ICollection<Atom> atoms, double[,] rot){
			List<Atom> result = new List<Atom> ();
			foreach (Atom a in atoms) {
				Vector v = Times (rot, a.Position());
				result.Add (new Atom (a.Serial, a.Name, a.AltLoc, a.ResName, a.ChainID, a.ResSeq, a.ICode, 
				                      v.X, v.Y, v.Z, a.Occupancy, a.TempFactor, a.Element, a.Charge, a.IsHetAtm));
			}
			return result;
		}

		public static List<Atom> MoveAll (ICollection<Atom> atoms, Vector shift){
			List<Atom> result = new List<Atom> ();
			foreach (Atom a in atoms) {
				Vector v = a.Position () + shift;
				result.Add (new Atom (a.Serial, a.Name, a.AltLoc, a.ResName, a.ChainID, a.ResSeq, a.ICode, 
				                      v.X, v.Y, v.Z, a.Occupancy, a.TempFactor, a.Element, a.Charge, a.IsHetAtm));
			}
			return result;
		}

		public static Tuple<double,double> ThetaPhi(Vector v){
			//Asin returns in <-pi/2, pi/2>
			//Acos returns in <0, pi>

			double xy = Math.Sqrt (v.X * v.X + v.Y * v.Y);
			double r = v.Size;
			double phi = Math.Abs (v.X) > Math.Abs (v.Y) ?
				v.X > 0 ?
					Math.Asin (v.Y / xy)
					: Math.PI - Math.Asin (v.Y / xy)
				: v.Y > 0 ?
					Math.Acos (v.X / xy)
					: 2 * Math.PI - Math.Acos (v.X / xy);
			if (phi < 0)
				phi += 2 * Math.PI;
			double theta = Math.Abs (v.Z) > xy ?
					v.Z > 0 ?
						Math.Asin (xy / r)
						: Math.PI - Math.Asin (xy / r)
					: Math.Acos (v.Z / r);

			return new Tuple<double,double> (theta, phi);
		}

		/*Vector after extrinsic rotation around x, y, and then z-axis.*/
		/*private Vector Rotate(Vector v, double phiX,double phiY,double phiZ){
			//rotation around X
			double cos = Math.Cos (phiX);
			double sin = Math.Sin (phiX);
			v = new Vector (v.X,
			                cos * v.Y - sin * v.Z,
			                sin * v.Y + cos * v.Z);
			//rotation around Y
			cos = Math.Cos (phiY);
			sin = Math.Sin (phiY);
			v = new Vector (cos * v.X + sin * v.Z,
			                v.Y,
			                -sin * v.X + cos * v.Z);
			//rotation around Z
			cos = Math.Cos (phiZ);
			sin = Math.Sin (phiZ);
			v = new Vector (cos * v.X - sin * v.Y,
			                sin * v.X + cos * v.Y,
			                v.Z);
			return v;
		}
		*/

		private static double[,] RotateX(double phi){
			double cos=Math.Cos(phi);
			double sin = Math.Sin (phi);
			return new double[3, 3] { {1,0,0}, {0,cos,-sin}, {0,sin,cos} };
		}
		private static double[,] RotateY(double phi){
			double cos=Math.Cos(phi);
			double sin = Math.Sin (phi);
			return new double[3, 3] { {cos,0,sin}, {0,1,0}, {-sin,0,cos} };
		}
		private static double[,] RotateZ(double phi){
			double cos=Math.Cos(phi);
			double sin = Math.Sin (phi);
			return new double[3, 3] { {cos,-sin,0}, {sin,cos,0}, {0,0,1} };
		}
		private static double[,] Times(double[,] a, double[,] b){
			double[,] result=new double[,]{{0,0,0},{0,0,0},{0,0,0}};
			for (int i=0; i<=2; i++)
				for (int j=0; j<=2; j++)
					result [i, j] = a [i, 0] * b [0, j] + a [i, 1] * b [1, j] + a [i, 2] * b [2, j];
			return result;
		}
		public static Vector Times(double[,] a, Vector v){
			double[] result=new double[]{0,0,0};
			for (int i=0; i<=2; i++)
				result [i] = a [i, 0] * v.X + a [i, 1] * v.Y + a [i, 2] * v.Z;
			return new Vector (result [0], result [1], result [2]);
		}

		public static double[,] RandomRot(Random r){
			return Times (RotateX(2*Math.PI*r.NextDouble()), Times (RotateY(2*Math.PI*r.NextDouble()), RotateZ(2*Math.PI*r.NextDouble())));}

		/*Find such rotational matrix R that R*u1=v1 and R*u2=v2.
		 Assumes that angle u1-u2 = v1-v2 (otherwise no solution exists), but does not check it.
		 This function is magic, don't try to understand it.*/
		public static double[,] FindRotation(Vector u1,Vector u2, Vector v1, Vector v2){
			Tuple<double,double> tpU1=ThetaPhi(u1);
			double[,] R1 = Times (RotateY (-tpU1.Item1), RotateZ (-tpU1.Item2));
			Tuple<double,double> tpV1 = ThetaPhi (v1);
			double[,] R2 = Times (RotateY (-tpV1.Item1), RotateZ (-tpV1.Item2));
			double[,] R2inv = Times (RotateZ (tpV1.Item2), RotateY (tpV1.Item1));
			Tuple<double,double> tpS = ThetaPhi (Times (R1, u2));
			Tuple<double,double> tpT = ThetaPhi (Times (R2, v2));
			double[,] R3 = RotateZ (tpT.Item2 - tpS.Item2);
			return Times (R2inv, Times (R3, R1));
		}

		public static List<Atom> Superimpose1 (ICollection<Atom> atoms, char chainID, 
		ICollection<Atom> reference, char chainIDref, 
		String name1, String name2, String name3){

			Vector a1 = atoms.First (delegate(Atom a) {
				return a.Name == name1 & a.ChainID == chainID;
			}).Position ();
			Vector a2 = atoms.First (delegate(Atom a) {
				return a.Name == name2 & a.ChainID == chainID;
			}).Position ();
			Vector a3 = atoms.First (delegate(Atom a) {
				return a.Name == name3 & a.ChainID == chainID;
			}).Position ();
			Vector r1 = reference.First (delegate(Atom a) {
				return a.Name == name1 & a.ChainID == chainIDref;
			}).Position ();
			Vector r2 = reference.First (delegate(Atom a) {
				return a.Name == name2 & a.ChainID == chainIDref;
			}).Position ();
			Vector r3 = reference.First (delegate(Atom a) {
				return a.Name == name3 & a.ChainID == chainIDref;
			}).Position ();

			return MoveAll (RotateAll (MoveAll (atoms, -a1), FindRotation (a2 - a1, a3 - a1, r2 - r1, r3 - r1)), r1);
		}
	}
}


using System;
using System.Collections.Generic;
using System.Linq;

using protein.Geometry;

namespace protein.Libraries
{
	public static class LibAlgebra
	{
		public static double Epsilon { get { return 1e-9; } }
		public static bool IsZero(double x){
			return Math.Abs (x) < Epsilon;
		}

		/**Solves a cubic equation c3*x^3 + c2*x^2 + c1*x + c0 = 0. 
		 * Returns 1, 2, or 3 REAL roots.
		 * Uses the method described at http://en.citizendium.org/wiki/Cubic_equation/Proofs */
		public static List<double> SolveCubic(double c3, double c2, double c1, double c0){
			// 

			//Console.WriteLine ("Solving equation: {0} x^3 + {1} x^2 + {2} x + {3} = 0", c3, c2, c1, c0); 

			try {
				double _ = 1/c3;
			} catch (DivideByZeroException e) {
				throw new DivideByZeroException ("Solving cubic equation with coefficient of cubic term == 0.", e);
			}
			double a = c2 / c3;
			double b = c1 / c3;
			double c = c0 / c3;

			double P = (3 * b - a * a) / 3;
			double Q = (27 * c - 9 * a * b + 2 * a * a * a) / 27;

			//Console.WriteLine ("P = {0}; Q = {1}", P, Q);

			double D = Q * Q + 4 * P * P * P / 27;
			
			double[] w3Re = new double[2];
			double[] w3Im = new double[2];
			double[] wRe = new double[6];
			double[] wIm = new double[6];

			if (D >= 0) {
				w3Re [0] = (-Q + Math.Sqrt (D)) / 2;
				w3Im [0] = 0;
				w3Re [1] = (-Q - Math.Sqrt (D)) / 2;
				w3Im [1] = 0;
			} else {
				w3Re [0] = (-Q) / 2;
				w3Im [0] = (Math.Sqrt (-D)) / 2;
				w3Re [1] = (-Q) / 2;
				w3Im [1] = (-Math.Sqrt (-D)) / 2;
			}
			
			//Console.WriteLine ("w3_1 = {0} + {1} i", w3Re[0], w3Im[0]);
			//Console.WriteLine ("w3_2 = {0} + {1} i", w3Re[1], w3Im[1]);

			//Error is somewhere here.

			for (int j=0; j<=1; j++) {
				double abs = Math.Pow (w3Re [j] * w3Re [j] + w3Im [j] * w3Im [j], 1.0 / 6.0);
				double phase = Math.Sign (w3Im [j]) * LibGeometry.AngleInRadians (new Vector (w3Re [j], w3Im [j], 0), new Vector (1.0, 0, 0)) / 3;
				wRe [3 * j] = abs * Math.Cos (phase);	wIm [3 * j] = abs * Math.Sin (phase);
				wRe [3 * j + 1] = -0.5 * wRe [3 * j] - 0.5 * Math.Sqrt (3.0) * wIm [3 * j];
				wIm [3 * j + 1] = -0.5 * wIm [3 * j] + 0.5 * Math.Sqrt (3.0) * wRe [3 * j];
				wRe [3 * j + 2] = -0.5 * wRe [3 * j] + 0.5 * Math.Sqrt (3.0) * wIm [3 * j];
				wIm [3 * j + 2] = -0.5 * wIm [3 * j] - 0.5 * Math.Sqrt (3.0) * wRe [3 * j];
			}

			double[] xRe = new double[6];
			double[] xIm = new double[6];

			for (int j=0; j<=5; j++) {
				//Console.WriteLine ("w_{0} = {1} + {2} i", j, wRe [j], wIm [j]);
				double denom = 3.0 * (wRe [j] * wRe [j] + wIm [j] * wIm [j]);
				xRe [j] = wRe [j] - P * wRe [j] / denom - a / 3.0;
				xIm [j] = wIm [j] + P * wIm [j] / denom;
			}

			bool[] output = new bool[6] { true, true, true, true, true, true };
			List<double> result = new List<double> ();
			for (int j=0; j<=5; j++) {
				//Console.WriteLine ("x_{0} = {1} + {2} i", j, xRe [j], xIm [j]);
				if (output [j] && IsZero(xIm [j])) {
					result.Add (xRe [j]);
					for (int i=j+1; i<=5; i++) {
						if (IsZero(xRe [i] - xRe [j]))
							output [i] = false;
					}
				}
			}

			if (result.Count == 0)
				throw new Exception ("SolveCubic: No real roots were found, but every cubic equation must have a real root. Something very wrong must have happened with the universe.");

			/*Console.WriteLine ("Found {0} roots.", result.Count);
			for (int i=0; i<result.Count; i++) {
				double value = c3 * result [i] * result [i] * result [i] + c2 * result [i] * result [i] + c1 * result [i] + c0;
				Console.WriteLine ("x = {0}; f(x) = {1}", result [i], value);
			}*/
			//Console.WriteLine ("---------------------------------------------------");
			

			return result;
		}


		/**Returns the eigenvectors and their corresponding eigenvalues of 3x3 matrix A with rows v1, v2, v3. */
		public static List<(Vector, double)> Eigenvectors(Matrix A){

			// Scaling so that we will work with numbers close to 1
			double scale = Math.Sqrt (A.SumSq); 
			A = (1.0 / scale) * A;

			//determinant of the matrix A-lambda*E is c3*lambda^3 + c2*lambda^2 + c1*lambda + c0 = 0
			double c3 = -1.0;
			double c2 = A[0,0] + A[1,1] + A[2,2];
			double c1 = -A[0,0] * A[1,1] - A[0,0] * A[2,2] - A[1,1] * A[2,2] + A[1,2] * A[2,1] + A[0,2] * A[2,0] + A[0,1] * A[1,0];
			double c0 = A[0,0] * A[1,1] * A[2,2] + A[0,2] * A[1,0] * A[2,1] + A[0,1] * A[1,2] * A[2,0]
				- A[0,0] * A[1,2] * A[2,1] - A[0,2] * A[1,1] * A[2,0] - A[0,1] * A[1,0] * A[2,2];
			
			List<double> lambdas = SolveCubic (c3, c2, c1, c0);
			lambdas.Sort ((double a, double b) => Math.Sign (b * b - a * a));
			//lambdas.ForEach (x => Console.WriteLine ("lambda = {0}", x));


			List<(Vector, double)> result = new List<(Vector, double)> ();

			//rows of the matrix B = A-lambda*E
			for (int i = 0; i < 3; i++) {
				Vector[] B = A.ToRowVectors ().ToArray ();
				B [0] -= new Vector (lambdas [i], 0, 0);
				B [1] -= new Vector (0, lambdas [i], 0);
				B [2] -= new Vector (0, 0, lambdas [i]);
				Vector[] eigenVecs = SolveHomo3 (B);
				//Console.WriteLine ("For lambda {0} we found {1} eigenvectors.", lambdas [i], eigenVecs.Length);
				foreach (Vector v in eigenVecs) {
					result.Add ((v, lambdas [i]));
				}
			}

			return result.Select (x => (x.Item1, scale * x.Item2)).ToList(); //Reversing the scaling
		}

		/**Solves homogenous system of 3 linear equations with 3 variables. */
		public static Vector[] SolveHomo3(Vector[] w){
			if (w.Length != 3)
				throw new ArgumentException ("SolveHomo3: Argument must be an array of 3 Vectors.");

			if (w [0].IsZero () && w [1].IsZero () && w [2].IsZero ()) {
				//Console.WriteLine ("Three solutions.");
				// rank 0, three solutions
				return new Vector[3] {
					new Vector (1, 0, 0),
					new Vector (0, 1, 0),
					new Vector (0, 0, 1)
				};
			} else if ((w [0] % w [1]).IsZero () && (w [0] % w [2]).IsZero () && (w [1] % w [2]).IsZero ()) {
				//Console.WriteLine ("Two solutions.");
				// rank 1, two solutions
				Vector p = w [0];
				if (w [1].SqSize > p.SqSize)
					p = w [1];
				if (w [2].SqSize > p.SqSize)
					p = w [2];
				p = p.Normalize ();
				Vector[] res = new Vector[2];
				if (p.Z * p.Z < p.X * p.X + p.Y * p.Y) {
					res[0] = (new Vector (0, 0, 1) % p).Normalize ();
					res[1] = (res[0] % p).Normalize ();
				} else {
					res[0] = (new Vector (1, 0, 0) % p).Normalize ();
					res[1] = (res[0] % p).Normalize ();					
				}
				return res;
			} else if (IsZero ((w [0] % w [1]) * w [2])) {
				//Console.WriteLine ("One solution.");
				//rank 2, one solution
				Vector r = w [0] % w [1];
				if ((w [0] % w [2]).SqSize > r.SqSize)
					r = w [0] % w [2];
				if ((w [1] % w [2]).SqSize > r.SqSize)
					r = w [1] % w [2];
				return new Vector[1]{r.Normalize ()};
			} else {
				//Console.WriteLine ("No solutions.");
				// rank 3 (regular) matrix, no solutions
				return new Vector[0]{};
			}
		}

		/**Performs Singular Value Decomposition of a n*3 matrix Y and returns U, S, and V such that U*S*V=Y. */
		public static (Matrix, Matrix.DiagMatrix, Matrix) SVD(Matrix Y){
			List<(Vector, double)> eigenVecs = Eigenvectors (Y.Transpose()*Y);
			if (eigenVecs.Any(x => x.Item2 < 0)) {
				Lib.WriteWarning($"SVD: Negative eigenvalue of covariance matrix: {eigenVecs.Select(x => x.Item2).EnumerateWithCommas()}");
			}
			Matrix.DiagMatrix S = new Matrix.DiagMatrix (eigenVecs.Select (x => Math.Sqrt(Math.Abs(x.Item2))));
			Matrix V = Matrix.FromRowVectors(eigenVecs.Select (x => x.Item1).ToList());
			// so far so good
			Matrix U = Y * V.Transpose () * S.Inverse();
			// Console.WriteLine ($"V = {V}\nV transpose = {V.Transpose()}\nS = {S}\nS inverse = {S.Inverse()}\nU = {U}\nEigenvecs = {eigenVecs.EnumerateWithCommas()}");
			return (U, S, V);
		}

		public static Matrix FitRotation (Matrix Mobile, Matrix Target){
			if (Mobile.Columns != 3 || Target.Columns != 3)
				throw new ArgumentException ("Matrices must have 3 columns.");
			if (Mobile.Rows != Target.Rows)
				throw new ArgumentException ("Matrices must have the same number of rows.");
			
			Matrix H = Mobile.Transpose () * Target;
			(Matrix U, Matrix.DiagMatrix S, Matrix V) = SVD (H);
			Matrix R = U*V;
			// Lib.WriteLineDebug($"U: {U}, V: {V}");
			if (R.Determinant3 () < 0) { 
				//Correction for reflexion matrix
				R = (U * new Matrix.DiagMatrix (new double[]{ 1, 1, -1 })) * V;
			}
			return R;
		}

		public static double RMSDPerRow(Matrix A, Matrix B){
			if (A.Columns != B.Columns)
				throw new ArgumentException ("Matrices must have the same number of columns.");
			if (A.Rows != B.Rows)
				throw new ArgumentException ("Matrices must have the same number of rows.");
			return Math.Sqrt ((A - B).SumSq / A.Rows);
		}

		/**Finds the optimal translation and rotation of mobile to get aligned with target.
		 * In outparams returns translation vector, rotation matrix, and final RMSD.
		 * Expects target to be mean-centered (i.e. centroid = [0,0,0]). */
		public static void Align (Matrix mobile, Matrix target, out Matrix translation, out Matrix rotation, out double RMSD){
			Matrix centeredMobile = mobile.Copy ();
			centeredMobile.MeanCenterVertically (out translation);
			translation = -translation;
			rotation = FitRotation (centeredMobile, target);
			RMSD = RMSDPerRow (centeredMobile * rotation, target);
		}

	}
}


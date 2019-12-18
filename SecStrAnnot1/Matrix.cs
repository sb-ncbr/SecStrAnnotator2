using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

using protein.Geometry;

namespace protein
{
    public class Matrix
    {
        private double[,] matrix;
        public double this[int i, int j] { get { return matrix[i, j]; } set { matrix[i, j] = value; } }
        public int Rows { get { return matrix.GetLength(0); } }
        public int Columns { get { return matrix.GetLength(1); } }

        public double SumSq
        {
            get
            {
                double sum = 0;
                for (int i = 0; i < this.Rows; i++)
                    for (int j = 0; j < this.Columns; j++)
                        sum += this[i, j] * this[i, j];
                return sum;
            }
        }

        private Matrix(int rows, int cols)
        {
            matrix = new double[rows, cols];
            // Hope that auto-initialization works.
            if (rows != Rows || cols != Columns)
                throw new Exception("WTF exception");
        }

        public static Matrix Zeros(int rows, int columns)
        {
            return new Matrix(rows, columns);
        }



        public static Matrix FromRowVectors(IEnumerable<Vector> vectors)
        {
            Matrix result = new Matrix(vectors.Count(), 3);
            int i = 0;
            foreach (Vector v in vectors)
            {
                result[i, 0] = v.X;
                result[i, 1] = v.Y;
                result[i, 2] = v.Z;
                i++;
            }
            return result;
        }

        public static Matrix CreateByColumns(int rows, int columns, ICollection<double> values)
        {
            if (rows * columns != values.Count)
                throw new ArgumentException("Number of values must be rows*columns.");
            Matrix result = new Matrix(rows, columns);
            int i = 0;
            int j = 0;
            foreach (double x in values)
            {
                result[i, j] = x;
                if (i == rows - 1)
                {
                    i = 0;
                    j++;
                }
                else
                {
                    i++;
                }
            }
            return result;
        }

        public static Matrix CreateByRows(int rows, int columns, ICollection<double> values)
        {
            if (rows * columns != values.Count)
                throw new ArgumentException("Number of values must be rows*columns.");
            Matrix result = new Matrix(rows, columns);
            int i = 0;
            int j = 0;
            foreach (double x in values)
            {
                result[i, j] = x;
                if (j == columns - 1)
                {
                    j = 0;
                    i++;
                }
                else
                {
                    j++;
                }
            }
            return result;
        }

        public void NormalizeRows()
        {
            for (int i = 0; i < Rows; i++)
            {
                double sumSq = 0;
                for (int j = 0; j < Columns; j++)
                {
                    sumSq += (this[i, j] * this[i, j]);
                }
                double q = 1.0 / Math.Sqrt(sumSq);
                for (int j = 0; j < Columns; j++)
                {
                    this[i, j] *= q;
                }
            }
        }

        /**Performs the mean-centering of the matrix (subtracts mean column from each column) and returns the mean column of the original matrix.*/
        public void MeanCenterHorizontally(out Matrix meanColumn)
        {
            meanColumn = new Matrix(this.Rows, 1);
            for (int i = 0; i < Rows; i++)
            {
                double sum = 0;
                for (int j = 0; j < Columns; j++)
                {
                    sum += this[i, j];
                }
                double mean = sum / Columns;
                for (int j = 0; j < Columns; j++)
                {
                    this[i, j] -= mean;
                }
                meanColumn[i, 0] = mean;
            }
        }
        public void MeanCenterHorizontally() { Matrix dontCare; MeanCenterHorizontally(out dontCare); }


        /**Performs the mean-centering of the matrix (subtracts mean row from each row) and returns the mean row of the original matrix.*/
        public void MeanCenterVertically(out Matrix meanRow)
        {
            meanRow = new Matrix(1, this.Columns);
            for (int j = 0; j < Columns; j++)
            {
                double sum = 0;
                for (int i = 0; i < Rows; i++)
                {
                    sum += this[i, j];
                }
                double mean = sum / Rows;
                for (int i = 0; i < Rows; i++)
                {
                    this[i, j] -= mean;
                }
                meanRow[0, j] = mean;
            }
        }
        public void MeanCenterVertically() { Matrix dontCare; MeanCenterVertically(out dontCare); }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.Columns != B.Rows)
                throw new ArgumentException("A.Columns!=B.Rows");
            Matrix C = new Matrix(A.Rows, B.Columns);
            for (int i = 0; i < C.Rows; i++)
                for (int j = 0; j < C.Columns; j++)
                    for (int k = 0; k < A.Columns; k++)
                        C[i, j] += A[i, k] * B[k, j];
            return C;
        }

        public static Matrix operator *(double a, Matrix B)
        {
            Matrix C = new Matrix(B.Rows, B.Columns);
            for (int i = 0; i < B.Rows; i++)
                for (int j = 0; j < B.Columns; j++)
                    C[i, j] = a * B[i, j];
            return C;
        }

        public void Add(Matrix B)
        {
            if (this.Columns != B.Columns)
                throw new ArgumentException("Matrices must have the same number of columns.");
            if (this.Rows != B.Rows)
                throw new ArgumentException("Matrices must have the same number of rows.");
            for (int i = 0; i < this.Rows; i++)
                for (int j = 0; j < this.Columns; j++)
                    this[i, j] += B[i, j];
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.Columns != B.Columns)
                throw new ArgumentException("Matrices must have the same number of columns.");
            if (A.Rows != B.Rows)
                throw new ArgumentException("Matrices must have the same number of rows.");
            Matrix C = new Matrix(A.Rows, A.Columns);
            for (int i = 0; i < A.Rows; i++)
                for (int j = 0; j < A.Columns; j++)
                    C[i, j] = A[i, j] - B[i, j];
            return C;
        }

        public static Matrix operator -(Matrix B)
        {
            Matrix C = new Matrix(B.Rows, B.Columns);
            for (int i = 0; i < B.Rows; i++)
                for (int j = 0; j < B.Columns; j++)
                    C[i, j] = -B[i, j];
            return C;
        }

        public List<Vector> ToRowVectors()
        {
            if (this.Columns != 3)
                throw new InvalidOperationException("Columns!=3");
            List<Vector> result = new List<Vector>();
            for (int i = 0; i < this.Rows; i++)
                result.Add(new Vector(this[i, 0], this[i, 1], this[i, 2]));
            return result;
        }

        public Matrix Transpose()
        {
            Matrix C = new Matrix(this.Columns, this.Rows);
            for (int i = 0; i < this.Rows; i++)
                for (int j = 0; j < this.Columns; j++)
                    C[j, i] = this[i, j];
            return C;
        }

        public Matrix Copy()
        {
            Matrix C = new Matrix(this.Rows, this.Columns);
            for (int i = 0; i < this.Rows; i++)
                for (int j = 0; j < this.Columns; j++)
                    C[i, j] = this[i, j];
            return C;
        }

        public double[] GetRow(int i)
        {
            double[] result = new double[this.Columns];
            for (int j = 0; j < this.Columns; j++)
                result[j] = this[i, j];
            return result;
        }

        public Matrix GetRowAsMatrix(int i)
        {
            return Matrix.CreateByRows(1, this.Columns, this.GetRow(i));
        }

        public double Determinant3()
        {
            if (this.Rows != 3 || this.Columns != 3)
                throw new ArgumentException("The matrix must be 3x3.");
            List<Vector> v = this.ToRowVectors();
            return (v[0] % v[1]) * v[2];
        }

        public override String ToString()
        {
            return this.ToRowVectors().Aggregate("Matrix:\n", (String x, Vector y) => x + y.ToString() + "\n");
        }

        public class DiagMatrix
        {
            private double[] diag;
            public double this[int i, int j] { get { if (i == j) return diag[i]; else return 0; } }
            public int Rows { get { return diag.Length; } }
            public int Columns { get { return diag.Length; } }
            public double[] Diagonal { get { return diag.ToArray(); } }

            public DiagMatrix(IEnumerable<double> diagonalElements)
            {
                diag = diagonalElements.ToArray();
            }
            public static Matrix operator *(Matrix A, DiagMatrix B)
            {
                if (A.Columns != B.Rows)
                    throw new ArgumentException("A.Columns!=B.Rows");
                Matrix C = new Matrix(A.Rows, B.Columns);
                for (int i = 0; i < C.Rows; i++)
                    for (int j = 0; j < C.Columns; j++)
                        C[i, j] += A[i, j] * B[j, j];
                return C;
            }

            public DiagMatrix Inverse()
            {
                return new DiagMatrix(this.Diagonal.Select(x => 1.0 / x));
            }

            public List<Vector> ToRowVectors()
            {
                if (this.Columns != 3)
                    throw new InvalidOperationException("Columns!=3");
                List<Vector> result = new List<Vector>();
                for (int i = 0; i < this.Rows; i++)
                    result.Add(new Vector(this[i, 0], this[i, 1], this[i, 2]));
                return result;
            }

            public override String ToString()
            {
                return this.diag.Aggregate("Diagonal matrix with diagonal elements:\n[", (String x, double y) => x + y.ToString() + ", ", w => w + "]");
            }

        }

    }
}


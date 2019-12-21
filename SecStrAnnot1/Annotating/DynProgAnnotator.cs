using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.Annotating
{
    class DynProgAnnotator : IAnnotator
    {
        public AnnotationContext Context { get; private set; }
        public DynProgAnnotator(AnnotationContext context)
        {
            context.ValidateOrdering();
            this.Context = context;
        }
        public static DynProgAnnotator New(AnnotationContext context)
        {
            return new DynProgAnnotator(context);
        }

        public List<(int, int)> GetMatching()
        {
            int m = Context.Templates.Length;
            int n = Context.Candidates.Length;
            double[,] metricMatrix = new double[m, n];
            double[,] dynProgMatrix = new double[m + 1, n + 1];
            int[,] reconstrMatrix = new int[m + 1, n + 1];
            double[] skipTemplatePen = Context.Templates.Select(sse => Context.SkipTemplatePenalty(sse)).ToArray();
            double[] skipCandidatePen = Context.Candidates.Select(sse => Context.SkipCandidatePenalty(sse)).ToArray();

            //coding reconstrMatrix: 
            //  -1: optimum was obtained from the left neighbor (skip candidate SSE)
            //   1: optimum was obtained from upper neighbor (skip template SSE)
            //   0: optimum was obtained from diagonal neighbor (pair a template and candidate SSE)

            // Calculation of the matrix of metric.
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    metricMatrix[i, j] = Context.MetricToMin(Context.Templates[i], Context.Candidates[j]);
                }
            }

            // Calculation of dynamic programming matrix.
            dynProgMatrix[0, 0] = 0.0;
            for (int i = 1; i < m + 1; i++)
            {
                reconstrMatrix[i, 0] = 1;
                dynProgMatrix[i, 0] = dynProgMatrix[i - 1, 0] + skipTemplatePen[i - 1];
            }
            for (int j = 1; j < n + 1; j++)
            {
                reconstrMatrix[0, j] = -1;
                dynProgMatrix[0, j] = dynProgMatrix[0, j - 1] + skipCandidatePen[j - 1];
            }
            for (int i = 1; i < m + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    SseInSpace temp = Context.Templates[i - 1];
                    SseInSpace cand = Context.Candidates[j - 1];

                    double valueUp = dynProgMatrix[i - 1, j] + skipTemplatePen[i - 1];
                    double valueLeft = dynProgMatrix[i, j - 1] + skipCandidatePen[j - 1];

                    if (valueUp <= valueLeft)
                    {
                        //skip template
                        reconstrMatrix[i, j] = 1;
                        dynProgMatrix[i, j] = valueUp;
                    }
                    else
                    {
                        //skip candidate
                        reconstrMatrix[i, j] = -1;
                        dynProgMatrix[i, j] = valueLeft;
                    }

                    if (Context.TypeMatching(temp.Type, cand.Type) && (Context.GuideDiscriminator == null || Context.GuideDiscriminator[i - 1, j - 1]))
                    {
                        double valueDiag = dynProgMatrix[i - 1, j - 1] + metricMatrix[i - 1, j - 1];
                        if (valueDiag <= dynProgMatrix[i, j])
                        {
                            //pair template with candidate
                            reconstrMatrix[i, j] = 0;
                            dynProgMatrix[i, j] = valueDiag;
                        }
                    }
                }
            }

            // print dynamic programming matrix and metric matrix
            if (Lib.DoWriteDebug)
            {
                PrintDynProgAndOrReconstructionMatrix(dynProgMatrix, reconstrMatrix);
                AnnotationHelper.PrintMetricMatrix(Context, metricMatrix, "metric_matrix.tsv");
            }

            // Reconstruction of the best solution.
            List<(int, int)> matching = new List<(int, int)>();
            int row = m;
            int col = n;
            while (row != 0 || col != 0)
            {
                switch (reconstrMatrix[row, col])
                {
                    case -1:
                        // skipped candidate SSE
                        col--;
                        break;
                    case 1:
                        // skipped template SSE
                        row--;
                        break;
                    case 0:
                        // paired template and candidate SSE
                        matching.Add((row - 1, col - 1));
                        row--;
                        col--;
                        break;
                }
            }

            return matching.AsEnumerable().Reverse().ToList();
        }

        private void PrintDynProgAndOrReconstructionMatrix(double[,] dynProgMatrix, int[,] reconstructionMatrix)
        {
            TextWriter w = new StreamWriter(Path.Combine(Setting.Directory, "dyn_prog_matrix.tsv"));
            w.Write("\t");
            for (int j = 0; j < Context.Candidates.Length + 1; j++)
            {
                w.Write(j == 0 ? "-\t" : Context.Candidates[j - 1].Label + "\t");
            }
            w.WriteLine();
            for (int i = 0; i < Context.Templates.Length + 1; i++)
            {
                if (reconstructionMatrix != null)
                {
                    w.Write("\t");
                    for (int j = 0; j < Context.Candidates.Length + 1; j++)
                    {
                        w.Write(reconstructionMatrix[i, j] == -1 ? ">" : reconstructionMatrix[i, j] == 1 ? "v" : "\\");
                        w.Write("\t");
                    }
                    w.WriteLine();
                }
                if (dynProgMatrix != null)
                {
                    w.Write(i == 0 ? "-\t" : Context.Templates[i - 1].Label + "\t");
                    for (int j = 0; j < Context.Candidates.Length + 1; j++)
                    {
                        w.Write(dynProgMatrix[i, j].ToString("0.0") + "\t");
                    }
                    w.WriteLine();
                }
            }
            w.Close();
        }
    }

}
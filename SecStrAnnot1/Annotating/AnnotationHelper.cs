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
    static class AnnotationHelper
    {
        public static bool CheckConnectivity(AnnotationContext context, List<(int, int)> matching)
        {
            return matching.All(m1 => matching.All(m2 => context.TemplateConnectivity[m1.Item1, m2.Item1] == context.CandidateConnectivity[m1.Item2, m2.Item2]));
        }

        /**Return true if sheet IDs are consistent.*/
        public static bool CheckSheetIDConsistency(AnnotationContext context, List<(int, int)> matching, out List<(int, int)> sheetIDMatching)
        {
            if (context.Templates.Any(sse => sse.IsSheet && sse.SheetId == null)
                || context.Candidates.Any(sse => sse.IsSheet && sse.SheetId == null))
            {
                //skip checking (some strand miss sheet ID)
                Lib.WriteWarning("Cannot check sheet ID consistency because some strand do not have sheet ID.");
                sheetIDMatching = null;
                return true;
            }
            else
            {
                Dictionary<int?, List<int?>> dictSheetIdByTemplate = new Dictionary<int?, List<int?>>();
                Dictionary<int?, List<int?>> dictSheetIdByCandidate = new Dictionary<int?, List<int?>>();
                foreach (var m in matching)
                {
                    Sse template = context.Templates[m.Item1];
                    Sse candidate = context.Candidates[m.Item2];
                    if (template.IsSheet && candidate.IsSheet)
                    {
                        dictSheetIdByTemplate.MultidictionaryAdd(template.SheetId, candidate.SheetId);
                        dictSheetIdByCandidate.MultidictionaryAdd(candidate.SheetId, template.SheetId);
                    }
                }
                if (dictSheetIdByTemplate.Values.Any(l => l.Distinct().Count() > 1)
                    || dictSheetIdByCandidate.Values.Any(l => l.Distinct().Count() > 1))
                {
                    // sheet ID mismatch found
                    Lib.WriteWarning("Some beta-strands caused a sheet ID mismatch!");
                    foreach (var kv in dictSheetIdByTemplate.Where(kv => kv.Value.Distinct().Count() > 1))
                    {
                        Lib.WriteWarning("Template sheet {0} matched to query sheets {1}. ", kv.Key, kv.Value.Distinct().EnumerateWithCommas());
                    }
                    foreach (var kv in dictSheetIdByCandidate.Where(kv => kv.Value.Distinct().Count() > 1))
                    {
                        Lib.WriteWarning("Query sheet {0} matched to template sheets {1}. ", kv.Key, kv.Value.Distinct().EnumerateWithCommas());
                    }
                    sheetIDMatching = null;
                    return false;
                }
                else
                {
                    sheetIDMatching = dictSheetIdByTemplate.Select(kv => (kv.Key.Value, kv.Value[0].Value)).ToList();
                    return true;
                }
            }
        }

        public static void PrintMetricMatrix<T>(AnnotationContext context, T[,] matrix, String filename)
        {
            PrintMatrix(context.Templates.Select(sse => sse.Label).ToArray(),
                context.Candidates.Select(sse => sse.Label).ToArray(),
                context.Templates.Select(sse => (sse.EndPoint - sse.StartPoint).Size).ToArray(),
                matrix, filename);
        }

        public static void PrintTemplateMatrix<T>(AnnotationContext context, T[,] matrix, String filename)
        {
            String[] names = context.Templates.Select(sse => sse.Label).ToArray();
            PrintMatrix(names, names, null as object[], matrix, filename);
        }

        public static void PrintCandidateMatrix<T>(AnnotationContext context, T[,] matrix, String filename)
        {
            String[] names = context.Candidates.Select(sse => sse.Label).ToArray();
            PrintMatrix(names, names, null as object[], matrix, filename);
        }

        public static void PrintCrossMatrix<T>(AnnotationContext context, T[,] matrix, String filename)
        {
            String[] tNames = context.Templates.Select(sse => sse.Label).ToArray();
            String[] cNames = context.Candidates.Select(sse => sse.Label).ToArray();
            PrintMatrix(tNames, cNames, null as object[], matrix, filename);
        }

        public static void PrintMatrix<T, T2>(String[] rowNames, String[] colNames, T2[] rowExtras, T[,] matrix, String filename)
        {
            TextWriter w = new StreamWriter(Path.Combine(Setting.Directory, filename));
            int m = rowNames.Length;
            int n = colNames.Length;
            bool extras = rowExtras != null;
            if (matrix.GetLength(0) != m || matrix.GetLength(1) != n || (extras && rowExtras.Length != m))
                throw new Exception("Wrong dimensions.");
            w.Write("\t");
            if (extras)
                w.Write("<extra_info>\t");
            for (int j = 0; j < n; j++)
            {
                w.Write(colNames[j] + "\t");
            }
            w.WriteLine();
            for (int i = 0; i < m; i++)
            {
                w.Write(rowNames[i] + "\t");
                if (extras)
                    w.Write(rowExtras[i].ToString() + "\t");
                for (int j = 0; j < n; j++)
                {
                    w.Write(matrix[i, j].ToString(/*"0.0"*/) + "\t");
                }
                w.WriteLine();
            }
            w.WriteLine();
            w.WriteLine();
            w.Close();
        }
    }
}


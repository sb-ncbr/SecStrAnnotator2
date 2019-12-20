using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.SSEs;
using protein.Libraries;

namespace protein.Annotating
{
    class MOMAnnotator : IAnnotator
    {
        public MOMAnnotationContext MContext { get; private set; }
        public bool SoftOrderConsistency { get; private set; }

        public MOMAnnotator(AnnotationContext context, bool softOrderConsistency)
        {
            MContext = new MOMAnnotationContext(context);
            SoftOrderConsistency = softOrderConsistency;
        }

        public List<(int, int)> GetMatching()
        {
            MContext.Context.ValidateOrdering();
            int m = MContext.TemplateFST.Length;
            int n = MContext.CandidateFST.Length;

            Func<int, int, double> score = (i, j) =>
                  MContext.Context.SkipTemplatePenalty(MContext.Context.Templates[i]) + MContext.Context.SkipCandidatePenalty(MContext.Context.Candidates[j]) - MContext.Context.MetricToMin(MContext.Context.Templates[i], MContext.Context.Candidates[j]);
            Func<int, int, bool> typeMatch = (i, j) =>
                  MContext.Context.TypeMatching(MContext.Context.Templates[i].Type, MContext.Context.Candidates[j].Type);
            Func<int, bool> templateIsHelix = (i) =>
                 MContext.Context.Templates[MContext.TemplateFST[i].Item1].IsHelix;
            Func<int, int, bool> guideDis = (i, j) =>
                  MContext.Context.GuideDiscriminator == null || MContext.Context.GuideDiscriminator[i, j];

            double LADDER_SCORE_SCALE = 0.5;

            double[,] scores = new double[m, n].Fill((i, j) =>
               MContext.TemplateFST[i].Item3 == MContext.CandidateFST[j].Item3
                              && typeMatch(MContext.TemplateFST[i].Item1, MContext.CandidateFST[j].Item1)
                              && typeMatch(MContext.TemplateFST[i].Item2, MContext.CandidateFST[j].Item2)
                              && guideDis(MContext.TemplateFST[i].Item1, MContext.CandidateFST[j].Item1)
                              && guideDis(MContext.TemplateFST[i].Item2, MContext.CandidateFST[j].Item2) ?
               (templateIsHelix(i) ?
                   score(MContext.TemplateFST[i].Item1, MContext.CandidateFST[j].Item1)
                   : LADDER_SCORE_SCALE * (score(MContext.TemplateFST[i].Item1, MContext.CandidateFST[j].Item1) + score(MContext.TemplateFST[i].Item2, MContext.CandidateFST[j].Item2))
               )
               : 0);
            // double[,] scoreScales = new double[m, n].Fill ((i, j) => 
            // MContext.TemplateFST [i].Item3 == MContext.CandidateFST [j].Item3
            //                && typeMatch (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1)
            //                && typeMatch (MContext.TemplateFST [i].Item2, MContext.CandidateFST [j].Item2)
            //                && guideDis (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1)
            //                && guideDis (MContext.TemplateFST [i].Item2, MContext.CandidateFST [j].Item2) ?
            // (templateIsHelix(i) ? 
            // 	1
            // 	: LADDER_SCORE_SCALE
            // )
            // : 1);
            if (Lib.DoWriteDebug)
            {
                AnnotationHelper.PrintMatrix(
                    MContext.TemplateFST.Select(t => MContext.Context.Templates[t.Item1].Label + "-" + MContext.Context.Templates[t.Item2].Label).ToArray(),
                    MContext.CandidateFST.Select(t => MContext.Context.Candidates[t.Item1].Label + "-" + MContext.Context.Candidates[t.Item2].Label).ToArray(),
                    null as object[],
                    scores,
                    "score_matrix.tsv");
            }
            /*PrintMetricMatrix (Context, scores);
            if (Context.CandidateExclusivity != null)
                PrintCandidateMatrix (Context, Context.CandidateExclusivity, "c_exclusivity.tsv");
            if (Context.CandidateConnectivity != null)
                PrintCandidateMatrix (Context, Context.CandidateConnectivity, "c_connectivity.tsv");*/

            DateTime stamp = DateTime.Now;
            Lib.WriteLineDebug("{3}.GetAnnotation(): initialized - {0} vs. {1} vertices ({2})", m, n, stamp, this.GetType().Name);

            List<(int, int)> bestMOM = LibAnnotation.MaxWeightMixedOrderedMatching(m, n,
                MContext.TemplateFST.Select(t => t.Item1).ToArray(),
                MContext.TemplateFST.Select(t => t.Item2).ToArray(),
                MContext.CandidateFST.Select(t => t.Item1).ToArray(),
                MContext.CandidateFST.Select(t => t.Item2).ToArray(),
                scores, SoftOrderConsistency).OrderBy(x => x).ToList();

            Lib.WriteLineDebug("MOMAnnotator.GetAnnotation(): found matching ({0})", DateTime.Now);
            Lib.WriteLineDebug("MOMAnnotator.GetAnnotation(): time: {0}", DateTime.Now - stamp);

            List<(int, int)> bestMatching = new List<(int, int)>();
            foreach ((int ti, int qi) in bestMOM)
            {
                (int, int, int) tLadder = MContext.TemplateFST[ti];
                (int, int, int) qLadder = MContext.CandidateFST[qi];
                if (tLadder.Item3 == 0)
                {//helix
                    bestMatching.Add((tLadder.Item1, qLadder.Item1));
                }
                else
                {//ladder
                    bestMatching.Add((tLadder.Item1, qLadder.Item1));
                    bestMatching.Add((tLadder.Item2, qLadder.Item2));
                }
            }

            Lib.WriteLineDebug("Template FST: {0}", MContext.TemplateFST.Select(t => "(" + MContext.Context.Templates[t.Item1].Label + "," + MContext.Context.Templates[t.Item2].Label + ")").EnumerateWithCommas());
            Lib.WriteLineDebug("Candidate FST: {0}", MContext.CandidateFST.Select(t => "(" + MContext.Context.Candidates[t.Item1].Label + "," + MContext.Context.Candidates[t.Item2].Label + ")").EnumerateWithCommas());
            Lib.WriteLineDebug("MOM matching: {0}", bestMOM.EnumerateWithCommas());
            Lib.WriteLineDebug("Matching: {0}", bestMatching.Select(t => "(" + MContext.Context.Templates[t.Item1].Label + "," + MContext.Context.Candidates[t.Item2].Label + ")").EnumerateWithCommas());

            return bestMatching;
        }

    }

}
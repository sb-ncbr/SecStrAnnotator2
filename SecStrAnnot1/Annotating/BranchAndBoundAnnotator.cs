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
    class BranchAndBoundAnnotator : IAnnotator
    {
        public AnnotationContext Context { get; private set; }
        public BranchAndBoundAnnotator(AnnotationContext context)
        {
            context.ValidateOrdering();
            context.ValidateBetaGraph();
            this.Context = context;
        }
        public static BranchAndBoundAnnotator New(AnnotationContext context)
        {
            return new BranchAndBoundAnnotator(context);
        }

        public List<(int, int)> GetMatching()
        {
            int m = Context.Templates.Length;
            int n = Context.Candidates.Length;
            double[,] scores = new double[m, n].Fill((i, j) =>
                Context.TypeMatching(Context.Templates[i].Type, Context.Candidates[j].Type) && (Context.GuideDiscriminator == null || Context.GuideDiscriminator[i, j]) ?
                Context.SkipTemplatePenalty(Context.Templates[i]) + Context.SkipCandidatePenalty(Context.Candidates[j]) - Context.MetricToMin(Context.Templates[i], Context.Candidates[j])
                : 0);
            if (Lib.DoWriteDebug)
            {
                AnnotationHelper.PrintMetricMatrix(Context, scores, "score_matrix.tsv");
            }
            if (Context.CandidateExclusivity != null)
                AnnotationHelper.PrintCandidateMatrix(Context, Context.CandidateExclusivity, "c_exclusivity.tsv");
            if (Context.CandidateConnectivity != null)
                AnnotationHelper.PrintCandidateMatrix(Context, Context.CandidateConnectivity, "c_connectivity.tsv");

            /*Lib.WriteLineDebug ("Candidate exclusivity:");
            for (int i = 0; i < n; i++) {
                for (int j = i+1; j < n; j++) {
                    if (Context.CandidateExclusivity [i, j])
                        Lib.WriteLineDebug ("Ex: {0}-{1} / {2}-{3}", Context.Candidates [i].Start, Context.Candidates [i].End, Context.Candidates [j].Start, Context.Candidates [j].End);
                }
            }*/

            DateTime stamp = DateTime.Now;
            Lib.WriteLineDebug("BranchAndBoundAnnotator.GetAnnotation(): initialized - {0} vs. {1} vertices ({2})", m, n, stamp);

            List<(int, int)> bestMatching = LibAnnotation.MaxWeightOrderedMatching(m, n, Context.TemplateConnectivity, Context.CandidateConnectivity, scores, Context.TemplateExclusivity, Context.CandidateExclusivity);

            Lib.WriteLineDebug("BranchAndBoundAnnotator.GetAnnotation(): found matching ({0})", DateTime.Now);
            Lib.WriteLineDebug("BranchAndBoundAnnotator.GetAnnotation(): time: {0}", DateTime.Now - stamp);


            return bestMatching;
        }
    }

}
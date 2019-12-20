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
    class CombinedAnnotator : IAnnotator
    {
        public AnnotationContext Context { get; private set; }
        public CombinedAnnotator(AnnotationContext context)
        {
            context.ValidateOrdering();
            context.ValidateBetaGraph();
            this.Context = context;
        }
        public static CombinedAnnotator New(AnnotationContext context)
        {
            return new CombinedAnnotator(context);
        }

        public List<(int, int)> GetMatching()
        {
            DynProgAnnotator dynProgAn = new DynProgAnnotator(Context);
            List<(int, int)> dynProgMatching = dynProgAn.GetMatching();
            if (AnnotationHelper.CheckConnectivity(Context, dynProgMatching))
            {
                return dynProgMatching;
            }
            else
            {
                Lib.WriteLineDebug("CombinedAnnotator.GetMatching(): Dynamic programming caused connectivity issues. Trying dynamic programming with guide strand matching from branch-and-bound.");
                Lib.Shuffler tShuffler;
                Lib.Shuffler qShuffler;
                IEnumerable<SSEInSpace> filteredTemplates = Context.Templates.WhereAndGetShuffler(sse => sse.IsSheet, out tShuffler);
                IEnumerable<SSEInSpace> filteredCandidates = Context.Candidates.WhereAndGetShuffler(sse => sse.IsSheet, out qShuffler);
                int[,] filteredTemplateConnectivity = tShuffler.ShuffleColumns(tShuffler.ShuffleRows(Context.TemplateConnectivity));
                int[,] filteredCandidateConnectivity = qShuffler.ShuffleColumns(qShuffler.ShuffleRows(Context.CandidateConnectivity));

                AnnotationContext bbContext = new AnnotationContext(Context.MetricToMin, Context.TypeMatching,
                                                  Context.SkipTemplatePenalty, Context.SkipCandidatePenalty,
                                                  filteredTemplates, filteredCandidates);
                bbContext.TemplateConnectivity = filteredTemplateConnectivity;
                bbContext.CandidateConnectivity = filteredCandidateConnectivity;

                BranchAndBoundAnnotator bb = new BranchAndBoundAnnotator(bbContext);
                List<(int, int)> guideMatching = bb.GetMatching().Select(m => (tShuffler.OldIndex(m.Item1), qShuffler.OldIndex(m.Item2))).ToList();
                bool[,] guideDiscriminator = new bool[Context.Templates.Length, Context.Candidates.Length]
                    .Fill((i, j) => !Context.Templates[i].IsSheet && !Context.Candidates[j].IsSheet);
                foreach (var m in guideMatching)
                {
                    guideDiscriminator[m.Item1, m.Item2] = true;
                }
                AnnotationContext dpContext = new AnnotationContext(Context.MetricToMin, Context.TypeMatching,
                                                  Context.SkipTemplatePenalty, Context.SkipCandidatePenalty,
                                                  Context.Templates, Context.Candidates);
                dpContext.GuideDiscriminator = guideDiscriminator;
                DynProgAnnotator dp = new DynProgAnnotator(dpContext);
                return dp.GetMatching();
            }
        }
    }
}
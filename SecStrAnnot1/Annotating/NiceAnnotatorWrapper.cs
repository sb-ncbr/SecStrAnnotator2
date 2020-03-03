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
    class NiceAnnotatorWrapper
    {
        const double SUSPICIOUSNESS_THRESHOLD = 1.0;
        const int SUSPICIOUS_JOINING_GAP = 5;

        private IAnnotator inner;
        public AnnotationContext Context { get; private set; }

        public Lib.Shuffler TemplateMapping { get; private set; }
        public Lib.Shuffler CandidateMapping { get; private set; }
        private List<(int, int)> rememberedMatching;
        private SseInSpace[] rememberedAnnotatedCandidates;
        public bool DoCheckSheetIDConsistency { get; set; }
        public bool DoRenameSheetIDs { get; set; }

        public NiceAnnotatorWrapper(AnnotationContext context, Func<AnnotationContext, IAnnotator> innerAnnotatorConstructor)
        {
            Context = context;
            Lib.Shuffler m;
            Context.Templates.OrderAndGetShuffler(out m);
            TemplateMapping = m;
            Context.Candidates.OrderAndGetShuffler(out m);
            CandidateMapping = m;
            DoCheckSheetIDConsistency = true;
            DoRenameSheetIDs = true;

            AnnotationContext innerContext = Context.Copy();
            innerContext.ApplyShufflers(TemplateMapping, CandidateMapping);
            inner = innerAnnotatorConstructor(innerContext);

            Lib.WriteLineDebug("Templates: {0}", Context.Templates.Select(s => s.Start).EnumerateWithCommas());
            Lib.WriteLineDebug("Template mapping: {0}", TemplateMapping);
            Lib.WriteLineDebug("Mapped templates: {0}", innerContext.Templates.Select(s => s.Start).EnumerateWithCommas());

        }

        public List<(int, int)> GetMatching()
        {
            if (inner == null)
                throw new Exception("Inner annotator has not been initialized.");
            if (rememberedMatching == null)
                rememberedMatching = inner.GetMatching().Select(m => (TemplateMapping.OldIndex(m.Item1), CandidateMapping.OldIndex(m.Item2))).ToList();
            return rememberedMatching;
        }

        public virtual IEnumerable<SseInSpace> GetAnnotatedCandidates()
        {
            if (rememberedAnnotatedCandidates == null)
            {
                List<(int, int)> matching = GetMatching();
                /*Lib.WriteLineDebug ("Matching: {0}",matching.EnumerateWithCommas ());
                Lib.WriteLineDebug ("temps: {0}", matching.Select (t => t.Item1).OrderBy (x => x).EnumerateWithCommas ());
                Lib.WriteLineDebug ("cands: {0}", matching.Select (t => t.Item2).OrderBy (x => x).EnumerateWithCommas ());*/
                if (matching.Select(m => m.Item1).Distinct().Count() == matching.Count && matching.Select(m => m.Item2).Distinct().Count() == matching.Count)
                {
                    // each vertex is matched at most once
                    Lib.Shuffler annotShuffler = Lib.Shuffler.FromMatching(matching);
                    rememberedAnnotatedCandidates = annotShuffler.ShuffleBack(Context.Candidates, () => SseInSpace.NewNotFound(null))
                        .Select((sse, i) => sse.RelabeledCopy(Context.Templates[i].Label, Context.Templates[i].Color, Context.Templates[i].Rainbow)).ToArray();
                }
                else
                {
                    // some vertices are matched more than once
                    matching = matching.OrderBy(m => m.Item1).ToList();
                    Dictionary<int, List<int>> multiMatching = new Dictionary<int, List<int>>();
                    foreach (var m in matching)
                        multiMatching.MultidictionaryAdd(m.Item1, m.Item2);
                    rememberedAnnotatedCandidates = new SseInSpace[Context.Templates.Length];
                    for (int i = 0; i < Context.Templates.Length; i++)
                    {
                        if (!multiMatching.ContainsKey(i))
                        {
                            rememberedAnnotatedCandidates[i] = SseInSpace.NewNotFound(Context.Templates[i].Label);
                        }
                        else
                        {
                            SseInSpace[] all = multiMatching[i].Select(j => Context.Candidates[j]).OrderBy(x => x).Distinct().ToArray();
                            if (all.Length == 1)
                            {
                                rememberedAnnotatedCandidates[i] = all[0].RelabeledCopy(Context.Templates[i].Label, Context.Templates[i].Color, Context.Templates[i].Rainbow);
                            }
                            else if (all.Length > 1)
                            {
                                if (all.Select(sse => sse.ChainID).Distinct().Count() > 1)
                                    throw new Exception("Strands from different chains have been annotated as parts of one strand.");
                                for (int j = 0; j < all.Length - 1; j++)
                                {
                                    if (all[j + 1].Start - all[j].End - 1 >= SUSPICIOUS_JOINING_GAP)
                                    {
                                        Lib.WriteWarning("Suspicious joining in {0}. Gap between joined SSEs = {1}", Context.Templates[i].Label, all[j + 1].Start - all[j].End - 1);
                                    }
                                }
                                rememberedAnnotatedCandidates[i] = SseInSpace.Join(all).RelabeledCopy(Context.Templates[i].Label);
                            }
                            else
                                throw new Exception("This should never happen!");
                        }
                    }
                }
                if (DoCheckSheetIDConsistency || DoRenameSheetIDs)
                {
                    List<(int, int)> sheetIDMatching;
                    AnnotationHelper.CheckSheetIDConsistency(Context, matching, out sheetIDMatching);
                    if (DoRenameSheetIDs && sheetIDMatching != null)
                    {
                        RenameSheetIDs(rememberedAnnotatedCandidates, sheetIDMatching);
                    }
                }
            }
            return rememberedAnnotatedCandidates;
        }

        public virtual List<(int, int, int)> GetAnnotatedConnectivity(List<(int, int, int)> connectivity)
        {
            Lib.Shuffler annotShuffler = Lib.Shuffler.FromMatching(GetMatching()).Inverted();
            List<(int, int, int)> annotConnectivity = annotShuffler.UpdateIndices(connectivity).ToList();
            return annotConnectivity;
        }

        public IEnumerable<T> SelectFromAnnotated<T>(Func<SseInSpace, SseInSpace, T> selector)
        {
            return GetAnnotatedCandidates()
                .Select((sse, i) => selector(Context.Templates[i], sse));
        }
        public IEnumerable<double> GetMetricList()
        {
            return SelectFromAnnotated((t, q) => (q != null && !q.IsNotFound()) ? Context.MetricToMin(t, q) : 0);
        }

        public void RenameSheetIDs(IEnumerable<SseInSpace> sses, List<(int, int)> sheetIDMatching)
        {
            Lib.Shuffler shuffler = Lib.Shuffler.FromMatching(sheetIDMatching);
            foreach (SseInSpace sse in sses)
            {
                if (sse.SheetId != null)
                {
                    sse.SheetId = shuffler.OldIndex((int)sse.SheetId);
                }
            }
        }

        public List<double> GetSuspiciousnessList()
        {
            List<(int, int)> matching = GetMatching();
            double[] result = new double[Context.Templates.Length];
            foreach (var m in matching)
            {
                SseInSpace template = Context.Templates[m.Item1];
                SseInSpace candidate = Context.Candidates[m.Item2];
                double otherMinInRow = Context.Candidates
                    .Where((c, i) => i != m.Item2 && Context.TypeMatching(template.Type, c.Type))
                    .Select(c => (double?)Context.MetricToMin(template, c))
                    .Min() ?? Double.PositiveInfinity;
                double otherMinInColumn = Context.Templates
                    .Where((t, i) => i != m.Item1 && Context.TypeMatching(t.Type, candidate.Type))
                    .Select(t => (double?)Context.MetricToMin(t, candidate))
                    .Min() ?? Double.PositiveInfinity;
                double otherMin = Math.Min(otherMinInRow, otherMinInColumn);
                double annotVal = Context.MetricToMin(template, candidate);
                double suspiciousness = annotVal / (annotVal + otherMin);
                result[m.Item1] = suspiciousness;
                if (suspiciousness >= SUSPICIOUSNESS_THRESHOLD)
                {
                    Lib.WriteWarning("Possibly ambiguous annotation around \"{0}\" (annotated / (alternative+annotated) = {1}).", template.Label, suspiciousness.ToString("0.00"));
                }
            }
            return result.ToList();
        }
    }
}
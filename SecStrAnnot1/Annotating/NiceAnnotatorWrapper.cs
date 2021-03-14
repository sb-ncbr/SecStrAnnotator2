using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.Annotating {
    class NiceAnnotatorWrapper {
        const double SUSPICIOUSNESS_THRESHOLD = 1.0;
        const int SUSPICIOUS_JOINING_GAP = 5;

        private IAnnotator inner;
        public AnnotationContext Context { get; private set; }

        public Lib.Shuffler TemplateMapping { get; private set; }
        public Lib.Shuffler CandidateMapping { get; private set; }
        private List<(int, int)> matching;
        private SseInSpace[] annotatedCandidates;
        // private SseInSpace[] unannotatedCandidates;
        public readonly bool DoCheckSheetIDConsistency;
        public readonly bool DoRenameSheetIDs;
        public readonly bool IncludeUnannotatedSSEs;

        public NiceAnnotatorWrapper(AnnotationContext context, Func<AnnotationContext, IAnnotator> innerAnnotatorConstructor, bool includeUnannotatedSSEs = false) {
            Context = context;
            Lib.Shuffler m;
            Context.Templates.OrderAndGetShuffler(out m);
            TemplateMapping = m;
            Context.Candidates.OrderAndGetShuffler(out m);
            CandidateMapping = m;
            DoCheckSheetIDConsistency = true;
            DoRenameSheetIDs = true;
            IncludeUnannotatedSSEs = includeUnannotatedSSEs;

            AnnotationContext innerContext = Context.Copy();
            innerContext.ApplyShufflers(TemplateMapping, CandidateMapping);
            inner = innerAnnotatorConstructor(innerContext);

            // Lib.WriteLineDebug("Templates: {0}", Context.Templates.Select(s => s.Start).EnumerateWithCommas());
            // Lib.WriteLineDebug("Template mapping: {0}", TemplateMapping);
            // Lib.WriteLineDebug("Mapped templates: {0}", innerContext.Templates.Select(s => s.Start).EnumerateWithCommas());

        }

        public List<(int, int)> GetMatching() {
            if (inner == null)
                throw new Exception("Inner annotator has not been initialized.");
            if (matching == null)
                matching = inner.GetMatching().Select(m => (TemplateMapping.OldIndex(m.Item1), CandidateMapping.OldIndex(m.Item2))).ToList();
            return matching;
        }

        protected virtual SseInSpace[] CalculateAnnotatedCandidates() {
            SseInSpace[] result;
            List<(int, int)> matching = GetMatching();
            /*Lib.WriteLineDebug ("Matching: {0}",matching.EnumerateWithCommas ());
            Lib.WriteLineDebug ("temps: {0}", matching.Select (t => t.Item1).OrderBy (x => x).EnumerateWithCommas ());
            Lib.WriteLineDebug ("cands: {0}", matching.Select (t => t.Item2).OrderBy (x => x).EnumerateWithCommas ());*/
            if (matching.Select(m => m.Item1).Distinct().Count() == matching.Count && matching.Select(m => m.Item2).Distinct().Count() == matching.Count) {
                // each vertex is matched at most once
                Lib.Shuffler annotShuffler = Lib.Shuffler.FromMatching(matching);
                result = annotShuffler.ShuffleBack(Context.Candidates, () => SseInSpace.NewNotFound(null))
                    .Select((sse, i) => sse.RelabeledCopy(Context.Templates[i].Label, Context.Templates[i].Color, Context.Templates[i].Rainbow)).ToArray();
            } else {
                // some vertices are matched more than once
                matching = matching.OrderBy(m => m.Item1).ToList();
                Dictionary<int, List<int>> multiMatching = new Dictionary<int, List<int>>();
                foreach (var m in matching)
                    multiMatching.MultidictionaryAdd(m.Item1, m.Item2);
                result = new SseInSpace[Context.Templates.Length];
                for (int i = 0; i < Context.Templates.Length; i++) {
                    if (!multiMatching.ContainsKey(i)) {
                        result[i] = SseInSpace.NewNotFound(Context.Templates[i].Label);
                    } else {
                        SseInSpace[] all = multiMatching[i].Select(j => Context.Candidates[j]).OrderBy(x => x).Distinct().ToArray();
                        if (all.Length == 1) {
                            result[i] = all[0].RelabeledCopy(Context.Templates[i].Label, Context.Templates[i].Color, Context.Templates[i].Rainbow);
                        } else if (all.Length > 1) {
                            if (all.Select(sse => sse.ChainID).Distinct().Count() > 1)
                                throw new Exception("Strands from different chains have been annotated as parts of one strand.");
                            for (int j = 0; j < all.Length - 1; j++) {
                                if (all[j + 1].Start - all[j].End - 1 >= SUSPICIOUS_JOINING_GAP) {
                                    Lib.WriteWarning("Suspicious joining in {0}. Gap between joined SSEs = {1}", Context.Templates[i].Label, all[j + 1].Start - all[j].End - 1);
                                }
                            }
                            result[i] = SseInSpace.Join(all).RelabeledCopy(Context.Templates[i].Label, Context.Templates[i].Color, Context.Templates[i].Rainbow);
                        } else
                            throw new Exception("This should never happen!");
                    }
                }
            }
            if (IncludeUnannotatedSSEs) {
                var matchedCandidates = (from m in matching select m.Item2).ToHashSet();
                var unannotatedCandidates = (
                    from i in Enumerable.Range(0, Context.Candidates.Length) 
                    where !matchedCandidates.Contains(i) 
                    select Context.Candidates[i]);
                unannotatedCandidates = unannotatedCandidates.Select(sse => sse.RelabeledCopy("_" + sse.Label));
                result = result.Concat(unannotatedCandidates).ToArray();
            }
            if (DoCheckSheetIDConsistency || DoRenameSheetIDs) {
                List<(int, int)> sheetIDMatching;
                AnnotationHelper.CheckSheetIDConsistency(Context, matching, out sheetIDMatching);
                if (DoRenameSheetIDs && sheetIDMatching != null) {
                    RenameSheetIDs(result, sheetIDMatching);
                }
            }
            return result;
        }

        public SseInSpace[] GetAnnotatedCandidates() {
            if (annotatedCandidates == null) {
                annotatedCandidates = CalculateAnnotatedCandidates();
            }
            return annotatedCandidates;
        }

        public virtual List<(int, int, int)> GetAnnotatedConnectivity(List<(int, int, int)> connectivity) {
            var indexMapping = new Dictionary<int,int>();
            foreach(var match in matching){
                indexMapping[match.Item2] = match.Item1;
            }
            if (IncludeUnannotatedSSEs){
                int nextAnnotIndex = Context.Templates.Length;
                for (int i = 0; i < Context.Candidates.Length; i++){
                    if (!indexMapping.ContainsKey(i)){
                        indexMapping[i] = nextAnnotIndex++;
                    }
                }
            }
            var annotConnectivity = new List<(int, int, int)>();
            foreach((int first, int second, int dir) in connectivity){
                if (indexMapping.ContainsKey(first) && indexMapping.ContainsKey(second)){
                    annotConnectivity.Add((indexMapping[first], indexMapping[second], dir));
                }
            } 
            // connectivity.Select(conn => (indexMapping[conn.Item1], indexMapping[conn.Item2], conn.Item3)).ToList();
            return annotConnectivity;
        }

        public T[] SelectFromAnnotated<T>(Func<SseInSpace, SseInSpace, T> selector) {
            var annotCands = GetAnnotatedCandidates();
            T[] result = new T[annotCands.Length];
            for (int i = 0; i < annotCands.Length; i++)
            {
                var temp = (i < Context.Templates.Length) ? Context.Templates[i] : null;
                var cand = annotCands[i];
                result[i] = selector(temp, cand);
            }
            return result;
            // return GetAnnotatedCandidates().Select((sse, i) => (i < Context.Templates.Length) ? selector(Context.Templates[i], sse)).ToArray();
        }
        public double[] GetMetricList() {
            return SelectFromAnnotated((t, q) => (t != null && q != null && !q.IsNotFound()) ? Context.MetricToMin(t, q) : 0.0);
        }

        public void RenameSheetIDs(IEnumerable<SseInSpace> sses, List<(int, int)> sheetIDMatching) {
            Lib.Shuffler shuffler = Lib.Shuffler.FromMatching(sheetIDMatching);
            IEnumerable<int?> sheetIds = sses.Select(sse => sse.SheetId).Distinct();
            List<int> nonmatchedSheetIds = sheetIds.Where(shid => shid.HasValue && !shuffler.HasOldIndex(shid.Value)).Select(shid => shid.Value).ToList();
            nonmatchedSheetIds.Sort();
            int lastOldIndex = shuffler.MaxOldIndex;
            var nonmatchedMapping = new Dictionary<int,int>();
            foreach (int shid in nonmatchedSheetIds){
                nonmatchedMapping[shid] = ++lastOldIndex;
            }
            foreach (SseInSpace sse in sses) {
                if (sse.SheetId != null) {
                    int shid = sse.SheetId.Value;
                    sse.SheetId = shuffler.HasOldIndex(shid) ? shuffler.OldIndex(shid) : nonmatchedMapping[shid];
                }
            }
        }

        public double[] GetSuspiciousnessList() {
            List<(int, int)> matching = GetMatching();
            double[] result = new double[GetAnnotatedCandidates().Length];
            foreach (var m in matching) {
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
                if (suspiciousness >= SUSPICIOUSNESS_THRESHOLD) {
                    Lib.WriteWarning("Possibly ambiguous annotation around \"{0}\" (annotated / (alternative+annotated) = {1}).", template.Label, suspiciousness.ToString("0.00"));
                }
            }
            return result.ToArray();
        }
    }
}
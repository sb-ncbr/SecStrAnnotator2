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
    class AnnotationContext
    {
        //General settings
        public Func<SSEInSpace, SSEInSpace, double> MetricToMin { get; private set; }
        public Func<SSEType, SSEType, bool> TypeMatching { get; private set; }
        public Func<SSEInSpace, double> SkipTemplatePenalty { get; private set; }
        public Func<SSEInSpace, double> SkipCandidatePenalty { get; private set; }

        //Templates and candidates
        public SSEInSpace[] Templates { get; private set; }
        public SSEInSpace[] Candidates { get; private set; }
        public bool[,] GuideDiscriminator { get; set; }

        //Beta-graph related context
        public int[,] TemplateConnectivity { get; set; } //1 = parallel beta-ladder, -1 = antiparallel beta-ladder, 0 = no ladder
        public int[,] CandidateConnectivity { get; set; } //1 = parallel beta-ladder, -1 = antiparallel beta-ladder, 0 = no ladder
        public bool[,] TemplateExclusivity { get; set; }
        public bool[,] CandidateExclusivity { get; set; }

        //Constructors
        public AnnotationContext(
            Func<SSEInSpace, SSEInSpace, double> metricToMin,
            Func<SSEType, SSEType, bool> typeMatching,
            Func<SSEInSpace, double> skipTemplatePenalty,
            Func<SSEInSpace, double> skipCandidatePenalty,
            IEnumerable<SSEInSpace> templates,
            IEnumerable<SSEInSpace> candidates
        )
        {
            this.MetricToMin = metricToMin;
            this.TypeMatching = typeMatching;
            this.SkipTemplatePenalty = skipTemplatePenalty;
            this.SkipCandidatePenalty = skipCandidatePenalty;

            this.Templates = templates.ToArray();
            this.Candidates = candidates.ToArray();
        }

        public void InitializeTemplateConnectivity(IEnumerable<(int, int, int)> templateConnections)
        {
            TemplateConnectivity = new int[Templates.Length, Templates.Length];
            foreach (var connection in templateConnections)
            {
                TemplateConnectivity[connection.Item1, connection.Item2] = connection.Item3;
                TemplateConnectivity[connection.Item2, connection.Item1] = connection.Item3;
            }
        }

        public void InitializeCandidateConnectivity(IEnumerable<(int, int, int)> candidateConnections)
        {
            CandidateConnectivity = new int[Candidates.Length, Candidates.Length];
            foreach (var connection in candidateConnections)
            {
                CandidateConnectivity[connection.Item1, connection.Item2] = connection.Item3;
                CandidateConnectivity[connection.Item2, connection.Item1] = connection.Item3;
            }
        }

        public AnnotationContext Copy()
        {
            AnnotationContext copy = new AnnotationContext(this.MetricToMin, this.TypeMatching,
                this.SkipTemplatePenalty, this.SkipCandidatePenalty,
                this.Templates, this.Candidates);
            copy.GuideDiscriminator = this.GuideDiscriminator;
            copy.TemplateConnectivity = this.TemplateConnectivity;
            copy.CandidateConnectivity = this.CandidateConnectivity;
            copy.TemplateExclusivity = this.TemplateExclusivity;
            copy.CandidateExclusivity = this.CandidateExclusivity;
            return copy;
        }

        public void ApplyShufflers(Lib.Shuffler templateShuffler, Lib.Shuffler candidateShuffler)
        {
            Templates = templateShuffler.Shuffle(Templates).ToArray();
            Candidates = candidateShuffler.Shuffle(Candidates).ToArray();
            if (GuideDiscriminator != null)
            {
                GuideDiscriminator = templateShuffler.ShuffleRows(candidateShuffler.ShuffleColumns(GuideDiscriminator));
            }

            if (TemplateConnectivity != null)
            {
                TemplateConnectivity = templateShuffler.ShuffleRowsAndColumns(TemplateConnectivity);
            }

            if (CandidateConnectivity != null)
            {
                CandidateConnectivity = candidateShuffler.ShuffleRowsAndColumns(CandidateConnectivity);
            }

            if (TemplateExclusivity != null)
            {
                TemplateExclusivity = templateShuffler.ShuffleRowsAndColumns(TemplateExclusivity);
            }

            if (CandidateExclusivity != null)
            {
                CandidateExclusivity = candidateShuffler.ShuffleRowsAndColumns(CandidateExclusivity);
            }
        }

        public void ValidateOrdering()
        {
            if (Templates.Select(x => x.ChainID).Distinct().Count() > 1)
                Lib.WriteWarning("{0}: passed template SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod().Name);
            if (Candidates.Select(x => x.ChainID).Distinct().Count() > 1)
                Lib.WriteWarning("{0}: passed candidate SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod().Name);

            if (!Templates.IsSorted())
                throw new Exception("Template SSEs are not ordered.");
            if (!Candidates.IsSorted())
                throw new Exception("Candidate SSEs are not ordered.");
        }

        public void ValidateBetaGraph()
        {
            if (TemplateConnectivity == null)
                throw new Exception("Template connectivity is null.");
            if (CandidateConnectivity == null)
                throw new Exception("Candidate connectivity is null.");
        }

        public AnnotationContext Ordered()
        {
            Lib.Shuffler tShuffler;
            Lib.Shuffler cShuffler;
            return Ordered(out tShuffler, out cShuffler);
        }
        public AnnotationContext Ordered(out Lib.Shuffler tShuffler, out Lib.Shuffler cShuffler)
        {
            AnnotationContext result = new AnnotationContext(this.MetricToMin, this.TypeMatching,
                this.SkipTemplatePenalty, this.SkipCandidatePenalty,
                this.Templates.OrderAndGetShuffler(out tShuffler), this.Candidates.OrderAndGetShuffler(out cShuffler));
            result.GuideDiscriminator = tShuffler.ShuffleRows(cShuffler.ShuffleColumns(this.GuideDiscriminator));
            result.TemplateConnectivity = tShuffler.ShuffleRowsAndColumns(this.TemplateConnectivity);
            result.CandidateConnectivity = cShuffler.ShuffleRowsAndColumns(this.CandidateConnectivity);
            result.TemplateExclusivity = tShuffler.ShuffleRowsAndColumns(this.TemplateExclusivity);
            result.CandidateExclusivity = cShuffler.ShuffleRowsAndColumns(this.CandidateExclusivity);
            return result;
        }

        public AnnotationContext Softened(int maxGap)
        {
            //soft matching - joining candidates 
            ValidateOrdering();
            ValidateBetaGraph();

            List<SSEInSpace> jointCandidates = new List<SSEInSpace>();
            List<(int, int, int)> jointConnections = new List<(int, int, int)>();
            List<(int, int)> conflicts = new List<(int, int)>();
            List<(int, int)> jointGuideAllowedMatches = new List<(int, int)>();
            for (int i = 0; i < Candidates.Length - 1; i++)
            {
                SSEInSpace s1 = Candidates[i];
                SSEInSpace s2 = Candidates[i + 1];
                if (s1.ChainID == s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap)
                {
                    int index = Candidates.Length + jointCandidates.Count;
                    // SSEInSpace newSSE = new SSEInSpace (new SSE (s1.Label + "+" + s2.Label, s1.ChainID, s1.Start, s2.End, 'E', null), s1.StartVector, s2.EndVector);
                    // newSSE.AddComment ("Created by joining " + s1.Label + " and " + s2.Label + " with gap " + (s2.Start - s1.End - 1) + ".");
                    // newSSE.AddNestedSSE (s1);
                    // newSSE.AddNestedSSE (s2);
                    SSEInSpace newSSE = SSEInSpace.Join(s1, s2, null);
                    jointCandidates.Add(newSSE);
                    if (GuideDiscriminator != null)
                    {
                        for (int j = 0; j < Templates.Length; j++)
                        {
                            if (GuideDiscriminator[j, i] || GuideDiscriminator[j, i + 1])
                                jointGuideAllowedMatches.Add((j, index));
                        }
                    }
                    for (int j = 0; j < Candidates.Length; j++)
                    {
                        if (CandidateConnectivity[i, j] != 0 && CandidateConnectivity[i + 1, j] != 0 && CandidateConnectivity[i, j] != CandidateConnectivity[i, j])
                            throw new Exception("Joint strand would have both parallel and antiparallel ladder to another strand.");
                        if (CandidateConnectivity[i, j] == 1 || CandidateConnectivity[i + 1, j] == 1)
                            jointConnections.Add((index, j, 1));
                        if (CandidateConnectivity[i, j] == -1 || CandidateConnectivity[i + 1, j] == -1)
                            jointConnections.Add((index, j, -1));
                    }
                    conflicts.Add((i, index));
                    conflicts.Add((i + 1, index));
                }
            }
            SSEInSpace[] newCandidates = Candidates.Concat(jointCandidates).ToArray();
            bool[,] newGuideDiscriminator;
            if (GuideDiscriminator != null)
            {
                newGuideDiscriminator = GuideDiscriminator?.Resized2D(Templates.Length, newCandidates.Length);
                foreach (var allowed in jointGuideAllowedMatches)
                {
                    newGuideDiscriminator[allowed.Item1, allowed.Item2] = true;
                }
            }
            else
            {
                newGuideDiscriminator = new bool[Templates.Length, newCandidates.Length].Fill(true);
            }
            int[,] newCandidateConnectivity = CandidateConnectivity.Resized2D(newCandidates.Length, newCandidates.Length);
            foreach (var conn in jointConnections)
            {
                newCandidateConnectivity[conn.Item1, conn.Item2] = conn.Item3;
                newCandidateConnectivity[conn.Item2, conn.Item1] = conn.Item3;
            }
            bool[,] newCandidateExclusivity = CandidateExclusivity?.Resized2D(newCandidates.Length, newCandidates.Length) ?? new bool[newCandidates.Length, newCandidates.Length];
            foreach (var conf in conflicts)
            {
                newCandidateExclusivity[conf.Item1, conf.Item2] = true;
                newCandidateExclusivity[conf.Item2, conf.Item1] = true;
            }

            //Extra limitations to soft matching to disable useless joining
            int[] tDegreePar = Enumerable.Range(0, Templates.Length).Select(i => Enumerable.Range(0, Templates.Length).Select(j => TemplateConnectivity[i, j]).Count(x => x == 1)).ToArray();
            int[] tDegreeAnti = Enumerable.Range(0, Templates.Length).Select(i => Enumerable.Range(0, Templates.Length).Select(j => TemplateConnectivity[i, j]).Count(x => x == -1)).ToArray();
            int[] jcDegreePar = Enumerable.Range(Candidates.Length, jointCandidates.Count).Select(i => Enumerable.Range(0, Candidates.Length).Select(j => newCandidateConnectivity[i, j]).Count(x => x == 1)).ToArray();
            int[] jcDegreeAnti = Enumerable.Range(Candidates.Length, jointCandidates.Count).Select(i => Enumerable.Range(0, Candidates.Length).Select(j => newCandidateConnectivity[i, j]).Count(x => x == -1)).ToArray();
            for (int i = 0; i < Templates.Length; i++)
            {
                for (int j = 0; j < jointCandidates.Count; j++)
                {
                    if (jcDegreePar[j] > tDegreePar[i] || jcDegreeAnti[j] > tDegreeAnti[i])
                        newGuideDiscriminator[i, Candidates.Length + j] = false;
                }
            }

            Lib.Shuffler shuffler;
            newCandidates = newCandidates.OrderAndGetShuffler(out shuffler).ToArray();
            AnnotationContext result = new AnnotationContext(MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, Templates, newCandidates);
            result.GuideDiscriminator = newGuideDiscriminator != null ? shuffler.ShuffleColumns(newGuideDiscriminator) : null;
            result.TemplateConnectivity = TemplateConnectivity;
            result.CandidateConnectivity = shuffler.ShuffleRowsAndColumns(newCandidateConnectivity);
            result.TemplateExclusivity = TemplateExclusivity;
            result.CandidateExclusivity = shuffler.ShuffleRowsAndColumns(newCandidateExclusivity);

            return result;
        }

        public AnnotationContext Softened_New(int maxGap)
        {
            //soft matching - joining candidates 
            ValidateOrdering();
            ValidateBetaGraph();

            List<(int, int)> joint = new List<(int, int)>();
            //Allow joining only 2 strands
            /*for (int i = 0; i < Candidates.Length - 1; i++) {
                int j = i + 1;
                SSEInSpace s1 = Candidates [i];
                SSEInSpace s2 = Candidates [j];
                if (s1.ChainID == s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap && CandidateConnectivity[i,j]==0) {
                    joint.Add ((i, j));
                }
            }*/
            //Allow joining more than 2 strands
            for (int i = 0; i < Candidates.Length - 1; i++)
            {
                for (int j = i + 1; j < Candidates.Length; j++)
                {
                    SSEInSpace s1 = Candidates[i];
                    SSEInSpace s2 = Candidates[j];
                    if (s1.ChainID == s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap && Enumerable.Range(i, j - i - 1).All(x => CandidateConnectivity[x, j] == 0))
                        joint.Add((i, j));
                    else
                        break;
                }
            }

            return SoftenedTemplatesOrCandidates(joint, 1).Ordered();
        }

        public AnnotationContext WithAlternativeTemplates(List<(String, int, int)> alternatives)
        {
            //soft matching - joining templates 
            ValidateBetaGraph();

            List<String> labels = alternatives.Select(a => a.Item1).ToList();
            List<(int, int)> ranges = alternatives.Select(a => (a.Item2, a.Item3)).ToList();

            int mOrig = this.Templates.Length;

            AnnotationContext result = this.SoftenedTemplatesOrCandidates(ranges, 0);
            for (int i = 0; i < labels.Count; i++)
            {
                result.Templates[mOrig + i] = result.Templates[mOrig + i].RelabeledCopy(labels[i]);
            }
            return result.Ordered();
        }

        /*, whichGraph==0 --> soften templates, whichGraph==1 --> soften candidates; no ordering apllied! */
        private AnnotationContext SoftenedTemplatesOrCandidates(List<(int, int)> joinedRanges, int whichGraph)
        {
            bool STRICT = false;
            SSEInSpace[] gSSEs; //unchanged
            SSEInSpace[] hSSEs; //to be softened
            bool[,] ghGuideDiscriminator;
            int[,] gConnectivity;
            int[,] hConnectivity;
            bool[,] gExclusivity;
            bool[,] hExclusivity;

            if (whichGraph == 0)
            {
                gSSEs = Candidates;
                hSSEs = Templates;
                ghGuideDiscriminator = GuideDiscriminator?.Transposed();
                gConnectivity = CandidateConnectivity;
                hConnectivity = TemplateConnectivity;
                gExclusivity = CandidateExclusivity;
                hExclusivity = TemplateExclusivity;
            }
            else if (whichGraph == 1)
            {
                gSSEs = Templates;
                hSSEs = Candidates;
                ghGuideDiscriminator = GuideDiscriminator;
                gConnectivity = TemplateConnectivity;
                hConnectivity = CandidateConnectivity;
                gExclusivity = TemplateExclusivity;
                hExclusivity = CandidateExclusivity;
            }
            else
            {
                throw new Exception("Invalid argument whichGraph.");
            }

            int m = gSSEs.Length;
            int nOrig = hSSEs.Length;
            int nNew = nOrig + joinedRanges.Count;
            hSSEs = hSSEs.Resized(nNew);
            ghGuideDiscriminator = ghGuideDiscriminator?.Resized2D(m, nNew) ?? new bool[m, nNew].Fill(true);
            hConnectivity = hConnectivity.Resized2D(nNew, nNew);
            hExclusivity = hExclusivity?.Resized2D(nNew, nNew) ?? new bool[nNew, nNew];

            List<(int, int, int)> jointConnections = new List<(int, int, int)>();
            List<(int, int)> conflicts = new List<(int, int)>();
            for (int k = 0; k < joinedRanges.Count; k++)
            {
                int i = joinedRanges[k].Item1;
                int j = joinedRanges[k].Item2;
                int index = nOrig + k;
                int[] partIndices = Enumerable.Range(i, j - i + 1).ToArray();
                IEnumerable<SSEInSpace> parts = hSSEs.Skip(i).Take(j - i + 1);
                SSEInSpace newSSE = new SSEInSpace(new SSE(parts.Select(p => p.Label).EnumerateWithSeparators("_"), parts.First().ChainID, parts.First().Start, parts.Last().End, SSEType.SHEET_TYPE, null),
                    parts.First().StartPoint, parts.Last().EndPoint);
                newSSE.AddComment("Created by joining " + parts.EnumerateWithSeparators(" and ") + ".");
                foreach (var p in parts)
                {
                    newSSE.AddNestedSSE(p);
                }
                hSSEs[index] = newSSE;
                if (ghGuideDiscriminator != null)
                {
                    for (int t = 0; t < m; t++)
                    {
                        ghGuideDiscriminator[t, index] = partIndices.Any(q => ghGuideDiscriminator[t, q]);
                    }
                }
                //Connectivity original-joint
                for (int q = 0; q < nOrig; q++)
                {
                    bool parConn = partIndices.Any(p => hConnectivity[q, p] == 1);
                    bool antiConn = partIndices.Any(p => hConnectivity[q, p] == -1);
                    if (parConn && antiConn)
                        throw new Exception("Joint strand would have both parallel and antiparallel ladder to another strand.");
                    if (parConn || antiConn)
                        jointConnections.Add((index, q, parConn ? 1 : -1));
                }
                conflicts.AddRange(partIndices.Select(p => (p, index)));
            }

            //Connectivity joint-joint, exclusivity joint-joint
            for (int k = 0; k < joinedRanges.Count; k++)
            {
                for (int l = k + 1; l < joinedRanges.Count; l++)
                {
                    int kIndex = nOrig + k;
                    int lIndex = nOrig + l;
                    IEnumerable<int> kPartIndices = Enumerable.Range(joinedRanges[k].Item1, joinedRanges[k].Item2 - joinedRanges[k].Item1 + 1);
                    IEnumerable<int> lPartIndices = Enumerable.Range(joinedRanges[l].Item1, joinedRanges[l].Item2 - joinedRanges[l].Item1 + 1);
                    bool parConn = kPartIndices.Any(kp => lPartIndices.Any(lp => hConnectivity[kp, lp] == 1));
                    bool antiConn = kPartIndices.Any(kp => lPartIndices.Any(lp => hConnectivity[kp, lp] == -1));
                    if (parConn && antiConn)
                        throw new Exception("Joint strand would have both parallel and antiparallel ladder to another strand.");
                    if (parConn || antiConn)
                        jointConnections.Add((kIndex, lIndex, parConn ? 1 : -1));
                    if (joinedRanges[k].Item2 >= joinedRanges[l].Item1)
                        conflicts.Add((kIndex, lIndex));
                }
            }

            foreach (var conn in jointConnections)
            {
                hConnectivity[conn.Item1, conn.Item2] = conn.Item3;
                hConnectivity[conn.Item2, conn.Item1] = conn.Item3;
            }
            foreach (var conf in conflicts)
            {
                hExclusivity[conf.Item1, conf.Item2] = true;
                hExclusivity[conf.Item2, conf.Item1] = true;
            }

            //Extra limitations to soft matching to disable useless joining
            if (STRICT)
            {
                int[] tDegreePar = Enumerable.Range(0, m).Select(i => Enumerable.Range(0, m).Select(j => gConnectivity[i, j]).Count(x => x == 1)).ToArray();
                int[] tDegreeAnti = Enumerable.Range(0, m).Select(i => Enumerable.Range(0, m).Select(j => gConnectivity[i, j]).Count(x => x == -1)).ToArray();
                int[] jcDegreePar = Enumerable.Range(nOrig, joinedRanges.Count).Select(i => Enumerable.Range(0, nOrig).Select(j => hConnectivity[i, j]).Count(x => x == 1)).ToArray();
                int[] jcDegreeAnti = Enumerable.Range(nOrig, joinedRanges.Count).Select(i => Enumerable.Range(0, nOrig).Select(j => hConnectivity[i, j]).Count(x => x == -1)).ToArray();
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < joinedRanges.Count; j++)
                    {
                        if (jcDegreePar[j] > tDegreePar[i] || jcDegreeAnti[j] > tDegreeAnti[i])
                            ghGuideDiscriminator[i, nOrig + j] = false;
                    }
                }
            }

            /*Lib.Shuffler shuffler;
            hSSEs = hSSEs.OrderAndGetShuffler (out shuffler).ToArray ();
            ghGuideDiscriminator = shuffler.ShuffleColumns (ghGuideDiscriminator);
            hConnectivity = shuffler.ShuffleRowsAndColumns (hConnectivity);
            hExclusivity = shuffler.ShuffleRowsAndColumns (hExclusivity);*/

            if (whichGraph == 0)
            {
                AnnotationContext result = new AnnotationContext(MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, hSSEs, gSSEs);
                result.GuideDiscriminator = ghGuideDiscriminator.Transposed();
                result.TemplateConnectivity = hConnectivity;
                result.CandidateConnectivity = gConnectivity;
                result.TemplateExclusivity = hExclusivity;
                result.CandidateExclusivity = gExclusivity;
                return result;
            }
            else if (whichGraph == 1)
            {
                AnnotationContext result = new AnnotationContext(MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, gSSEs, hSSEs);
                result.GuideDiscriminator = ghGuideDiscriminator;
                result.TemplateConnectivity = gConnectivity;
                result.CandidateConnectivity = hConnectivity;
                result.TemplateExclusivity = gExclusivity;
                result.CandidateExclusivity = hExclusivity;
                return result;
            }
            else
            {
                throw new Exception("Invalid argument whichGraph.");
            }
        }
    }
}
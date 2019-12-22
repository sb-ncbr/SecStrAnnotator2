using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using SecStrAnnotator2.Utils;
using Cif.Components;
using Cif.Tables;
using protein.Libraries;
using protein.Sses;
using protein.HydrogenAdding;
using protein.SecStrAssigning.Helpers;

namespace protein.SecStrAssigning
{
    public class HBondSecStrAssigner2 : ISecStrAssigner
    {
        private const int NONEXISTING_RESIDUE_INDEX = -1;

        private Model model;
        private Residue[] residues;
        private int nResidues;
        private HashSet<(int donor, int acceptor)> hBonds;
        private Dictionary<int, List<int>> acceptorsOf;
        public bool DetectHelices;
        public bool DetectSheets;

        public HBondSecStrAssigner2(Protein protein, double hBondEnergyCutoff)
        {
            this.DetectSheets = true;
            this.DetectHelices = true;
            MyStopwatch watch = new MyStopwatch();
            var hydrogenAdder = new DsspLikeAmideHydrogenAdder();
            protein = hydrogenAdder.AddHydrogens(protein); // removing residues without C-alpha is probably not needed here
            this.model = protein.Model;
            this.residues = protein.GetResidues().ToArray();
            this.nResidues = model.Residues.Count;
            MyStopwatch watchHBF = new MyStopwatch();
            var hBondFinder = new BoxingHBondFinder(this.residues, hBondEnergyCutoff);
            watchHBF.Stop("new BoxingHbondFinder");
            this.hBonds = new HashSet<(int donor, int acceptor)>();
            this.acceptorsOf = new Dictionary<int, List<int>>();
            for (int don = 0; don < nResidues; don++)
            {
                var acceptors = hBondFinder.FindHAcceptors(don);
                acceptorsOf[don] = acceptors;
                foreach (int acc in acceptors)
                {
                    hBonds.Add((don, acc));
                }
            }
            watchHBF.Stop("Found H-bonds");
            watch.Stop("new HBondSecStrAssigner2()");
        }

        public String GetDescription()
        {
            return "hydrogen-bond-based method (similar to DSSP)";
        }

        public SecStrAssignment GetSecStrAssignment()
        {
            SecStrAssignment assignment = new SecStrAssignment(new List<Sse>());
            assignment.HBonds = HBondsAsResidueTuples(this.hBonds);
            if (DetectSheets)
            {
                (List<Sse> strands, var connectivity) = GetStrands();
                assignment.SSEs.AddRange(strands);
                assignment.Connectivity = connectivity;
            }
            if (DetectHelices)
            {
                List<Sse> helices = GetHelices();
                assignment.SSEs.AddRange(helices);
            }
            return assignment;
        }

        private List<Sse> GetHelices()
        {
            List<Ladder> ladders = new List<Ladder>();
            for (int don = 0; don < nResidues; don++)
            {
                foreach (int acc in acceptorsOf[don])
                {
                    MicrostrandType type;
                    if (don == Move(acc, +3))
                        type = MicrostrandType.REGULAR_G_HELIX;
                    else if (don == Move(acc, +4))
                        type = MicrostrandType.REGULAR_H_HELIX;
                    else if (don == Move(acc, +5))
                        type = MicrostrandType.REGULAR_I_HELIX;
                    else
                        continue;
                    Hbond current = new Hbond(don, acc, false);
                    Hbond previous = PreviousInHelix(current);
                    Hbond next = NextInHelix(current);
                    bool startsHelix = Exists(next) && !Exists(previous);
                    if (startsHelix)
                    {
                        Ladder newLadder = ExtendHelix(current, next, type);
                        ladders.Add(newLadder);
                    }
                }
            }
            List<Sse> helices = new List<Sse>(capacity: ladders.Count);
            for (int i = 0; i < ladders.Count; i++)
            {
                Ladder ladder = ladders[i];
                int startZeta = ladder.firstHbond.Zeta1;
                int endZeta = ladder.lastHbond.Zeta0;
                (int startResIdx, int endResIdx) = ZRangeToResidueIndexRange(startZeta, endZeta, HBondSSAConstants.HELICES_BY_ALPHA);
                (string chainId, int startSeqId, int endSeqId) = ChainStartEnd(startResIdx, endResIdx);
                SseType type =
                    ladder.type0 == MicrostrandType.REGULAR_H_HELIX ? SseType.HELIX_H_TYPE
                    : ladder.type0 == MicrostrandType.REGULAR_G_HELIX ? SseType.HELIX_G_TYPE
                    : SseType.HELIX_I_TYPE;
                Sse sse = new Sse(null/*$"{type}{i}"*/, chainId, startSeqId, endSeqId, type, null);
                helices.Add(sse);
            }
            return helices;
        }

        private (List<Sse> strands, List<(int strand0, int strand1, int type)> connectivity) GetStrands()
        {
            MyStopwatch watch = new MyStopwatch();
            BetaGraph betaGraph = BuildBetaGraph();
            watch.Stop("BuildBetaGraph()");

            int nStrands = betaGraph.macrostrands.Count;
            List<Sse> strands = new List<Sse>(capacity: nStrands);
            var antiparallelBulgeTypes = new HashSet<MicrostrandType>{
                MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_SHORT_SIDE,
                MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_LONG_SIDE,
                MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_EQUAL_SIDE };
            var parallelBulgeTypes = new HashSet<MicrostrandType>{
                MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_SHORT_SIDE,
                MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_LONG_SIDE,
                MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_EQUAL_SIDE };
            for (int i = 0; i < nStrands; i++)
            {
                var strand = betaGraph.macrostrands[i];
                (int startResIdx, int endResIdx) = ZRangeToResidueIndexRange(strand.startZeta, strand.endZeta, HBondSSAConstants.STRANDS_BY_ALPHA);
                (string chainId, int startSeqId, int endSeqId) = ChainStartEnd(startResIdx, endResIdx);
                SseType type = (startResIdx == endResIdx) ? SseType.ISOLATED_BETA_BRIDGE_TYPE : SseType.SHEET_TYPE;
                int sheetId = HBondSSAConstants.SHEET_NUMBERING_FROM + strand.sheet;
                Sse sse = new Sse(null/*$"{type}{i}"*/, chainId, startSeqId, endSeqId, type, sheetId);
                if (Lib.DoWriteDebug) {
                    foreach (int iMis in strand.microstrands) {
                        var mis = betaGraph.microstrands[iMis];
                        int iMil = mis.microladder;
                        var mil = betaGraph.microladders[iMil];
                        if (parallelBulgeTypes.Contains(mil.type0)){
                            string bulgeCode = ParallelBulgeCode(mil, invert: mis.side == 1);
                            (startResIdx, endResIdx) = ZRangeToResidueIndexRange(mis.startZeta, mis.endZeta, HBondSSAConstants.BULGES_BY_ALPHA);
                            (chainId, startSeqId, endSeqId) = ChainStartEnd(startResIdx, endResIdx);
                            type = SseType.BULGE_PARALLEL_UNSPECIFIED;
                            Sse bulge = new Sse(null, chainId, startSeqId, endSeqId, type, -iMil);
                            bulge.AddComment(bulgeCode);
                            sse.AddNestedSSE(bulge);
                        }
                        if (antiparallelBulgeTypes.Contains(mil.type0)){
                            string bulgeCode = AntiparallelBulgeCode(mil, invert: mis.side == 1);
                            (startResIdx, endResIdx) = ZRangeToResidueIndexRange(mis.startZeta, mis.endZeta, HBondSSAConstants.BULGES_BY_ALPHA);
                            (chainId, startSeqId, endSeqId) = ChainStartEnd(startResIdx, endResIdx);
                            type = SseType.BULGE_ANTIPARALLEL_UNSPECIFIED;
                            Sse bulge = new Sse(null, chainId, startSeqId, endSeqId, type, -iMil);
                            bulge.AddComment(bulgeCode);
                            sse.AddNestedSSE(bulge);
                        }
                    }
                }
                strands.Add(sse);
            }
            var connectivity = betaGraph.macroladders.Select(t => (t.macrostrand0, t.macrostrand1, t.type == MacroladderType.PARALLEL ? 1 : -1)).ToList();
            return (strands, connectivity);
        }

        private BetaGraph BuildBetaGraph()
        {

            MyStopwatch watch = new MyStopwatch();
            List<Ladder> ladders = FindLadders();
            watch.Stop("FindLadders()");
            Ladder[] parallelLadders = ladders.Where(l => l.type0 == MicrostrandType.REGULAR_PARALEL).ToArray();
            Ladder[] antiparallelLadders = ladders.Where(l => l.type0 == MicrostrandType.REGULAR_ANTIPARALEL).ToArray();
            List<Ladder> parallelLaddersAndBulges = IncludeParallelBulges(parallelLadders);
            List<Ladder> antiparallelLaddersAndBulges = IncludeAntiparallelBulges(antiparallelLadders);
            List<Ladder> microladders = parallelLaddersAndBulges.Concat(antiparallelLaddersAndBulges).ToList();

            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug("Microladders:");
                foreach (Ladder ladder in microladders)
                {
                    Lib.WriteLineDebug(RepresentLadder(ladder));
                }
            }

            var microstrands = new List<(int startZeta, int endZeta, int ladderIndex, int side)>(capacity: 2 * (parallelLaddersAndBulges.Count + antiparallelLaddersAndBulges.Count));
            int nParallel = parallelLaddersAndBulges.Count;
            int nAntiparallel = antiparallelLaddersAndBulges.Count;
            int nMicroladders = nParallel + nAntiparallel;
            // Adding microstrands of parallel ladders and bulges
            for (int i = 0; i < nParallel; i++)
            {
                Ladder ladder = microladders[i];
                microstrands.Add((ladder.firstHbond.Zeta0, ladder.lastHbond.Zeta0, i, 0));
                microstrands.Add((ladder.firstHbond.Zeta1, ladder.lastHbond.Zeta1, i, 1));
            }
            // Adding microstrands of antiparallel ladders and bulges
            for (int i = nParallel; i < nParallel + nAntiparallel; i++)
            {
                Ladder ladder = microladders[i];
                microstrands.Add((ladder.firstHbond.Zeta0, ladder.lastHbond.Zeta0, i, 0));
                microstrands.Add((ladder.lastHbond.Zeta1, ladder.firstHbond.Zeta1, i, 1));
            }
            int nMicrostrands = microstrands.Count;
            var sortedMicrostrandIndices = Enumerable.Range(0, nMicrostrands).OrderBy(i => microstrands[i]).ToArray();
            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug("Microstrands:");
                for (int i = 0; i < nMicrostrands; i++)
                {
                    var ms = microstrands[sortedMicrostrandIndices[i]];
                    Lib.WriteLineDebug(ms.ToString());
                }
            }

            // Building macrostrands
            var macrostrands = new List<(int startZeta, int endZeta, List<int> microstrands)>(capacity: nMicrostrands);
            for (int i = 0; i < nMicrostrands; i++)
            {
                bool starting = (i == 0);
                int microstrandIndex = sortedMicrostrandIndices[i];
                var microstrand = microstrands[microstrandIndex];
                var lastMacrostrand = starting ? (0, 0, null) : macrostrands[macrostrands.Count - 1];
                if (!starting && DoOverlap(lastMacrostrand.startZeta, lastMacrostrand.endZeta, microstrand.startZeta, microstrand.endZeta))
                {
                    Lib.WriteLineDebug($"Joining macro {lastMacrostrand} and micro {microstrand}");
                    lastMacrostrand.microstrands.Add(microstrandIndex);
                    if (microstrand.endZeta > lastMacrostrand.endZeta)
                    {
                        macrostrands[macrostrands.Count - 1] = (lastMacrostrand.startZeta, microstrand.endZeta, lastMacrostrand.microstrands);
                    }
                }
                else
                {
                    macrostrands.Add((microstrand.startZeta, microstrand.endZeta, new List<int> { microstrandIndex }));
                }
            }
            int nMacrostrands = macrostrands.Count;
            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug("Macrostrands:");
                foreach (var ms in macrostrands)
                {
                    int startResIdx = ZToResidueIndex(ms.startZeta);
                    int endResIdx = ZToResidueIndex(ms.endZeta);
                    string chainId = model.Chains.Id[model.Residues.ChainIndex[startResIdx]];
                    int startSeqId = model.Residues.SeqNumber[startResIdx];
                    int endSeqId = model.Residues.SeqNumber[endResIdx];
                    Lib.WriteLineDebug($"{chainId} {startSeqId}-{endSeqId} {ms.ToString()} " + Lib.EnumerateWithCommas(ms.microstrands));
                }
            }
            int[] microstrand2macrostrand = new int[nMicrostrands];
            for (int i = 0; i < nMacrostrands; i++)
            {
                foreach (int microstrandIndex in macrostrands[i].microstrands)
                {
                    microstrand2macrostrand[microstrandIndex] = i;
                }
            }
            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug(Lib.EnumerateWithCommas(microstrand2macrostrand));
            }

            var macroladders = new List<(int macrostrand0, int macrostrand1, MacroladderType type)>();
            for (int i = 0; i < nMicroladders; i++)
            {
                bool parallel = i < nParallel;
                MacroladderType type = parallel ? MacroladderType.PARALLEL : MacroladderType.ANTIPARALLEL;
                int mis0 = 2 * i;
                int mis1 = 2 * i + 1;
                int mas0 = microstrand2macrostrand[mis0];
                int mas1 = microstrand2macrostrand[mis1];
                macroladders.Add((mas0, mas1, type));
            }
            macroladders = macroladders.Distinct().OrderBy(mal => mal).ToList();
            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug("Macroladders:");
                Lib.WriteLineDebug(Lib.EnumerateWithSeparators(macroladders, "\n"));
            }
            var edges = macroladders.Select(mal => (mal.macrostrand0, mal.macrostrand1)).ToList();
            MyStopwatch watchComp = new MyStopwatch();
            var sheets = ConnectedComponents(edges);
            int nSheets = sheets.Count;
            watchComp.Stop("ConnecteComponents()");
            var macrostrand2sheet = new int[nMacrostrands];
            for (int iSheet = 0; iSheet < nSheets; iSheet++)
            {
                foreach (int iMacrostrand in sheets[iSheet])
                {
                    macrostrand2sheet[iMacrostrand] = iSheet;
                }
            }
            if (Lib.DoWriteDebug)
            {
                Lib.WriteLineDebug("Sheets:");
                foreach (var sheet in sheets)
                {
                    Lib.WriteLineDebug(Lib.EnumerateWithCommas(sheet));
                }
                Lib.WriteLineDebug("Macrostrand2sheet:");
                Lib.WriteLineDebug(Lib.EnumerateWithCommas(macrostrand2sheet));
                Lib.WriteLineDebug("");
            }
            watch.Stop("rest of BuildBetaGraph()");

            BetaGraph bg = new BetaGraph();
            bg.microstrands = Enumerable.Zip(microstrands, microstrand2macrostrand, (m, m2m) => (m.startZeta, m.endZeta, m.ladderIndex, m.side, m2m)).ToList();
            bg.microladders = microladders;
            bg.nParallelMicroladders = nParallel;
            bg.macrostrands = Enumerable.Zip(macrostrands, macrostrand2sheet, (m, m2s) => (m.startZeta, m.endZeta, m.microstrands, m2s)).ToList();
            bg.macroladders = macroladders;
            bg.sheets = sheets;

            return bg;
        }

        private List<Ladder> FindLadders()
        {
            List<Ladder> ladders = new List<Ladder>();
            for (int don = 0; don < nResidues; don++)
            {
                foreach (int acc in acceptorsOf[don])
                {
                    if (acc <= don)
                    {
                        continue;
                        // this avoids: 1) recognizing H-bond i->i-2 as start of a ladder, 2) recognizing each ladder twice
                    }
                    // antiparallel
                    Hbond current = new Hbond(don, acc, false);
                    Hbond previous = PreviousInAntiparallel(current);
                    Hbond preprevious = PreviousInAntiparallel(previous);
                    Hbond next = NextInAntiparallel(current);
                    bool startsMotifA = Exists(next) && !Exists(previous);
                    bool startsMotifB = Exists(previous) && !Exists(preprevious);
                    if (startsMotifA)
                    {
                        Ladder newLadder = ExtendAntiparallelLadder(current, next);
                        ladders.Add(newLadder);
                    }
                    else if (startsMotifB)
                    {
                        Ladder newLadder = ExtendAntiparallelLadder(previous, current);
                        ladders.Add(newLadder);
                    }
                    // parallel
                    previous = PreviousInParallel(current);
                    preprevious = PreviousInParallel(previous);
                    next = NextInParallel(current);
                    bool startsMotifC = Exists(next) && !Exists(previous);
                    bool startsMotifD = Exists(previous) && !Exists(preprevious);
                    if (startsMotifC)
                    {
                        Ladder newLadder = ExtendParallelLadder(current, next);
                        ladders.Add(newLadder);
                    }
                    else if (startsMotifD)
                    {
                        Ladder newLadder = ExtendParallelLadder(previous, current);
                        ladders.Add(newLadder);
                    }
                }
            }
            return ladders;
        }

        private string RepresentLadder(Ladder ladder)
        {
            string type = ladder.type0 == MicrostrandType.REGULAR_PARALEL ? "P "
                : ladder.type0 == MicrostrandType.REGULAR_ANTIPARALEL ? " A"
                : ladder.type0.ToString();
            string chain0 = residues[ladder.firstHbond.r0].ChainId;
            string chain1 = residues[ladder.firstHbond.r1].ChainId;
            int start0 = residues[ladder.firstHbond.r0].SeqNumber;
            int end0 = residues[ladder.lastHbond.r0].SeqNumber;
            int start1 = residues[ladder.firstHbond.r1].SeqNumber;
            int end1 = residues[ladder.lastHbond.r1].SeqNumber;
            return $"{type}[{chain0} {start0}-{end0} : {chain1} {start1}-{end1}]";
        }

        private List<(Residue, Residue)> HBondsAsResidueTuples(IEnumerable<(int donor, int acceptor)> hbonds)
        {
            return hbonds.OrderBy(hb => hb).Select(hb => (residues[hb.donor], residues[hb.acceptor])).ToList();
        }
        private List<Ladder> IncludeParallelBulges(Ladder[] parallelLadders)
        {  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int n = parallelLadders.Length;
            List<Ladder> result = new List<Ladder>(capacity: n);
            if (n == 0)
            {
                return result;
            }
            result.Add(parallelLadders[0]);
            for (int i = 1; i < n; i++)
            {
                Ladder? potentialBulge = TryMakeParallelBulge(parallelLadders[i - 1], parallelLadders[i]);
                if (potentialBulge != null)
                {
                    result.Add(potentialBulge.Value);
                }
                result.Add(parallelLadders[i]);
            }
            return result;
        }
        private List<Ladder> IncludeAntiparallelBulges(Ladder[] antiparallelLadders)
        {  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int n = antiparallelLadders.Length;
            List<Ladder> result = new List<Ladder>(capacity: n);
            if (n == 0)
            {
                return result;
            }
            result.Add(antiparallelLadders[0]);
            for (int i = 1; i < n; i++)
            {
                Ladder? potentialBulge = TryMakeAntiparallelBulge(antiparallelLadders[i - 1], antiparallelLadders[i]);
                if (potentialBulge != null)
                {
                    result.Add(potentialBulge.Value);
                    Lib.WriteLineDebug("Found bulge");
                }
                result.Add(antiparallelLadders[i]);
            }
            return result;
        }

        private Ladder? TryMakeParallelBulge(Ladder firstLadder, Ladder secondLadder)
        {  // Ladders must have strand0 in lower residue index and must be sorted wrt strand0 residue index!
            int zetaShift0 = secondLadder.firstHbond.Zeta0 - firstLadder.lastHbond.Zeta0;
            int zetaShift1 = secondLadder.firstHbond.Zeta1 - firstLadder.lastHbond.Zeta1;
            var fragmentIndex = model.Residues.FragmentIndex;
            bool sameFragments = fragmentIndex[firstLadder.firstHbond.r0] == fragmentIndex[secondLadder.firstHbond.r0] && fragmentIndex[firstLadder.firstHbond.r1] == fragmentIndex[secondLadder.firstHbond.r1];
            if (!sameFragments || zetaShift1 < 0)
            {
                return null;
            }
            Ladder? result = null;
            if (zetaShift0 < zetaShift1 && zetaShift0 <= 4 && zetaShift1 <= 13)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_SHORT_SIDE);
            }
            else if (zetaShift0 > zetaShift1 && zetaShift0 <= 13 && zetaShift1 <= 4)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_LONG_SIDE);
            }
            else if (zetaShift0 == zetaShift1 && zetaShift0 <= 4)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_EQUAL_SIDE);
            }
            // TODO distinguish bulge type
            if (result.HasValue)
            {
                Lib.WriteLineDebug("BULGE: " + ParallelBulgeCode(result.Value));
                Lib.WriteLineDebug($"{zetaShift0} {zetaShift1}");
            }
            return result;
        }

        private Ladder? TryMakeAntiparallelBulge(Ladder firstLadder, Ladder secondLadder)
        {  // Ladders must have strand0 in lower residue index and must be sorted wrt strand0 residue index!
            int zetaShift0 = secondLadder.firstHbond.Zeta0 - firstLadder.lastHbond.Zeta0;
            int zetaShift1 = firstLadder.lastHbond.Zeta1 - secondLadder.firstHbond.Zeta1;
            var fragmentIndex = model.Residues.FragmentIndex;
            bool sameFragments = fragmentIndex[firstLadder.firstHbond.r0] == fragmentIndex[secondLadder.firstHbond.r0] && fragmentIndex[firstLadder.firstHbond.r1] == fragmentIndex[secondLadder.firstHbond.r1];
            if (!sameFragments || zetaShift1 < 0)
            {
                return null;
            }
            Ladder? result = null;
            if (zetaShift0 < zetaShift1 && zetaShift0 <= 4 && zetaShift1 <= 13)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_SHORT_SIDE);
            }
            else if (zetaShift0 > zetaShift1 && zetaShift0 <= 13 && zetaShift1 <= 4)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_LONG_SIDE);
            }
            else if (zetaShift0 == zetaShift1 && zetaShift0 <= 4)
            {
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_EQUAL_SIDE);
            }
            // TODO distinguish bulge type
            if (result.HasValue)
            {
                Lib.WriteLineDebug("BULGE: " + AntiparallelBulgeCode(result.Value));
                Lib.WriteLineDebug($"{zetaShift0} {zetaShift1}");
            }
            return result;
        }

        private List<List<int>> ConnectedComponents(List<(int, int)> edges)
        {
            var components = new List<List<int>>();
            var neighbourLists = new Dictionary<int, List<int>>(capacity: 2 * edges.Count);
            foreach ((int u, int v) in edges)
            {
                if (!neighbourLists.ContainsKey(u)) neighbourLists[u] = new List<int>();
                if (!neighbourLists.ContainsKey(v)) neighbourLists[v] = new List<int>();
                neighbourLists[u].Add(v);
                neighbourLists[v].Add(u);
            }
            var remainingVertices = new HashSet<int>(neighbourLists.Keys);
            while (remainingVertices.Count > 0)
            {
                int seed = remainingVertices.First();
                var outList = new List<int>();
                FindComponentDFS(seed, neighbourLists, remainingVertices, outList);
                outList.Sort();
                components.Add(outList);
            }
            components = components.OrderBy(comp => comp[0]).ToList();
            return components;
        }

        private void FindComponentDFS(int seed, Dictionary<int, List<int>> neighbourLists, HashSet<int> remainingVertices, List<int> outList)
        {
            remainingVertices.Remove(seed);
            outList.Add(seed);
            foreach (int neighbour in neighbourLists[seed])
            {
                if (remainingVertices.Contains(neighbour))
                {
                    FindComponentDFS(neighbour, neighbourLists, remainingVertices, outList);
                }
            }
        }

        private const int DONOR_Z_CHANGE = -1;
        private const int ACCEPTOR_Z_CHANGE = 1;

        private static int ZToResidueIndex(int z) => (z - DONOR_Z_CHANGE) / 3;
        private static bool ZIsDonor(int z) => z % 3 != ACCEPTOR_Z_CHANGE; // Doesn't take into account alpha-carbon (not donor nor acceptor)!
        private static bool ZIsAcceptor(int z) => z % 3 == ACCEPTOR_Z_CHANGE; // Doesn't take into account alpha-carbon (not donor nor acceptor)!

        private (int, int) ZRangeToResidueIndexRange(int startZeta, int endZeta, bool byAlpha)
        {
            if (byAlpha)
            {
                int startResIdx = ZIsAcceptor(startZeta) ? ZToResidueIndex(startZeta) + 1 : ZToResidueIndex(startZeta);
                int endResIdx = ZIsAcceptor(endZeta) ? ZToResidueIndex(endZeta) : ZToResidueIndex(endZeta) - 1;
                return (startResIdx, endResIdx);
            }
            else
            {
                return (ZToResidueIndex(startZeta), ZToResidueIndex(endZeta));
            }
        }

        private (string chainId, int startSeqId, int endSeqId) ChainStartEnd(int startResidueIndex, int endResidueIndex)
        {
            string chainId = model.Chains.Id[model.Residues.ChainIndex[startResidueIndex]];
            int startSeqId = model.Residues.SeqNumber[startResidueIndex];
            int endSeqId = model.Residues.SeqNumber[endResidueIndex];
            return (chainId, startSeqId, endSeqId);
        }

        private int Move(int residueIndex, int move)
        {
            int candidateIndex = residueIndex + move;
            if (candidateIndex >= 0 && candidateIndex < nResidues && model.Residues.FragmentIndex[residueIndex] == model.Residues.FragmentIndex[candidateIndex])
            {
                return candidateIndex;
            }
            else
            {
                return NONEXISTING_RESIDUE_INDEX;
            }
        }

        private struct Hbond
        {
            public readonly int r0;  // residue index
            public readonly int r1;  // residue index
            public readonly bool opposite;  // if r1 is donor

            public Hbond(int r0, int r1, bool opposite)
            {
                this.r0 = r0;
                this.r1 = r1;
                this.opposite = opposite;
            }
            public (int donor, int acceptor) DonorAcceptor => opposite ? (r1, r0) : (r0, r1);
            public int Zeta0 => 3 * r0 + (opposite ? ACCEPTOR_Z_CHANGE : DONOR_Z_CHANGE);
            public int Zeta1 => 3 * r1 + (opposite ? DONOR_Z_CHANGE : ACCEPTOR_Z_CHANGE);
            public Hbond Inverted => new Hbond(r1, r0, !opposite);
        }

        private struct Ladder
        {
            public readonly Hbond firstHbond;
            public readonly Hbond lastHbond;
            public readonly MicrostrandType type0;

            public Ladder(Hbond firstHbond, Hbond lastHbond, MicrostrandType type0)
            {
                this.firstHbond = firstHbond;
                this.lastHbond = lastHbond;
                this.type0 = type0;
            }
        }

        private struct BetaGraph
        {
            public List<(int startZeta, int endZeta, int microladder, int side, int macrostrand)> microstrands;
            public List<Ladder> microladders;
            public int nParallelMicroladders;  // First nParallelMicroladders microladders are parallel, the rest is antiparallel
            public List<(int startZeta, int endZeta, List<int> microstrands, int sheet)> macrostrands;
            public List<(int macrostrand0, int macrostrand1, MacroladderType type)> macroladders;
            public List<List<int>> sheets;
        }

        private Hbond NextInAntiparallel(Hbond hbond)
        {
            if (hbond.opposite)
            {
                return new Hbond(Move(hbond.r0, +2), Move(hbond.r1, -2), false);
            }
            else
            {
                return new Hbond(hbond.r0, hbond.r1, true);
            }
        }
        private Hbond PreviousInAntiparallel(Hbond hbond)
        {
            if (hbond.opposite)
            {
                return new Hbond(hbond.r0, hbond.r1, false);
            }
            else
            {
                return new Hbond(Move(hbond.r0, -2), Move(hbond.r1, +2), true);
            }
        }
        private Hbond NextInParallel(Hbond hbond)
        {
            if (hbond.opposite)
            {
                return new Hbond(Move(hbond.r0, +2), hbond.r1, false);
            }
            else
            {
                return new Hbond(hbond.r0, Move(hbond.r1, +2), true);
            }
        }
        private Hbond PreviousInParallel(Hbond hbond)
        {
            if (hbond.opposite)
            {
                return new Hbond(hbond.r0, Move(hbond.r1, -2), false);
            }
            else
            {
                return new Hbond(Move(hbond.r0, -2), hbond.r1, true);
            }
        }
        private Hbond NextInHelix(Hbond hbond)
        {
            return new Hbond(Move(hbond.r0, +1), Move(hbond.r1, +1), hbond.opposite);
        }
        private Hbond PreviousInHelix(Hbond hbond)
        {
            return new Hbond(Move(hbond.r0, -1), Move(hbond.r1, -1), hbond.opposite);
        }
        private bool Exists(Hbond hbond)
        {
            return this.hBonds.Contains(hbond.DonorAcceptor);
        }
        private Ladder ExtendAntiparallelLadder(Hbond first, Hbond last)
        {
            while (true)
            {
                Hbond next = NextInAntiparallel(last);
                if (Exists(next) && next.r0 < next.r1)
                { // The second condition prevents beta-hairpin from extending over the end
                    last = next;
                }
                else
                {
                    break;
                }
            }
            return new Ladder(first, last, MicrostrandType.REGULAR_ANTIPARALEL);
        }
        private Ladder ExtendParallelLadder(Hbond first, Hbond last)
        {
            while (true)
            {
                Hbond next = NextInParallel(last);
                if (Exists(next))
                    last = next;
                else
                    break;
            }
            return new Ladder(first, last, MicrostrandType.REGULAR_PARALEL);
        }
        private Ladder ExtendHelix(Hbond first, Hbond last, MicrostrandType type)
        {
            while (true)
            {
                Hbond next = NextInHelix(last);
                if (Exists(next))
                    last = next;
                else
                    break;
            }
            return new Ladder(first, last, type);
        }

        private bool DoOverlap(int startZeta1, int endZeta1, int startZeta2, int endZeta2)
        {
            return startZeta2 <= endZeta1 - HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING
                && startZeta1 <= endZeta2 - HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING;
        }


        private string ParallelBulgeCode(Ladder bulge, bool invert = false){
            Hbond firstBond = bulge.firstHbond;
            Hbond lastBond = bulge.lastHbond;
            string side0 = BulgeSideCode(firstBond, lastBond);
            string side1 = BulgeSideCode(firstBond.Inverted, lastBond.Inverted);
            if (invert)
                return $"P_{side1}_{side0}_" ;
            else
                return $"P_{side0}_{side1}_";
        }

        private string AntiparallelBulgeCode(Ladder bulge, bool invert = false){
            Hbond firstBond = bulge.firstHbond;
            Hbond lastBond = bulge.lastHbond;
            string side0 = BulgeSideCode(firstBond, lastBond);
            string side1 = BulgeSideCode(lastBond.Inverted, firstBond.Inverted);
            if (invert)
                return $"A_{side1}_{side0}_" ;
            else
                return $"A_{side0}_{side1}_";
        }

        private string BulgeSideCode(Hbond firstHbond, Hbond lastHbond)
        {
            int riStart = firstHbond.r0;
            int riEnd = lastHbond.r0;
            if (firstHbond.opposite && lastHbond.opposite)
            {
                // starting with C, ending with C
                int cycles = riEnd - riStart;
                return "C" + "NcC".Repeat(cycles);
            }
            else if (firstHbond.opposite && !lastHbond.opposite)
            {
                // starting with C, ending with N
                int cycles = riEnd - riStart - 1;
                return "C" + "NcC".Repeat(cycles) + "N";
            }
            else if (!firstHbond.opposite && lastHbond.opposite)
            {
                // starting with N, ending with C
                int cycles = riEnd - riStart + 1;
                return "NcC".Repeat(cycles);
            }
            else
            {
                // starting with N, ending with N
                int cycles = riEnd - riStart;
                return "NcC".Repeat(cycles) + "N";
            }
        }
    }
}


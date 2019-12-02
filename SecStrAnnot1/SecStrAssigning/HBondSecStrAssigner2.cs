using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;
using protein.SecStrAssigning.Helpers;
using Cif.Tables;

namespace protein.SecStrAssigning
{
    public class HBondSecStrAssigner2 : ISecStrAssigner
    {
        private const int NONEXISTING_RESIDUE_INDEX = -1;

        private Model model;
        private Residue[] residues;
        private int nResidues;
        private HashSet<(int donor, int acceptor)> hBonds;
        private Dictionary<int,List<int>> acceptorsOf;

        public HBondSecStrAssigner2(Protein protein, double hBondEnergyCutoff){
            // this.DetectSheets=true;
            // this.DetectHelices=true;
            var hydrogenAdder = new HydrogenAdders.DsspLikeAmideHydrogenAdder();
            protein = hydrogenAdder.AddHydrogens(protein); // removing residues without C-alpha is probably not needed here
            this.model = protein.Model;
            this.residues = protein.GetResidues().ToArray();
            var hBondFinder = new BoxingHBondFinder(this.residues, hBondEnergyCutoff);
            this.nResidues = model.Residues.Count;
            this.hBonds = new HashSet<(int donor, int acceptor)>();
            this.acceptorsOf = new Dictionary<int,List<int>>();
            for (int don = 0; don < nResidues; don++){
                var acceptors = hBondFinder.FindHAcceptors(don);
                acceptorsOf[don] = acceptors;
                foreach (int acc in acceptors){
                    hBonds.Add((don, acc));
                }
            }
        }

        public String GetDescription(){
            return "hydrogen-bond-based method (similar to DSSP)";
        }

        public SecStrAssignment GetSecStrAssignment()
        {
            BuildBetaGraph();
            throw new NotImplementedException();
        }

        private void BuildBetaGraph(){
            var microstrand_firstZ = new List<int>();
            var microstrand_lastZ = new List<int>();
            var microstrand_fragment = new List<int>();
            var microstrand_microladder = new List<int>();
            var microstrand_type = new List<MicrostrandType>();
            var microstrand_macrostrand = new List<int>();
            int microstrand_count = 0;

            var microladder_microstrand0 = new List<int>();
            var microladder_microstrand1 = new List<int>();
            var microladder_type0 = new List<MicrostrandType>();
            var microladder_macroladder = new List<int>();
            var microladder_count = 0;

            var macrostrand_firstZ = new List<int>();
            var macrostrand_lastZ = new List<int>();
            var macrostrand_fragment = new List<int>();
            var macrostrand_microstrands = new List<List<int>>();
            var macrostrand_macroladders = new List<List<int>>();
            var macrostrand_sheet = new List<int>();
            var macrostrand_count = 0;
            
            var macroladder_microladders = new List<List<int>>();
            var macroladder_macrostrand0 = new List<int>();
            var macroladder_macrostrand1 = new List<int>();
            var macroladder_type = new List<MacroladderType>();
            var macroladder_count = 0;

            var sheet_macrostrands = new List<List<int>>();
            var sheet_macroladders = new List<List<int>>();
            int sheet_count = 0;

            List<Ladder> ladders = FindLadders();
            Ladder[] parallelLadders = ladders.Where(l => l.type0 == MicrostrandType.REGULAR_PARALEL).ToArray();
            Ladder[] antiparallelLadders = ladders.Where(l => l.type0 == MicrostrandType.REGULAR_ANTIPARALEL).ToArray();
            List<Ladder> parallelLaddersAndBulges = IncludeParallelBulges(parallelLadders);
            List<Ladder> antiparallelLaddersAndBulges = IncludeAntiparallelBulges(antiparallelLadders);
            Ladder[] microladders = parallelLaddersAndBulges.Concat(antiparallelLaddersAndBulges).ToArray();
            
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug("Microladders:");
                foreach (Ladder ladder in microladders){
                    Lib.WriteLineDebug(RepresentLadder(ladder));
                }
            }
            
            var microstrands = new List<(int startZeta, int endZeta, int ladderIndex, int side)>(capacity: 2 * (parallelLaddersAndBulges.Count + antiparallelLaddersAndBulges.Count));
            int nParallel = parallelLaddersAndBulges.Count;
            int nAntiparallel = antiparallelLaddersAndBulges.Count;
            int nMicroladders = nParallel + nAntiparallel;
            // Adding microstrands of parallel ladders and bulges
            for (int i = 0; i < nParallel; i++){
                Ladder ladder = microladders[i];
                microstrands.Add((ladder.firstHbond.Zeta0, ladder.lastHbond.Zeta0, i, 0));
                microstrands.Add((ladder.firstHbond.Zeta1, ladder.lastHbond.Zeta1, i, 1));
            }
            // Adding microstrands of antiparallel ladders and bulges
            for (int i = nParallel; i < nParallel+nAntiparallel; i++){
                Ladder ladder = microladders[i];
                microstrands.Add((ladder.firstHbond.Zeta0, ladder.lastHbond.Zeta0, i, 0));
                microstrands.Add((ladder.lastHbond.Zeta1, ladder.firstHbond.Zeta1, i, 1));
            }
            int nMicrostrands = microstrands.Count;
            var sortedMicrostrandIndices = Enumerable.Range(0, nMicrostrands).OrderBy(i => microstrands[i]).ToArray();
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug("Microstrands:");
                for (int i = 0; i < nMicrostrands; i++){
                    var ms = microstrands[sortedMicrostrandIndices[i]];
                    Lib.WriteLineDebug(ms.ToString());
                }
            }
            // Building macrostrands
            var macrostrands = new List<(int startZeta, int endZeta, List<int> microstrands)>(capacity: nMicrostrands);
            for (int i = 0; i < nMicrostrands; i++){
                bool starting = (i == 0);
                int microstrandIndex = sortedMicrostrandIndices[i];
                var microstrand = microstrands[microstrandIndex];
                var lastMacrostrand = starting ? (0, 0, null) : macrostrands[macrostrands.Count-1];
                if (!starting && DoOverlap(lastMacrostrand.startZeta, lastMacrostrand.endZeta, microstrand.startZeta, microstrand.endZeta)){
                    Lib.WriteLineDebug($"Joining macro {lastMacrostrand} and micro {microstrand}");
                    lastMacrostrand.microstrands.Add(microstrandIndex);
                    if (microstrand.endZeta > lastMacrostrand.endZeta){
                        macrostrands[macrostrands.Count-1] = (lastMacrostrand.startZeta, microstrand.endZeta, lastMacrostrand.microstrands);
                    }
                } else {
                    macrostrands.Add((microstrand.startZeta, microstrand.endZeta, new List<int>{microstrandIndex}));
                }
            }
            int nMacrostrands = macrostrands.Count;
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug("Macrostrands:");
                foreach (var ms in macrostrands){
                    int startResIdx = ZToResidueIndex(ms.startZeta);
                    int endResIdx = ZToResidueIndex(ms.endZeta);
                    string chainId = model.Chains.Id[model.Residues.ChainIndex[startResIdx]];
                    int startSeqId = model.Residues.SeqNumber[startResIdx];
                    int endSeqId = model.Residues.SeqNumber[endResIdx];
                    Lib.WriteLineDebug($"{chainId} {startSeqId}-{endSeqId} {ms.ToString()} " + Lib.EnumerateWithCommas(ms.microstrands));
                }
            }
            int[] microstrand2macrostrand = new int[nMicrostrands];
            for (int i = 0; i < nMacrostrands; i++){
                foreach (int microstrandIndex in macrostrands[i].microstrands){
                    microstrand2macrostrand[microstrandIndex] = i;
                }
            }
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug(Lib.EnumerateWithCommas(microstrand2macrostrand));
            }
            var macroladders = new List<(int macrostrand0, int macrostrand1, MacroladderType type)>();
            for (int i = 0; i < nMicroladders; i++){
                bool parallel = i < nParallel;
                MacroladderType type = parallel ? MacroladderType.PARALLEL : MacroladderType.ANTIPARALLEL;
                int mis0 = 2*i;
                int mis1 = 2*i + 1;
                int mas0 = microstrand2macrostrand[mis0];
                int mas1 = microstrand2macrostrand[mis1];
                macroladders.Add((mas0, mas1, type));
            }
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug("Macroladders:");
                Lib.WriteLineDebug(Lib.EnumerateWithSeparators(macroladders, "\n"));
            }
            var edges = macroladders.Select(mal => (mal.macrostrand0, mal.macrostrand1)).ToList();
            var sheets = ConnectedComponents(edges);
            if (Lib.DoWriteDebug){
                Lib.WriteLineDebug("Sheets:");
                foreach (var sheet in sheets){
                    Lib.WriteLineDebug(Lib.EnumerateWithCommas(sheet));
                }
                Lib.WriteLineDebug("");
            }
            // TODO continue here
            // TODO test on 2axt,EA
            
            


        }

        private List<Ladder> FindLadders(){
            List<Ladder> ladders = new List<Ladder>();
            for (int don = 0; don < nResidues; don++){
                foreach (int acc in acceptorsOf[don]){
                    if (acc <= don){
                        continue;
                        // this avoids: 1) recognizing H-bond i->i-2 as start of a ladder, 2) recognizing each ladder twice
                    }
                    // antiparallel
                    Hbond current = new Hbond(don, acc, false);
                    Hbond previous = PreviousAntiparallel(current);
                    Hbond preprevious = PreviousAntiparallel(previous);
                    Hbond next = NextAntiparallel(current);
                    bool startsMotifA = Exists(next) && !Exists(previous);
                    bool startsMotifB = Exists(previous) && !Exists(preprevious);
                    if (startsMotifA) {
                        Ladder newLadder = ExtendAntiparallelLadder(current, next);
                        ladders.Add(newLadder);
                    } else if (startsMotifB) {
                        Ladder newLadder = ExtendAntiparallelLadder(previous, current);
                        ladders.Add(newLadder);
                    }
                    // parallel
                    previous = PreviousParallel(current);
                    preprevious = PreviousParallel(previous);
                    next = NextParallel(current);
                    bool startsMotifC = Exists(next) && !Exists(previous);
                    bool startsMotifD = Exists(previous) && !Exists(preprevious);
                    if (startsMotifC) {
                        Ladder newLadder = ExtendParallelLadder(current, next);
                        ladders.Add(newLadder);
                    } else if (startsMotifD) {
                        Ladder newLadder = ExtendParallelLadder(previous, current);
                        ladders.Add(newLadder);
                    }
                }
            }
            return ladders;
        }
        
        private string RepresentLadder(Ladder ladder){
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

        private List<Ladder> IncludeParallelBulges(Ladder[] parallelLadders){  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int n = parallelLadders.Length;
            List<Ladder> result = new List<Ladder>(capacity: n);
            if (n == 0) {
                return result;
            }
            result.Add(parallelLadders[0]);
            for (int i = 1; i < n; i++){
                Ladder? potentialBulge = TryMakeParallelBulge(parallelLadders[i-1], parallelLadders[i]);
                if (potentialBulge != null){
                    result.Add(potentialBulge.Value);
                }
                result.Add(parallelLadders[i]);
            }
            return result;
        }
        private List<Ladder> IncludeAntiparallelBulges(Ladder[] antiparallelLadders){  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int n = antiparallelLadders.Length;
            List<Ladder> result = new List<Ladder>(capacity: n);
            if (n == 0) {
                return result;
            }
            result.Add(antiparallelLadders[0]);
            for (int i = 1; i < n; i++){
                Ladder? potentialBulge = TryMakeAntiparallelBulge(antiparallelLadders[i-1], antiparallelLadders[i]);
                if (potentialBulge != null){
                    result.Add(potentialBulge.Value);
                    Lib.WriteLineDebug("Found bulge");
                }
                result.Add(antiparallelLadders[i]);
            }
            return result;
        }

        private Ladder? TryMakeParallelBulge(Ladder firstLadder, Ladder secondLadder) {  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int zetaShift0 = secondLadder.firstHbond.Zeta0 - firstLadder.lastHbond.Zeta0;
            int zetaShift1 = secondLadder.firstHbond.Zeta1 - firstLadder.lastHbond.Zeta1;
            var fragmentIndex = model.Residues.FragmentIndex;
            bool sameFragments = fragmentIndex[firstLadder.firstHbond.r0] == fragmentIndex[secondLadder.firstHbond.r0] && fragmentIndex[firstLadder.firstHbond.r1] == fragmentIndex[secondLadder.firstHbond.r1];
            if (!sameFragments || zetaShift1 < 0){
                return null;
            }
            Ladder? result = null;
            if (zetaShift0 < zetaShift1 && zetaShift0 <= 4 && zetaShift1 <= 13){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_SHORT_SIDE);
            } else if (zetaShift0 > zetaShift1 && zetaShift0 <= 13 && zetaShift1 <= 4){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_LONG_SIDE);
            } else if (zetaShift0 == zetaShift1 && zetaShift0 <= 4){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_PARALLEL_EQUAL_SIDE);
            }
            // TODO distinguish bulge type
            return result;
        }
        private Ladder? TryMakeAntiparallelBulge(Ladder firstLadder, Ladder secondLadder) {  // Ladders must have strand0 in lower residue index and must sorted wrt strand0 residue index!
            int zetaShift0 = secondLadder.firstHbond.Zeta0 - firstLadder.lastHbond.Zeta0;
            int zetaShift1 = firstLadder.lastHbond.Zeta1 - secondLadder.firstHbond.Zeta1;
            var fragmentIndex = model.Residues.FragmentIndex;
            bool sameFragments = fragmentIndex[firstLadder.firstHbond.r0] == fragmentIndex[secondLadder.firstHbond.r0] && fragmentIndex[firstLadder.firstHbond.r1] == fragmentIndex[secondLadder.firstHbond.r1];
            if (!sameFragments || zetaShift1 < 0){
                return null;
            }
            Ladder? result = null;
            if (zetaShift0 < zetaShift1 && zetaShift0 <= 4 && zetaShift1 <= 13){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_SHORT_SIDE);
            } else if (zetaShift0 > zetaShift1 && zetaShift0 <= 13 && zetaShift1 <= 4){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_LONG_SIDE);
            } else if (zetaShift0 == zetaShift1 && zetaShift0 <= 4){
                result = new Ladder(firstLadder.lastHbond, secondLadder.firstHbond, MicrostrandType.BULGE_UNSPECIFIED_ANTIPARALLEL_EQUAL_SIDE);
            }
            // TODO distinguish bulge type
            if (result != null){
                Lib.WriteLineDebug($"{zetaShift0} {zetaShift1}");
            }
            return result;
        }
        
        private List<List<int>> ConnectedComponents(List<(int, int)> edges){
            var components = new List<List<int>>();
            var neighbourLists = new Dictionary<int, List<int>>(capacity: 2*edges.Count);
            foreach ((int u, int v) in edges){
                if (!neighbourLists.ContainsKey(u)) neighbourLists[u] = new List<int>();
                if (!neighbourLists.ContainsKey(v)) neighbourLists[v] = new List<int>();
                neighbourLists[u].Add(v);
                neighbourLists[v].Add(u);
            }
            var remainingVertices = new HashSet<int>(neighbourLists.Keys);
            while (remainingVertices.Count > 0){
                int seed = remainingVertices.First();
                var outList = new List<int>();
                FindComponentDFS(seed, neighbourLists, remainingVertices, outList);
                outList.Sort();
                components.Add(outList);
            }
            components = components.OrderBy(comp => comp[0]).ToList();
            return components;
        }

        private void FindComponentDFS(int seed, Dictionary<int, List<int>> neighbourLists, HashSet<int> remainingVertices, List<int> outList){
            remainingVertices.Remove(seed);
            outList.Add(seed);
            foreach (int neighbour in neighbourLists[seed]){
                if (remainingVertices.Contains(neighbour)){
                    FindComponentDFS(neighbour, neighbourLists, remainingVertices, outList);
                }
            }
        }
        
        private const int DONOR_Z_CHANGE = -1;
        private const int ACCEPTOR_Z_CHANGE = 1;

        // private int ZDonor(int donor) => 3*donor + DONOR_Z_CHANGE;
        // private int ZAcceptor(int acceptor) => 3*acceptor + ACCEPTOR_Z_CHANGE;

        // private bool ZIsDonor(int z) => z % 3 != ACCEPTOR_Z_CHANGE;
        // private bool ZIsAcceptor(int z) => z % 3 == ACCEPTOR_Z_CHANGE;
        private int ZToResidueIndex(int z) => (z-DONOR_Z_CHANGE) / 3;

        private int Move(int residueIndex, int move){
            int candidateIndex = residueIndex + move;
            if (candidateIndex >= 0 && candidateIndex < nResidues && model.Residues.FragmentIndex[residueIndex] == model.Residues.FragmentIndex[candidateIndex]){
                return candidateIndex;
            } else {
                return NONEXISTING_RESIDUE_INDEX;
            }
        }
        private struct Hbond {
            public int r0;  // residue index
            public int r1;  // residue index
            public bool opposite;  // if r1 is donor
            public Hbond(int r0, int r1, bool opposite){
                this.r0 = r0;
                this.r1 = r1;
                this.opposite = opposite;
            }
            public (int donor, int acceptor) Tuple => opposite ? (r1, r0) : (r0, r1);
            public int Zeta0 => 3*r0 + (opposite ? ACCEPTOR_Z_CHANGE : DONOR_Z_CHANGE);
            public int Zeta1 => 3*r1 + (opposite ? DONOR_Z_CHANGE : ACCEPTOR_Z_CHANGE);
        }
        private struct Ladder {
            public Hbond firstHbond;
            public Hbond lastHbond;
            public MicrostrandType type0;
            public Ladder(Hbond firstHbond, Hbond lastHbond, MicrostrandType type0){
                this.firstHbond = firstHbond;
                this.lastHbond = lastHbond;
                this.type0 = type0;
            }
        }
        private Hbond NextAntiparallel(Hbond hbond){
            if (hbond.opposite){
                return new Hbond(Move(hbond.r0, +2), Move(hbond.r1, -2), false);
            } else {
                return new Hbond(hbond.r0, hbond.r1, true);
            }
        }
        private Hbond PreviousAntiparallel(Hbond hbond){
            if (hbond.opposite){
                return new Hbond(hbond.r0, hbond.r1, false);
            } else {
                return new Hbond(Move(hbond.r0, -2), Move(hbond.r1, +2), true);
            }
        }
        private Hbond NextParallel(Hbond hbond){
            if (hbond.opposite){
                return new Hbond(Move(hbond.r0, +2), hbond.r1, false);
            } else {
                return new Hbond(hbond.r0, Move(hbond.r1, +2), true);
            }
        }
        private Hbond PreviousParallel(Hbond hbond){
            if (hbond.opposite){
                return new Hbond(hbond.r0, Move(hbond.r1, -2), false);
            } else {
                return new Hbond(Move(hbond.r0, -2), hbond.r1, true);
            }
        }
        private bool Exists(Hbond hbond){
            return this.hBonds.Contains(hbond.Tuple);
        }
        private Ladder ExtendAntiparallelLadder(Hbond first, Hbond last){
            while (true){
                Hbond next = NextAntiparallel(last);
                if (Exists(next) && next.r0 < next.r1){ // The second condition prevents beta-hairpin from extending over the end
                    last = next;   
                } else {
                    break;
                }
            }
            return new Ladder(first, last, MicrostrandType.REGULAR_ANTIPARALEL);
        }
        private Ladder ExtendParallelLadder(Hbond first, Hbond last){
            while (true){
                Hbond next = NextParallel(last);
                if (Exists(next))
                    last = next;
                else
                    break;
            }
            return new Ladder(first, last, MicrostrandType.REGULAR_PARALEL);
        }
        
        private bool DoOverlap(int startZeta1, int endZeta1, int startZeta2, int endZeta2){
            return startZeta2 <= endZeta1 - HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING 
                && startZeta1 <= endZeta2 - HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING;
        }
    }
}


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using SecStrAnnotator2.Utils;
using Cif.Components;
using protein.Libraries;
using protein.Sses;
using protein.HydrogenAdding;
using protein.SecStrAssigning.Helpers;

namespace protein.SecStrAssigning
{
    public class HBondSecStrAssigner : ISecStrAssigner
    {
        public bool DetectSheets { get; set; }
        public bool DetectHelices { get; set; }

        private Residue[] residues;
        List<int>[] ignoreHBondsTo_Antiparallel;
        List<int>[] ignoreHBondsFrom_Antiparallel;
        List<int>[] ignoreHBondsTo_Parallel;
        List<int>[] ignoreHBondsFrom_Parallel;

        private IHydrogenAdder hydrogenAdder;
        private IHBondFinder hBondFinder;

        public HBondSecStrAssigner(Protein protein, double hBondEnergyCutoff)
        {
            this.DetectSheets = true;
            this.DetectHelices = true;
            this.hydrogenAdder = new DsspLikeAmideHydrogenAdder();
            this.residues = hydrogenAdder.AddHydrogens(protein).GetResidues().ToArray(); // removing residues without C-alpha is probably not needed here
            // Cif.Libraries.Lib.WriteLineDebug($"HBondSecStrAssigner(): {this.residues.Length}");
            ignoreHBondsTo_Antiparallel = new List<int>[residues.Length];
            ignoreHBondsFrom_Antiparallel = new List<int>[residues.Length];
            ignoreHBondsTo_Parallel = new List<int>[residues.Length];
            ignoreHBondsFrom_Parallel = new List<int>[residues.Length];
            this.hBondFinder = new SparseBoxingHBondFinder(this.residues, hBondEnergyCutoff);
        }


        private String Ladder2String(BetaLadder ladder)
        {
            char fd = ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1 ? 'v' : '^';
            char ld = ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1 ? 'v' : '^';
            return String.Format("[{0} {1}-{2} : {3} {4}-{5} {6}{7}]", residues[ladder.Start0].ChainId, residues[ladder.Start0].SeqNumber, residues[ladder.End0].SeqNumber, residues[ladder.Start1].ChainId, residues[ladder.Start1].SeqNumber, residues[ladder.End1].SeqNumber, fd, ld);
        }

        private int Residue2After(int i)
        {
            return ResidueXAfter(i, 2);
        }

        private int Residue2Before(int i)
        {
            return ResidueXBefore(i, 2);
        }

        private int ResidueXAfter(int resIndex, int x)
        {
            if (x < 0)
                return ResidueXBefore(resIndex, -x);
            if (resIndex < 0)
                return -1; //does not exist
            int candidateIndex = resIndex + x;
            if (candidateIndex < residues.Length
                    && residues[candidateIndex].ChainId == residues[resIndex].ChainId
                    && residues[candidateIndex].SeqNumber == residues[resIndex].SeqNumber + x)
            {
                return candidateIndex;
            }
            else
            {
                return -1; //does not exist
            }
            // for (int i = x; i >= 0; i--) {
            //     if (resIndex + i < residues.Length 
            //             && residues[resIndex + i].ChainId == residues[resIndex].ChainId 
            //             && residues[resIndex + i].SeqNumber == residues[resIndex].SeqNumber + x){
            //         return resIndex + i;
            //     }
            // }
            // return -1; //does not exist
        }

        private int ResidueXBefore(int resIndex, int x)
        {
            if (x < 0)
                return ResidueXAfter(resIndex, -x);
            if (resIndex < 0)
                return -1; //does not exist
            int candidateIndex = resIndex - x;
            if (candidateIndex >= 0
                    && residues[candidateIndex].ChainId == residues[resIndex].ChainId
                    && residues[candidateIndex].SeqNumber == residues[resIndex].SeqNumber - x)
            {
                return candidateIndex;
            }
            else
            {
                return -1; //does not exist
            }
            // for (int i = x; i >= 0; i--) {
            //     if (resIndex - i >= 0 
            //             && residues[resIndex - i].ChainId == residues[resIndex].ChainId 
            //             && residues[resIndex - i].SeqNumber == residues[resIndex].SeqNumber - x){
            //         return resIndex - i;
            //     }
            // }
            // return -1; //does not exist
        }

        private bool CanExtendLadder(BetaLadder ladder)
        {
            if (ladder.Type == BetaLadder.LadderType.Antiparallel)
            {
                // Antiparallel
                if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
                {
                    int res2afterEnd0 = ResidueXAfter(ladder.End0, 2);
                    return hBondFinder.IsHBond(ladder.Start1, ladder.End0) && ladder.Start1 > res2afterEnd0 && res2afterEnd0 > 0;
                    // return hBondFinder.IsHBond(ladder.Start1, ladder.End0) && (ladder.Start1 > ResidueXAfter(ladder.End0, 2));
                }
                else
                {
                    int i = Residue2After(ladder.End0);
                    int j = Residue2Before(ladder.Start1);
                    return i >= 0 && j >= 0 && hBondFinder.IsHBond(i, j);
                }
            }
            else
            {
                // Parallel
                if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
                {
                    int j = Residue2After(ladder.End1);
                    return j >= 0 && hBondFinder.IsHBond(j, ladder.End0);
                }
                else
                {
                    int i = Residue2After(ladder.End0);
                    return i >= 0 && hBondFinder.IsHBond(i, ladder.End1);
                }
            }
        }

        private void AddSafely<T>(ref List<T> list, T newElement)
        {
            if (list == null)
                list = new List<T>();
            list.Add(newElement);
        }

        private void IgnoreRedundantMotifs(BetaLadder ladder)
        {
            if (ladder.Type == BetaLadder.LadderType.Antiparallel)
            {
                if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    AddSafely(ref ignoreHBondsTo_Antiparallel[ladder.Start0], ladder.End1);
                for (int i = 2; i <= ladder.End0 - ladder.Start0; i += 2)
                    AddSafely(ref ignoreHBondsTo_Antiparallel[ladder.Start0 + i], ladder.End1 - i);
                for (int i = 0; i < ladder.End0 - ladder.Start0; i += 2)
                    AddSafely(ref ignoreHBondsFrom_Antiparallel[ladder.Start0 + i], ladder.End1 - i);
                if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
                    AddSafely(ref ignoreHBondsFrom_Antiparallel[ladder.End0], ladder.Start1);
            }
            else
            { //Parallel
                if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                {
                    for (int i = 0; i <= ladder.End0 - ladder.Start0; i += 2)
                    {
                        AddSafely(ref ignoreHBondsTo_Parallel[ladder.Start0 + i], ladder.Start1 + i);
                    }
                    for (int i = 0; i < ladder.End0 - ladder.Start0; i += 2)
                    {
                        AddSafely(ref ignoreHBondsFrom_Parallel[ladder.Start0 + i], ladder.Start1 + 2 + i);
                    }
                    if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
                        AddSafely(ref ignoreHBondsFrom_Parallel[ladder.End0], ladder.End1);
                }
                else
                {
                    for (int i = 2; i <= ladder.End0 - ladder.Start0; i += 2)
                    {
                        AddSafely(ref ignoreHBondsTo_Parallel[ladder.Start0 + i], ladder.Start1 - 2 + i);
                    }
                    for (int i = 0; i < ladder.End0 - ladder.Start0; i += 2)
                    {
                        AddSafely(ref ignoreHBondsFrom_Parallel[ladder.Start0 + i], ladder.Start1 + i);
                    }
                    if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
                        AddSafely(ref ignoreHBondsFrom_Parallel[ladder.End0], ladder.End1);
                }
            }
        }

        private int CheckLadderAndCountHBonds(BetaLadder ladder)
        {
            if (ladder.Type == BetaLadder.LadderType.Antiparallel)
            {
                if (ladder.End0 - ladder.Start0 != ladder.End1 - ladder.Start1)
                    throw new SecStrAssignmentException("Antiparallel beta-ladder: Strands have different lengths.");
                if ((ladder.End0 - ladder.Start0) % 2 != 0)
                    throw new SecStrAssignmentException("Antiparallel beta-ladder: Strands have even lengths.");
                int numHBonds = (ladder.End0 - ladder.Start0 + 2)
                                - (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0 ? 1 : 0)
                                - (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1 ? 1 : 0);
                if (numHBonds < 2)
                    throw new SecStrAssignmentException("Antiparallel beta-ladder: Less than 2 hydrogen bonds.");
                return numHBonds;
            }
            else
            { //Parallel
                if (ladder.End0 - ladder.Start0 - ladder.End1 + ladder.Start1
                    != (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0 ? 2 : 0) - (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0 ? 2 : 0))
                    //Lib.WriteLineDebug("parallel[{0}-{1}|{2}-{3}|{4},{5}]",residues[ladder.Start0],residues[ladder.End0],residues[ladder.Start1],residues[ladder.End1],ladder.FirstHBondDirection,ladder.LastHBondDirection);
                    throw new SecStrAssignmentException("Parallel beta-ladder: Strands have incompatible lengths.");
                int numHBonds = 1 + (ladder.End0 - ladder.Start0 + ladder.End1 - ladder.Start1) / 2;
                if (numHBonds < 2)
                    throw new SecStrAssignmentException("Parallel beta-ladder: Less than 2 hydrogen bonds.");
                return numHBonds;
            }
        }

        private bool DoOverlap(BetaStrandInSheet a, BetaStrandInSheet b)
        {
            if (a.SSE.ChainID == b.SSE.ChainID)
            {
                foreach (BetaLadder la in a.UpLadders.Union(a.DownLadders))
                    foreach (BetaLadder lb in b.UpLadders.Union(b.DownLadders))
                        if (la.ZEnd0 - lb.ZStart0 >= HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING && lb.ZEnd0 - la.ZStart0 >= HBondSSAConstants.MIN_Z_OVERLAP_FOR_JOINING)
                            return true;
            }
            return false;
            //return a.SSE.ChainID==b.SSE.ChainID && a.SSE.End >= b.SSE.Start + MIN_OVERLAP_FOR_JOINING && b.SSE.End >= a.SSE.Start + MIN_OVERLAP_FOR_JOINING;
        }

        private bool AreShortSideOfBetaBulge(BetaStrandInSheet a, BetaStrandInSheet b)
        {
            if (a.SSE.ChainID != b.SSE.ChainID)
                return false;
            foreach (BetaLadder la in a.UpLadders.Union(a.DownLadders))
            {
                foreach (BetaLadder lb in b.UpLadders.Union(b.DownLadders))
                {
                    if (BuildBetaBulge(la, lb) != null)
                        return true;
                }
            }
            return false;
        }

        private bool AreLongSideOfBetaBulge(BetaStrandInSheet a, BetaStrandInSheet b)
        {
            if (a.SSE.ChainID != b.SSE.ChainID)
                return false;
            foreach (BetaLadder la in a.UpLadders.Union(a.DownLadders))
            {
                foreach (BetaLadder lb in b.UpLadders.Union(b.DownLadders))
                {
                    if (BuildBetaBulge(la.Inverted(), lb.Inverted()) != null)
                        return true;
                }
            }
            return false;
        }

        /** Tries to build a beta-bulge so than the la.Strand0 and lb.Strand0 form the SHORT side of the bulge. */
        private BetaBulge BuildBetaBulge(BetaLadder la, BetaLadder lb)
        {
            return BuildBetaBulgeWithFirstLadderFirst(la, lb) ?? BuildBetaBulgeWithFirstLadderFirst(lb, la);
        }
        /** Tries to build a beta-bulge so than the la.Strand0 and lb.Strand0 form the SHORT side of the bulge. */
        private BetaBulge BuildBetaBulgeWithFirstLadderFirst(BetaLadder la, BetaLadder lb)
        {
            if (residues[la.Start1].ChainId == residues[lb.Start1].ChainId)
            {
                if (la.Type == BetaLadder.LadderType.Antiparallel && lb.Type == BetaLadder.LadderType.Antiparallel)
                {
                    //Antiparallel types
                    if (lb.Start0 == la.End0 && la.Start1 == ResidueXAfter(lb.End1, 1))
                    {
                        // classic bulge a-b(type A12)
                        return new BetaBulge(BetaBulge.BulgeType.Classic, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                    if (lb.Start0 == ResidueXAfter(la.End0, 2) && la.Start1 == ResidueXAfter(lb.End1, 3)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // wide bulge a-b(type A34)
                        return new BetaBulge(BetaBulge.BulgeType.Wide, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                    if (lb.Start0 == ResidueXAfter(la.End0, 1) && la.Start1 == ResidueXAfter(lb.End1, 1)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // bulge-not-bulge a-b(type A22) found in 1gei
                        return new BetaBulge(BetaBulge.BulgeType.Antiparallel22, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                    if (HBondSSAConstants.ALLOW_BULGE_A33 && lb.Start0 == ResidueXAfter(la.End0, 2) && la.Start1 == ResidueXAfter(lb.End1, 2)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // bulge-not-bulge a-b(type A33)
                        return new BetaBulge(BetaBulge.BulgeType.Antiparallel33, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                    if (lb.Start0 == la.End0 && la.Start1 == ResidueXAfter(lb.End1, 4)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
                    {
                        // arch-bulge a-b(type A15) found in 1gjm
                        return new BetaBulge(BetaBulge.BulgeType.Antiparallel15, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                    if (lb.Start0 == ResidueXAfter(la.End0, 1) && la.Start1 == ResidueXAfter(lb.End1, 2)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // (type 23) found in 3dbg
                        return new BetaBulge(BetaBulge.BulgeType.Antiparallel23, la.End0, lb.Start0, lb.End1, la.Start1);
                    }
                }
                if (la.Type == BetaLadder.LadderType.Parallel && lb.Type == BetaLadder.LadderType.Parallel)
                {
                    //Parallel types
                    if (lb.Start0 == la.End0 && lb.Start1 == ResidueXAfter(la.End1, 3)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
                    {
                        // parallel bulge(type P14)
                        return new BetaBulge(BetaBulge.BulgeType.Parallel14, la.End0, lb.Start0, la.End1, lb.Start1);
                    }
                    if (lb.Start0 == ResidueXAfter(la.End0, 2) && lb.Start1 == ResidueXAfter(la.End1, 1)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // parallel bulge(type P32)
                        return new BetaBulge(BetaBulge.BulgeType.Parallel32, la.End0, lb.Start0, la.End1, lb.Start1);
                    }
                    if (lb.Start0 == la.End0 && lb.Start1 == ResidueXAfter(la.End1, 2)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // parallel bulge(type P13)
                        return new BetaBulge(BetaBulge.BulgeType.Parallel13, la.End0, lb.Start0, la.End1, lb.Start1);
                    }
                    if (lb.Start0 == ResidueXAfter(la.End0, 2) && lb.Start1 == ResidueXAfter(la.End1, 2)
                        && la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                    {
                        // parallel bulge(type P33)
                        return new BetaBulge(BetaBulge.BulgeType.Parallel33, la.End0, lb.Start0, la.End1, lb.Start1);
                    }
                }
            }
            return null;
        }

        private void InvertSheetLevels(BetaStrandInSheet strand)
        {
            Lib2.WriteLineDebug("Invert sheet: {0}.", strand.SheetId);
            strand.DFS(s =>
            {
                //s.Level = -s.Level;
                s.EvenUp = !s.EvenUp;
                var aux = s.DownNeighbours;
                s.DownNeighbours = s.UpNeighbours;
                s.UpNeighbours = aux;
                var aux2 = s.DownLadders;
                s.DownLadders = s.UpLadders;
                s.UpLadders = aux2;
            });
        }

        private void SetIdToSheet(BetaStrandInSheet strand, int newID)
        {
            Lib2.WriteLineDebug("Set ID: {0} to {1}.", strand.SheetId, newID);
            strand.DFS(s =>
            {
                s.SheetId = newID;
            });
            Lib2.WriteLineDebug("");
        }

        private void LinkStrands(BetaStrandInSheet lower, BetaStrandInSheet upper, BetaLadder ladder)
        {
            Lib2.WriteLineDebug("Linking strands {0}({1}) and {2}({3}) {4}.", lower.SSE.Label, lower.SheetId, upper.SSE.Label, upper.SheetId, Ladder2String(ladder));
            if (lower == upper)
                Lib2.WriteWarning("Trying to link a beta-strand to itself.");
            if (lower.SheetId != upper.SheetId)
                throw new ArgumentException();
            lower.UpNeighbours.Add(upper);
            lower.UpLadders.Add(ladder);
            upper.DownNeighbours.Add(lower);
            upper.DownLadders.Add(ladder.Inverted());
        }

        private void UnlinkStrands(BetaStrandInSheet lower, BetaStrandInSheet upper)
        {
            int i = lower.UpNeighbours.IndexOf(upper);
            int j = upper.DownNeighbours.IndexOf(lower);
            Lib2.WriteLineDebug("Unlinking: upper={0}, lower={1}", i, j);
            if (lower == upper)
                Lib2.WriteWarning("Trying to unlink a beta-strand from itself.");
            lower.UpNeighbours.RemoveAt(i);
            lower.UpLadders.RemoveAt(i);
            upper.DownNeighbours.RemoveAt(j);
            upper.DownLadders.RemoveAt(j);
        }

        private void AddVertexAndPossiblyMerge(List<BetaStrandInSheet> vertices, BetaStrandInSheet newVertex)
        {
            Lib2.WriteLineDebug("Adding vertex {0}.", newVertex);
            vertices.Add(newVertex);
            int iNew = vertices.Count - 1;
            for (int i = vertices.Count - 2; i >= 0; i--)
            {
                if (DoOverlap(vertices[i], vertices[iNew]) || AreShortSideOfBetaBulge(vertices[i], vertices[iNew]) || AreLongSideOfBetaBulge(vertices[i], vertices[iNew]))
                {
                    bool merged = TryMergeSheetsByStrands(vertices[i], vertices[iNew]);
                    if (merged)
                    {
                        vertices.RemoveAt(iNew);
                        iNew = i;
                    }
                }
            }
        }

        private bool TryMergeSheetsByStrands(BetaStrandInSheet extendedStrand, BetaStrandInSheet extendingStrand)
        {
            Lib2.WriteLineDebug("Merging strands {0}({1}) and {2}({3}).", extendedStrand.SSE.Label, extendedStrand.SheetId, extendingStrand.SSE.Label, extendingStrand.SheetId);
            if (extendedStrand == extendingStrand)
                throw new ArgumentException();
            string chain = extendedStrand.SSE.ChainID;
            int start = Math.Min(extendedStrand.SSE.Start, extendingStrand.SSE.Start);
            int end = Math.Max(extendedStrand.SSE.End, extendingStrand.SSE.End);

            if (extendedStrand.SheetId == extendingStrand.SheetId)
            {
                // in the same sheet
                if (extendedStrand.EvenUp != extendingStrand.EvenUp)
                {
                    Lib2.WriteLineDebug("Important! Cycle with inconsistent direction would arise in beta-strand graph(around chain {0} {1}-{2})! Skipping merging strands.", chain, start, end);
                    //Lib.WriteWarning("Important! Cycle with inconsistent direction would arise in beta-strand graph(around chain {0} {1}-{2})! Skipping merging strands.", chain, start, end);
                    // return false;
                }
                else
                {
                    Lib2.WriteLineDebug("Cycle detected in beta-strand graph(around chain {0} {1}-{2}).", chain, start, end);
                }
                /*if (extendedStrand.Level != extendingStrand.Level) {
                    Lib.WriteWarning("Cycle with inconsistent level(might be a barrel) detected in beta-strand graph(around chain {0} {1}-{2})! Skipping merging strands.", extendedStrand.SSE.ChainID, Math.Min(extendedStrand.SSE.Start, extendingStrand.SSE.Start), Math.Max(extendedStrand.SSE.End, extendingStrand.SSE.End));
                }*/
            }
            else
            {
                // in different sheets
                if (extendedStrand.EvenUp != extendingStrand.EvenUp)
                    InvertSheetLevels(extendingStrand);
                SetIdToSheet(extendingStrand, extendedStrand.SheetId);
            }

            extendedStrand.SSE.Start = start;
            extendedStrand.SSE.End = end;
            BetaStrandInSheet[] oldUpNeighbours = extendingStrand.UpNeighbours.ToArray();
            BetaLadder[] oldUpLadders = extendingStrand.UpLadders.ToArray();
            BetaStrandInSheet[] oldDownNeighbours = extendingStrand.DownNeighbours.ToArray();
            BetaLadder[] oldDownLadders = extendingStrand.DownLadders.ToArray();
            Lib2.WriteLineDebug("Extended: {0}", extendedStrand.SSE);
            Lib2.WriteLineDebug("Extending: {0}", extendingStrand.SSE);
            for (int i = 0; i < oldUpNeighbours.Length; i++)
            {
                UnlinkStrands(extendingStrand, oldUpNeighbours[i]);
                LinkStrands(extendedStrand, oldUpNeighbours[i], oldUpLadders[i]);
            }
            for (int i = 0; i < oldDownNeighbours.Length; i++)
            {
                UnlinkStrands(oldDownNeighbours[i], extendingStrand);
                LinkStrands(oldDownNeighbours[i], extendedStrand, oldDownLadders[i].Inverted());
            }
            return true;
        }

        private Sse GetStrand0(BetaLadder ladder)
        {
            string chainId = residues[ladder.Start0].ChainId;
            int start = residues[ladder.Start0].SeqNumber;
            int end = residues[ladder.End0].SeqNumber;
            if (HBondSSAConstants.STRANDS_BY_ALPHA)
            {
                if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
                    start++;
                if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
                    end--;
            }
            return new Sse(null, chainId, start, end, LadderSSEType(ladder), null);
        }

        private Sse GetStrand1(BetaLadder ladder)
        {
            string chainId = residues[ladder.Start1].ChainId;
            int start = residues[ladder.Start1].SeqNumber;
            int end = residues[ladder.End1].SeqNumber;
            if (HBondSSAConstants.STRANDS_BY_ALPHA)
            {
                if (ladder.Type == BetaLadder.LadderType.Antiparallel)
                {
                    if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
                        end--;
                    if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
                        start++;
                }
                if (ladder.Type == BetaLadder.LadderType.Parallel)
                {
                    if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
                        start++;
                    if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
                        end--;
                }
            }
            return new Sse(null, chainId, start, end, LadderSSEType(ladder), null);
        }

        private Sse GetShortStrand(BetaBulge bulge)
        {
            if (HBondSSAConstants.BULGES_BY_ALPHA)
            {
                throw new NotImplementedException();
            }
            else
            {
                switch (bulge.Type)
                {
                    case BetaBulge.BulgeType.Classic:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL_CLASSIC_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Wide:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL_WIDE_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel22:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL22_DONOR_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel33:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL33_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel15:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL15_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel23:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_ANTIPARALLEL23_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Parallel14:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_PARALLEL14_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Parallel32:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_PARALLEL32_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Parallel13:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_PARALLEL13_SHORT_SIDE, null);
                    case BetaBulge.BulgeType.Parallel33:
                        return new Sse(null, residues[bulge.StartShort].ChainId, residues[bulge.StartShort].SeqNumber, residues[bulge.EndShort].SeqNumber, SseType.BULGE_PARALLEL33_SHORT_SIDE, null);
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        private Sse GetLongStrand(BetaBulge bulge)
        {
            if (HBondSSAConstants.BULGES_BY_ALPHA)
            {
                throw new NotImplementedException();
            }
            else
            {
                switch (bulge.Type)
                {
                    case BetaBulge.BulgeType.Classic:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALLEL_CLASSIC_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Wide:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALEL_WIDE_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel22:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALLEL22_ACCEPTOR_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel33:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALLEL33_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel15:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALLEL15_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Antiparallel23:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_ANTIPARALLEL23_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Parallel14:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_PARALLEL14_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Parallel32:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_PARALLEL32_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Parallel13:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_PARALLEL13_LONG_SIDE, null);
                    case BetaBulge.BulgeType.Parallel33:
                        return new Sse(null, residues[bulge.StartLong].ChainId, residues[bulge.StartLong].SeqNumber, residues[bulge.EndLong].SeqNumber, SseType.BULGE_PARALLEL33_LONG_SIDE, null);
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        private SseType LadderSSEType(BetaLadder ladder)
        {
            if (ladder.Start0 == ladder.Start1)
            {
                // C7 motif
                if (residues[ladder.End0].SeqNumber - residues[ladder.Start0].SeqNumber != 2)
                    throw new SecStrAssignmentException("C7 turn with more than 1 stabilizing H-bond.");
                else
                    return SseType.TURN_C7_TYPE;
            }
            else
            {
                // real beta ladder
                return (CheckLadderAndCountHBonds(ladder) > 2) ? SseType.SHEET_TYPE : SseType.ISOLATED_BETA_BRIDGE_TYPE;
            }
        }

        public SecStrAssignment GetSecStrAssignment()
        {
            if (DetectSheets)
            {
                if (DetectHelices)
                {
                    return SecStrAssignment.Combine(GetSheets(), GetHelices());
                }
                else
                {
                    return GetSheets();
                }
            }
            else
            {
                if (DetectHelices)
                {
                    return GetHelices();
                }
                else
                {
                    return new SecStrAssignment(new Sse[] { });
                }
            }
        }

        private SecStrAssignment GetHelices()
        {
            var gHelices = GetXHelices(3, SseType.HELIX_G_TYPE);
            var hHelices = GetXHelices(4, SseType.HELIX_H_TYPE);
            var iHelices = GetXHelices(5, SseType.HELIX_I_TYPE);
            return new SecStrAssignment(gHelices.Concat(hHelices).Concat(iHelices));
        }

        /** Returns the list of helices formed by(i+x -> i) H bonds. */
        private List<Sse> GetXHelices(int x, SseType assignedType)
        {
            List<Sse> helices = new List<Sse>();
            int currentStart = -1; // -1 = currently not making helix
            int currentEnd = -1;
            for (int i = 0; i < residues.Length; i++)
            {
                int j = ResidueXAfter(i, x);
                if (currentStart == -1 && hBondFinder.IsHBond(j, i))
                {
                    //start helix
                    currentStart = i;
                    currentEnd = j;
                }
                else if (currentStart >= 0 && hBondFinder.IsHBond(j, i))
                {
                    //continue helix
                    currentEnd = j;
                }
                else if (currentStart >= 0 && !hBondFinder.IsHBond(j, i))
                {
                    //finish or reject helix
                    if (i - currentStart >= HBondSSAConstants.MIN_HBONDS_PER_HELIX)
                    {
                        helices.Add(new Sse(null,
                            residues[currentStart].ChainId,
                            residues[currentStart].SeqNumber + (HBondSSAConstants.HELICES_BY_ALPHA ? 1 : 0),
                            residues[currentEnd].SeqNumber - (HBondSSAConstants.HELICES_BY_ALPHA ? 1 : 0),
                            assignedType, null));
                    }
                    currentStart = -1;
                }
            }
            return helices;
        }

        private List<BetaBulge> FindBulgesEachToEach(List<BetaLadder> ladders)
        {
            List<BetaBulge> bulges = new List<BetaBulge>();
            BetaBulge bulge;
            for (int i = 0; i < ladders.Count; i++)
            {
                for (int j = i + i; j < ladders.Count; j++)
                {
                    bulge = BuildBetaBulge(ladders[i], ladders[j]);
                    if (bulge != null) bulges.Add(bulge);
                    bulge = BuildBetaBulge(ladders[i], ladders[j].Inverted());
                    if (bulge != null) bulges.Add(bulge);
                    bulge = BuildBetaBulge(ladders[i].Inverted(), ladders[j]);
                    if (bulge != null) bulges.Add(bulge);
                    bulge = BuildBetaBulge(ladders[i].Inverted(), ladders[j].Inverted());
                    if (bulge != null) bulges.Add(bulge);
                }
            }
            return bulges;
        }
        private void BuildBetaGraph_New(List<BetaLadder> ladders)
        {
            List<BetaBulge> bulges = FindBulgesEachToEach(ladders);
            List<(int startZ, int endZ)> microstrands = new List<(int startZ, int endZ)>();
            foreach (BetaLadder ladder in ladders)
            {
                microstrands.Add((ladder.ZStart0, ladder.ZEnd0));
                microstrands.Add((ladder.ZStart1, ladder.ZEnd1));
            }
            foreach ((int i, int j) in microstrands)
            {

            }
        }

        private SecStrAssignment GetSheets()
        {
            // Find H-binding motifs
            List<BetaLadder> ladders = new List<BetaLadder>();
            List<(Residue, Residue)> hBonds = new List<(Residue, Residue)>();
            for (int i = 0; i < residues.Length; i++)
            {
                int nLaddersBefore = ladders.Count;
                List<int> hAcceptors = hBondFinder.FindHAcceptors(i).Where(j => j > i)/*.Where(j => ignoreHBondsTo[i] == null || !ignoreHBondsTo[i].Contains(j))*/.ToList();
                foreach (int j in hAcceptors)
                {
                    // try motif A
                    if ((ignoreHBondsTo_Antiparallel[i] == null || !ignoreHBondsTo_Antiparallel[i].Contains(j)) && hBondFinder.IsHBond(j, i))
                    {
                        BetaLadder ladder = new BetaLadder(BetaLadder.LadderType.Antiparallel, i, j, BetaLadder.HBondDirection.From0To1);
                        ladder.AddOneHBond();
                        while (CanExtendLadder(ladder))
                            ladder.AddOneHBond();
                        CheckLadderAndCountHBonds(ladder);
                        IgnoreRedundantMotifs(ladder);
                        ladders.Add(ladder);
                    }
                    // try motif C
                    if ((ignoreHBondsTo_Parallel[i] == null || !ignoreHBondsTo_Parallel[i].Contains(j)) && hBondFinder.IsHBond(Residue2After(j), i))
                    {
                        BetaLadder ladder = new BetaLadder(BetaLadder.LadderType.Parallel, i, j, BetaLadder.HBondDirection.From0To1);
                        ladder.AddOneHBond();
                        while (CanExtendLadder(ladder))
                            ladder.AddOneHBond();
                        CheckLadderAndCountHBonds(ladder);
                        IgnoreRedundantMotifs(ladder);
                        ladders.Add(ladder);
                    }
                }
                List<int> hDonors = hBondFinder.FindHDonors(i).Where(j => j > i)/*.Where(j => ignoreHBondsFrom[i] == null || !ignoreHBondsFrom[i].Contains(j))*/.ToList();
                foreach (int j in hDonors)
                {
                    // try motif B
                    if ((ignoreHBondsFrom_Antiparallel[i] == null || !ignoreHBondsFrom_Antiparallel[i].Contains(j)) && hBondFinder.IsHBond(Residue2After(i), Residue2Before(j)))
                    {
                        BetaLadder ladder = new BetaLadder(BetaLadder.LadderType.Antiparallel, i, j, BetaLadder.HBondDirection.From1To0);
                        ladder.AddOneHBond();
                        while (CanExtendLadder(ladder))
                            ladder.AddOneHBond();
                        CheckLadderAndCountHBonds(ladder);
                        IgnoreRedundantMotifs(ladder);
                        ladders.Add(ladder);
                    }
                    // try motif D
                    if ((ignoreHBondsFrom_Parallel[i] == null || !ignoreHBondsFrom_Parallel[i].Contains(j)) && hBondFinder.IsHBond(Residue2After(i), j))
                    {
                        BetaLadder ladder = new BetaLadder(BetaLadder.LadderType.Parallel, i, j, BetaLadder.HBondDirection.From1To0);
                        ladder.AddOneHBond();
                        while (CanExtendLadder(ladder))
                            ladder.AddOneHBond();
                        CheckLadderAndCountHBonds(ladder);
                        IgnoreRedundantMotifs(ladder);
                        ladders.Add(ladder);
                    }
                }
                if (ladders.Count - nLaddersBefore > 1)
                    Lib2.WriteLineDebug("More than one beta-ladder motif found from residue {0}.", residues[i].ToString(true));
                hBonds.AddRange(hAcceptors.Select(a => (residues[i], residues[a])));
                hBonds.AddRange(hDonors.Select(d => (residues[d], residues[i])));
            }

            foreach (BetaLadder l in ladders)
                if (l.Start0 == l.Start1 && (l.End0 - l.Start0) > 2)
                    Lib2.WriteWarning("Strange secondary structure in chain {0} {1}-{2}.", residues[l.Start0].ChainId, residues[l.Start0].SeqNumber, residues[l.End0].SeqNumber);

            List<Sse> c7Turns = ladders
                .Where(l => l.Start0 == l.Start1)
                .OrderBy(l => l.Start0)
                .Select(l => GetStrand0(l)/*new SSE("x", residues[l.Start0].ChainID, residues[l.Start0].ResSeq, residues[l.End0].ResSeq, SSE.TURN_C7_TYPE,null)*/)
                .ToList();

            // Build C7 turns and wiggles
            List<Sse> c7TurnsAndWiggles = new List<Sse>();
            string currentChain = null;
            int? firstResi = null;
            int? lastResi = null;
            int strandCounter = 0;
            foreach (Sse t in c7Turns)
            {
                if (t.ChainID == currentChain && t.End == lastResi + 1)
                {
                    // elongation of existing SSE
                    lastResi++;
                }
                else
                {
                    if (currentChain != null)
                    {
                        // finishing the old SSE
                        SseType type = (lastResi - firstResi > 2) ? SseType.WIGGLE_C7_TYPE : SseType.TURN_C7_TYPE;
                        c7TurnsAndWiggles.Add(new Sse($"{type.AsString()}{strandCounter++}", currentChain, (int)firstResi, (int)lastResi, type, null));
                    }
                    // starting a new SSE
                    currentChain = t.ChainID;
                    firstResi = t.Start;
                    lastResi = t.End;
                }
            }

            // DEBUG: Write some info about the ladders
            if (Lib2.DoWriteDebug)
            {
                foreach (BetaLadder ladder in ladders)
                {
                    if (ladder.Start0 != ladder.Start1)
                        Lib2.WriteInColor(ConsoleColor.Cyan, "*");
                    else
                        Lib2.WriteInColor(ConsoleColor.Red, "*");
                    Lib2.WriteLineDebug(Ladder2String(ladder));
                }
            }

            // Build beta-sheet graph - NEW VERSION
            List<BetaLadder> realLadders1 = ladders.Where(l => l.Start0 != l.Start1).Where(l => CheckLadderAndCountHBonds(l) >= HBondSSAConstants.MIN_HBONDS_PER_LADDER).ToList();
            BuildBetaGraph_New(realLadders1);

            // Build beta-sheet graph
            List<BetaLadder> realLadders = ladders.Where(l => l.Start0 != l.Start1).Where(l => CheckLadderAndCountHBonds(l) >= HBondSSAConstants.MIN_HBONDS_PER_LADDER).ToList();
            List<BetaStrandInSheet> vertices = new List<BetaStrandInSheet>();
            int sheetCounter = 0;
            foreach (BetaLadder ladder in realLadders)
            {
                Sse strand0 = GetStrand0(ladder).RelabeledCopy("x" + strandCounter++);
                Sse strand1 = GetStrand1(ladder).RelabeledCopy("x" + strandCounter++);
                BetaStrandInSheet lower = new BetaStrandInSheet(strand0, sheetCounter, /*0, */strand0.Start % 2 == 0);
                BetaStrandInSheet upper = new BetaStrandInSheet(strand1, sheetCounter++, /*1, */strand1.Start % 2 != 0);
                LinkStrands(lower, upper, ladder);
                AddVertexAndPossiblyMerge(vertices, lower);
                AddVertexAndPossiblyMerge(vertices, upper);
            }

            // Get beta-sheets(i.e. connected components of the beta-sheet graph)
            List<BetaStrandInSheet> seeds = new List<BetaStrandInSheet>();
            foreach (BetaStrandInSheet v in vertices)
            {
                if (!seeds.Any(u => u.SheetId == v.SheetId))
                {
                    //SetIdToSheet(v, sheetCounter++);
                    seeds.Add(v);
                }
            }
            seeds.Sort((p, q) => Sse.Compare(p.SSE, q.SSE));
            List<List<BetaStrandInSheet>> sheets = new List<List<BetaStrandInSheet>>();
            sheetCounter = 1;
            foreach (BetaStrandInSheet v in seeds)
            {
                List<BetaStrandInSheet> sheet = new List<BetaStrandInSheet>();
                v.DFS(u =>
                {
                    u.SheetId = sheetCounter;
                    u.SSE.SheetId = sheetCounter;
                    IEnumerable<BetaLadder> uLadders = u.DownLadders.Union(u.UpLadders);
                    u.SSE.Type = (u.SSE.End - u.SSE.Start > 0 || uLadders.Any(l => LadderSSEType(l) == SseType.SHEET_TYPE)) ? SseType.SHEET_TYPE : SseType.ISOLATED_BETA_BRIDGE_TYPE;
                    sheet.Add(u);
                });
                sheetCounter++;
                sheets.Add(sheet);
            }

            // // DEBUG: Write out sheet info.
            // Lib.WriteLineDebug("From list of sheets:");
            // foreach (var sheet in sheets) {
            // 	Lib.WriteLineDebug("Sheet {0}:", sheet.First().SheetId);
            // 	foreach (var u in sheet) {
            // 		Lib.WriteLineDebug("Strand {0}({1} {2}-{3}): up: {4}, down: {5}",
            // 			u.SSE.Label,
            // 			u.SSE.ChainID,
            // 			u.SSE.Start,
            // 			u.SSE.End,
            // 			Lib.EnumerateWithCommas(u.UpNeighbours.Select(n => n.SSE.Label)),
            // 			Lib.EnumerateWithCommas(u.DownNeighbours.Select(n => n.SSE.Label)));
            // 	}
            // }

            List<Sse> realStrands = sheets
                .SelectMany(s => s)
                .Select(v => v.SSE)
                .OrderBy(sse => sse)
                .Select((sse, i) => sse.RelabeledCopy(sse.Type.AsString() + i))
                .ToList();

            // Get beta-bulges //TODO: Beta-bulges can be detected directly(as H-bond pattern) instead of this way.
            List<BetaBulge> bulges = new List<BetaBulge>();
            BetaLadder[] bothSideRealLadders = realLadders.Union(realLadders.Select(l => l.Inverted())).ToArray();
            for (int i = 0; i < bothSideRealLadders.Length; i++)
            {
                for (int j = i + 1; j < bothSideRealLadders.Length; j++)
                {
                    BetaBulge bulge = BuildBetaBulge(bothSideRealLadders[i], bothSideRealLadders[j]);
                    if (bulge != null)
                        bulges.Add(bulge);
                }
            }

            int bulgeCounter = 0;

            //List<SSE> bulgeSides=bulges.SelectMany(b=>new SSE[]{GetShortStrand(b),GetLongStrand(b)}).ToList();
            List<Sse> bulgeSides = new List<Sse>();
            //Lib.WriteLineDebug("Detected bulges:");
            foreach (BetaBulge bulge in bulges)
            {
                foreach (Sse b in new Sse[] { GetShortStrand(bulge), GetLongStrand(bulge) })
                {
                    Sse side = b.RelabeledCopy(b.Type.AsString() + bulgeCounter.ToString());
                    Sse[] includingStrands = realStrands.Where(s => s.ChainID == side.ChainID && s.Start <= side.Start && s.End >= side.End).ToArray();
                    if (includingStrands.Length != 1)
                    {
                        throw new InvalidOperationException("One side of a beta-bulge belongs to more than one beta-strand(this should never happen)!");
                    }
                    else
                    {
                        side.SheetId = includingStrands[0].SheetId;
                        includingStrands[0].AddNestedSSE(side);
                        bulgeSides.Add(side);
                    }
                }
                bulgeCounter++;
            }

            List<Sse> resultSSEs = c7TurnsAndWiggles.Union(realStrands)/*.Union(bulgeSides)*/.OrderBy(sse => sse).ToList();
            /*foreach (SSE sse in result) {
                int s = ResidueXBefore(sse.Start, 1);
                int e = ResidueXAfter(sse.End, 1);
                if (s == -1 || e == -1) {
                    sse.AddComment("Edge of chain. ");
                } else if (Enumerable.Range(s, e - s + 1).Any(i => residues[i].GetAtoms().Any(a => a.AltLoc != ' '))){
                    sse.AddComment("Alternative locations. ");
                }
            }*/

            Dictionary<(string, int, int), int> chainStartEnd2Index = new Dictionary<(string, int, int), int>();
            for (int i = 0; i < resultSSEs.Count; i++)
            {
                Sse sse = resultSSEs[i];
                chainStartEnd2Index[(sse.ChainID, sse.Start, sse.End)] = i;
            }
            List<(int, int, int)> edges = new List<(int, int, int)>();
            foreach (var seed in seeds)
            {
                seed.DFS(u =>
                {
                    int vertex1 = chainStartEnd2Index[(u.SSE.ChainID, u.SSE.Start, u.SSE.End)];
                    for (int i = 0; i < u.DownNeighbours.Count; i++)
                    {
                        var v = u.DownNeighbours[i];
                        var ladder = u.DownLadders[i];
                        int vertex2 = chainStartEnd2Index[(v.SSE.ChainID, v.SSE.Start, v.SSE.End)];
                        int ladderType = ladder.Type == BetaLadder.LadderType.Parallel ? 1 : -1;
                        edges.Add((Math.Min(vertex1, vertex2), Math.Max(vertex1, vertex2), ladderType));
                    }
                    /*foreach (var v in u.DownNeighbours) {
                        int vertex2=chainStartEnd2Index[(v.SSE.ChainID, v.SSE.Start, v.SSE.End)];
                        edges.Add((Math.Min(vertex1,vertex2),Math.Max(vertex1,vertex2)));
                    }*/
                });
            }
            edges = edges.Distinct().OrderBy(e => e).ToList();  // a ladder contained a bulge, there would be listed two edges connecting the strands; this removes the duplicite edges

            SecStrAssignment result = new SecStrAssignment(resultSSEs);
            result.Connectivity = edges;
            result.HBonds = hBonds.OrderBy(hb => (hb.Item1.ResidueIndex, hb.Item2.ResidueIndex)).ToList();
            return result;
        }

        public String GetDescription()
        {
            return "hydrogen-bond-based method (similar to DSSP)";
        }

    }
}


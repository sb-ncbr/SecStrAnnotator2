using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class SingleChainGeomSecStrAssigner : ISecStrAssigner {
        public Chain Chain { get; private set; }
        public double RmsdLimit { get; private set; }
        public bool HelicesAllowed { get; private set; }
        public bool SheetsAllowed { get; private set; }

        public SingleChainGeomSecStrAssigner(Chain chain, double rmsdLimit){
            this.Chain=chain;
            this.RmsdLimit=rmsdLimit;
            this.HelicesAllowed=true;
            this.SheetsAllowed=true;
        }

        public SingleChainGeomSecStrAssigner(Chain chain, double rmsdLimit, char[] acceptedSSETypes){
            this.Chain=chain;
            this.RmsdLimit=rmsdLimit;
            this.HelicesAllowed= acceptedSSETypes.Contains (SSE.MIXED_HELIX_TYPE) || acceptedSSETypes.Contains (SSE.HELIX_H_TYPE);
            this.SheetsAllowed=acceptedSSETypes.Contains (SSE.SHEET_TYPE);
        }

        public SecStrAssignment GetSecStrAssignment(){
            int firstAssignedInQuad = 1;
            int lastAssignedInQuad = 2;
            int minimumSSELength = 3;
            List<SSE> result = new List<SSE> ();

            List<Residue> residues = Chain.GetResidues ().Where (r => r.HasCAlpha ()).ToList ();
            bool makingHelix = false;
            bool makingSheet = false;
            char helixType = SSE.MIXED_HELIX_TYPE;
            char sheetType = SSE.SHEET_TYPE;
            int start = 0;
            int end = 0;
            double rmsd = 0;
            double maxRmsd = 0;
            for (int i = 0; i < residues.Count-3; i++) {
                List<Residue> quad = residues.GetRange (i, 4);
                //Console.WriteLine ("at {0}",residues[i].ResSeq);
                if (makingHelix) {
                    if (LibAnnotation.CheckGeometryOf1Unit (quad, helixType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
                        //continue helix
                        end = i + lastAssignedInQuad;
                        maxRmsd = Math.Max (maxRmsd, rmsd);
                    } else {
                        //finish helix
                        if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
                            result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, helixType,null));
                        makingHelix = false;
                    }
                } else if (makingSheet) {
                    if (LibAnnotation.CheckGeometryOf1Unit (quad, sheetType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
                        //continue sheet
                        end = i + lastAssignedInQuad;
                        maxRmsd = Math.Max (maxRmsd, rmsd);
                    } else {
                        //finish sheet
                        //Console.WriteLine ("Ending sheet at {0}", residues [i + lastAssignedInQuad].ResSeq);
                        if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
                            result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, sheetType,null));
                        makingSheet = false;
                    }
                } 
                if (!makingHelix && !makingSheet) {
                    if (HelicesAllowed && LibAnnotation.CheckGeometryOf1Unit (quad, helixType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
                        //start helix
                        //Console.WriteLine ("Starting helix at {0}", residues [i + firstAssignedInQuad].ResSeq);
                        makingHelix = true;
                        start = i + firstAssignedInQuad;
                        end = i + lastAssignedInQuad;
                        maxRmsd = rmsd;
                    } else if (SheetsAllowed && LibAnnotation.CheckGeometryOf1Unit (quad, sheetType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
                        //start sheet
                        //Console.WriteLine ("Starting sheet at {0}", residues [i + firstAssignedInQuad].ResSeq);
                        makingSheet = true;
                        start = i + firstAssignedInQuad;
                        end = i + lastAssignedInQuad;
                        maxRmsd = rmsd;
                    }
                }
            }
            if (makingHelix) {
                //finish helix
                if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
                    result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, helixType,null));
                makingHelix = false;
            }
            if (makingSheet){
                //finish sheet
                if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
                    result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, sheetType,null));
                makingSheet = false;
            }

            if (SheetsAllowed) {
                // Deleting sheets that have no hydrogen bonds
                List<List<Residue>> sheets = Chain.GetResidues (result.Where (sse => sse.Type == sheetType).Select (sse => (sse.Start, sse.End)));
                IEnumerable<Atom> atomsN = sheets.SelectMany (sse => sse.SelectMany (r => r.GetNAmides()));
                IEnumerable<Atom> atomsO = sheets.SelectMany (sse => sse.SelectMany (r => r.GetOCarbs()));
                double maxDistanceForHydrogenBond = 3.5; //TODO Get a rational value for this number (distance N-O). This value is just a guess.
                int minResiduesBetweenHydrogenBond = 2; // number of residues that must be between 2 residues that are forming a beta-sheet stabilizing hydrogen bond

                List<int> deletedSheetsByStart = new List<int> ();
                for (int i = sheets.Count - 1; i >= 0; i--) {
                    int firstResSeq = sheets [i] [0].SeqNumber;
                    int lastResSeq = sheets [i] [sheets [i].Count - 1].SeqNumber;

                    // Version without shortening sheets
                    /*var atomsNHere = sheets [i].SelectMany (r => r.GetAtoms ().Where (a => a.Name == " N  "));
            var atomsOHere = sheets [i].SelectMany (r => r.GetAtoms ().Where (a => a.Name == " O  "));
            if (!atomsNHere.Any (a1 => atomsO.Any (a2 => !(a2.ResSeq >= firstResSeq && a2.ResSeq <= lastResSeq) && (Math.Abs (a2.ResSeq-a1.ResSeq)>minResiduesBetweenHydrogenBond) && LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))
                && !atomsOHere.Any (a1 => atomsN.Any (a2 => !(a2.ResSeq >= firstResSeq && a2.ResSeq <= lastResSeq) && (Math.Abs (a2.ResSeq-a1.ResSeq)>minResiduesBetweenHydrogenBond) && LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))) {
                deletedSheetsByStart.Add (firstResSeq);
                Console.WriteLine ("Deleting sheet at {0}.",firstResSeq);
            }else
                Console.WriteLine ("Keeping sheet at {0}.",firstResSeq);
            result.RemoveAll (sse => sse.Type == sheetType && deletedSheetsByStart.Contains (sse.Start));*/

                    // Version with shortening sheets
                    List<int> hbN = sheets [i].Where (r => r.GetNAmides().Any (a1 => atomsO.Any (a2 => 
                        !(a2.ResidueSeqNumber >= firstResSeq && a2.ResidueSeqNumber <= lastResSeq)
                        && (Math.Abs (a2.ResidueSeqNumber - a1.ResidueSeqNumber) > minResiduesBetweenHydrogenBond)
                        && /*LibProtein.Distance (a1, a2)*/ (a1.Position() - a2.Position()).Size <= maxDistanceForHydrogenBond))).Select (r => r.SeqNumber).ToList ();
                    List<int> hbO = sheets [i].Where (r => r.GetOCarbs().Any (a1 => atomsN.Any (a2 => 
                        !(a2.ResidueSeqNumber >= firstResSeq && a2.ResidueSeqNumber <= lastResSeq)
                        && (Math.Abs (a2.ResidueSeqNumber - a1.ResidueSeqNumber) > minResiduesBetweenHydrogenBond)
                        && /*LibProtein.Distance (a1, a2)*/ (a1.Position() - a2.Position()).Size <= maxDistanceForHydrogenBond))).Select (r => r.SeqNumber).ToList ();
                    List<int> hb = hbN.Union (hbO).ToList ();
                    if (hb.Count == 0) {
                        result.RemoveAll (sse => sse.Start == firstResSeq);
                    } else {
                        int index = result.FindIndex (sse => sse.Start == firstResSeq);
                        SSE old = result [index];
                        //result [index] = new SSE (old.Label, old.ChainID, Math.Max (old.Start, hb.Min ()-1), Math.Min (old.End, hb.Max ()+1), result [index].Type);
                        result [index] = new SSE (old.Label, old.ChainID, hb.Min (), hb.Max (), result [index].Type,null);
                    }
                }
            }

            return new SecStrAssignment(result.Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ())).ToList ());
        }

        public String GetDescription(){
            return "geometrical method with RMSD limit " + RmsdLimit + " A (accepted types: " + (HelicesAllowed?SheetsAllowed?"H, E":"H":SheetsAllowed?"E":"" ) + ")";
        }
    }	
}


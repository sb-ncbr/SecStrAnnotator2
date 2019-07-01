using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class GeomDsspSecStrAssigner : ISecStrAssigner {
        public ISecStrAssigner HelixAssigner { get; private set; }
        public ISecStrAssigner SheetAssigner { get; private set; }

        public GeomDsspSecStrAssigner(ISecStrAssigner helixAssigner, ISecStrAssigner sheetAssigner){
            this.HelixAssigner=helixAssigner;
            this.SheetAssigner=sheetAssigner;
        }

        public GeomDsspSecStrAssigner(Protein protein, IEnumerable<string> chainIDs,/*IEnumerable<Chain> chains,*/ double rmsdLimit, String dsspExecutable, String renumberedPDBFile, String DSSPFile, char[] acceptedSSETypes){
            this.HelixAssigner = new GeomSecStrAssigner(chainIDs.Select(c => protein.GetChain(c)), rmsdLimit, acceptedSSETypes.Intersect (SSE.ALL_HELIX_TYPES).ToArray ());
            this.SheetAssigner = new DsspSecStrAssigner(protein, dsspExecutable, renumberedPDBFile, DSSPFile, chainIDs, acceptedSSETypes.Intersect (SSE.ALL_SHEET_TYPES).ToArray ());
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result = SecStrAssignment.Order (SecStrAssignment.Combine (SheetAssigner.GetSecStrAssignment (), HelixAssigner.GetSecStrAssignment ()));
            //reporting overlapping SSEs
            for (int i = 0; i < result.SSEs.Count - 1; i++) {
                if (result.SSEs [i].End >= result.SSEs [i + 1].Start) {
                    Lib.WriteLineDebug ("Overlapping SSEs: {0} and {1}", result.SSEs [i], result.SSEs [i + 1]);
                }
            }
            // truncating overlapping SSEs
            /*if (result.Select (sse=>sse.ChainID).Distinct ().Count () > 1) 
                throw new NotImplementedException ("This method is not implemented for SSEs coming from multiple chains.");
            for (int i = 0; i < result.Count-1; i++) {
                if (result [i].End >= result [i + 1].Start){
                    //TODO This solution might not be perfect.
                    result[i].End = result[i+1].Start-1;
                    Lib.WriteLineDebug ("Some SSEs were overlapping and therefore they have been truncated.");
                }
            }*/
            return result;
        }

        public String GetDescription(){
            return "mixed method (for helices: " + HelixAssigner.GetDescription () + ", for sheets: " + SheetAssigner.GetDescription ();
        }
    }	
}


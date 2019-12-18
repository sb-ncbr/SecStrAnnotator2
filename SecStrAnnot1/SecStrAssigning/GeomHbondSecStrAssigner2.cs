using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.SSEs;
using protein.Libraries;

namespace protein.SecStrAssigning
{
    public class GeomHbondSecStrAssigner2 : ISecStrAssigner {
        public ISecStrAssigner GeomAssigner { get; private set; }
        public HBondSecStrAssigner2 HBondAssigner2 { get; private set; }
        public const int MIN_HELIX_SUBHELIX_OVERLAP = 1;

        public GeomHbondSecStrAssigner2(Protein p, double rmsdLimit, double hBondEnergyCutoff){
            this.GeomAssigner = new GeomSecStrAssigner(p.GetChains (), rmsdLimit, SSE.ALL_HELIX_TYPES);
            this.HBondAssigner2 = new HBondSecStrAssigner2(p, hBondEnergyCutoff);
        }

        public SecStrAssignment GetSecStrAssignment(){
            HBondAssigner2.DetectSheets = true;
            HBondAssigner2.DetectHelices = false;
            SecStrAssignment sheetAssignment = HBondAssigner2.GetSecStrAssignment ();
            HBondAssigner2.DetectSheets = false;
            HBondAssigner2.DetectHelices = true;
            List<SSE> subhelices = HBondAssigner2.GetSecStrAssignment ().SSEs.OrderBy (sse => sse).ToList ();
            HBondAssigner2.DetectSheets = true;
            List<SSE> bigHelices = GeomAssigner.GetSecStrAssignment ().SSEs.OrderBy (sse => sse).ToList ();

            Func<SSE,SSE,int,bool> doOverlap = (a, b, mo) => a.ChainID == b.ChainID && a.End >= b.Start + mo && b.End >= a.Start + mo;

            List<SSE> processedHelices = new List<SSE> ();
            foreach (var helix in bigHelices) {
                List<SSE> subs = subhelices.Where (sub => doOverlap (helix, sub, MIN_HELIX_SUBHELIX_OVERLAP)).ToList ();
                if (subs.Count > 0) {
                    foreach (var sub in subs) {
                        helix.AddNestedSSE (sub);
                    }
                    var subTypes = subs.Select (sub => sub.Type).Distinct ().ToList();
                    helix.Type = subTypes.Count == 1 ? subTypes[0] : SSE.MIXED_HELIX_TYPE;
                    processedHelices.Add (helix);
                } else {
                    // ignore helix that contains no subhelices (no H-bonds)
                }
            }
            var result = SecStrAssignment.Order(SecStrAssignment.Combine (sheetAssignment, new SecStrAssignment (processedHelices)));

            //reporting overlapping SSEs
            var sses = result.SSEs.OrderBy (x=>x).ToArray ();
            for (int i = 0; i < sses.Length - 1; i++) {
                if (doOverlap(sses[i], sses[i+1], 0))
                    Lib.WriteLineDebug ("Overlapping SSEs: {0} and {1}", sses [i], sses [i + 1]);
            }

            return result;
        }

        public String GetDescription(){
            return "mixed method (for helices: " + GeomAssigner.GetDescription () + " with further division using " + GeomAssigner.GetDescription () + ", for sheets: " + HBondAssigner2.GetDescription ();
        }
    }	
}


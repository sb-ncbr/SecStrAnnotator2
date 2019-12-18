using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.SSEs;

namespace protein.SecStrAssigning
{
    public class GeomSecStrAssigner : ISecStrAssigner {
        public ISecStrAssigner[] Assigners { get; private set; }

        public GeomSecStrAssigner(IEnumerable<Chain> chains, double rmsdLimit){
            Assigners=chains.Select (c=>new SingleChainGeomSecStrAssigner (c,rmsdLimit)).ToArray ();
        }

        public GeomSecStrAssigner(IEnumerable<Chain> chains, double rmsdLimit, char[] acceptedSSETypes){
            Assigners=chains.Select (c=>new SingleChainGeomSecStrAssigner (c,rmsdLimit,acceptedSSETypes)).ToArray ();
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result = new SecStrAssignment (new List<SSE>());
            foreach(ISecStrAssigner assigner in Assigners){
                result = SecStrAssignment.Combine (result, assigner.GetSecStrAssignment ());
            }
            return SecStrAssignment.Order (result);
        }

        public String GetDescription(){
            return Assigners.Count () > 0 ? Assigners [0].GetDescription () : "";
        }
    }	
}


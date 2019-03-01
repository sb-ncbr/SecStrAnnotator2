using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class FilteringSecStrAssigner : ISecStrAssigner {
        /*Returns the same SSEs as its InnerAssigner, but relabels them with this scheme:
            <Prefix><SSE_type><sequentially_assigned_identifier>, e.g. 1tqn_H10 */
        public ISecStrAssigner InnerAssigner { get; private set; }
        public char[] AcceptedSSETypes { get; private set; }
        public string[] ChainIDs{ get; private set; }

        public FilteringSecStrAssigner(ISecStrAssigner innerAssigner, IEnumerable<char> acceptedTypes, IEnumerable<string> chainIDs){
            this.InnerAssigner = innerAssigner;
            this.AcceptedSSETypes=acceptedTypes.ToArray ();
            this.ChainIDs=chainIDs.ToArray ();
        }

        public SecStrAssignment GetSecStrAssignment(){
            /*Lib.Shuffler shuffler;
            List<SSE> result = InnerAssigner
                .GetSecStrAssignment (out connectivity)
                .WhereAndGetShuffler (sse => AcceptedSSETypes.Contains (sse.Type) && ChainIDs.Contains (sse.ChainID), out shuffler)
                .ToList ();
            connectivity=shuffler.UpdateIndices (connectivity).ToList ();*/
            return SecStrAssignment.Filter (InnerAssigner.GetSecStrAssignment (), sse => AcceptedSSETypes.Contains (sse.Type) && ChainIDs.Contains (sse.ChainID));
        }

        public String GetDescription(){
            return InnerAssigner.GetDescription ();
        }
    }	
}


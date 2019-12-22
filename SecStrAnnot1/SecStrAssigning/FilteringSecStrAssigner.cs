using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Sses;

namespace protein.SecStrAssigning
{
    public class FilteringSecStrAssigner : ISecStrAssigner
    {
        /*Returns the same SSEs as its InnerAssigner, but relabels them with this scheme:
            <Prefix><SSE_type><sequentially_assigned_identifier>, e.g. 1tqn_H10 */
        public ISecStrAssigner InnerAssigner { get; private set; }
        public SseType[] AcceptedSSETypes { get; private set; }
        public string[] ChainIDs { get; private set; }

        public FilteringSecStrAssigner(ISecStrAssigner innerAssigner, IEnumerable<SseType> acceptedTypes, IEnumerable<string> chainIDs)
        {
            this.InnerAssigner = innerAssigner;
            this.AcceptedSSETypes = acceptedTypes.ToArray();
            this.ChainIDs = chainIDs.ToArray();
        }

        public SecStrAssignment GetSecStrAssignment()
        {
            return SecStrAssignment.Filter(InnerAssigner.GetSecStrAssignment(), sse => AcceptedSSETypes.Contains(sse.Type) && ChainIDs.Contains(sse.ChainID));
        }

        public String GetDescription()
        {
            return InnerAssigner.GetDescription();
        }
    }
}


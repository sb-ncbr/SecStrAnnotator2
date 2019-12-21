using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Sses;

namespace protein.SecStrAssigning
{
    public class RelabellingSecStrAssigner : ISecStrAssigner {
        /*Returns the same SSEs as its InnerAssigner, but relabels them with this scheme:
            <Prefix><SSE_type><sequentially_assigned_identifier>, e.g. 1tqn_H10 */
        public ISecStrAssigner InnerAssigner { get; private set; }
        public String Prefix { get; private set; }
        public bool SetAllLabelsToNull{ get; private set; }

        public RelabellingSecStrAssigner(ISecStrAssigner innerAssigner, String prefix, bool setAllLabelToNull){
            this.InnerAssigner = innerAssigner;
            this.Prefix = prefix;
            this.SetAllLabelsToNull=setAllLabelToNull;
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
            result.SSEs=result.SSEs
                .Select ((sse, i) => sse.RelabeledCopy (SetAllLabelsToNull ? null: Prefix + sse.Type.AsString() + i.ToString ()))
                .ToList ();
            return result;
        }

        public String GetDescription(){
            return InnerAssigner.GetDescription ();
        }
    }	
}


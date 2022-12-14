using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Sses;

namespace protein.SecStrAssigning
{
    public class AuthFieldsAddingSecStrAssigner : ISecStrAssigner {
        public ISecStrAssigner InnerAssigner { get; private set; }
        public Protein Protein { get; private set; }
        public AuthFieldsAddingSecStrAssigner(ISecStrAssigner innerAssigner, Protein protein){
            this.InnerAssigner = innerAssigner;
            this.Protein = protein;
        }
        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment ass = InnerAssigner.GetSecStrAssignment();
            foreach (Sse sse in ass.SSEs){
                sse.AddAuthFields(Protein);
            }
            return ass;
        }
        public String GetDescription(){
            return InnerAssigner.GetDescription ();
        }
    
    }	
}


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Libraries;
using protein.SSEs;

namespace protein.SecStrAssigning
{
    public class JoiningSecStrAssigner : ISecStrAssigner
    {
        public ISecStrAssigner InnerAssigner { get; private set; }
        public Protein Protein { get; private set; }
        public double JoiningRmsdLimit { get; private set; }
        public Func<SSEType, SSEType, SSEType?> JoiningTypeCombining { get; private set; }

        public JoiningSecStrAssigner(ISecStrAssigner innerAssigner, Protein protein, double joiningRmsdLimit, Func<SSEType, SSEType, SSEType?> joiningTypeCombining)
        {
            this.InnerAssigner = innerAssigner;
            this.Protein = protein;
            this.JoiningRmsdLimit = joiningRmsdLimit;
            this.JoiningTypeCombining = joiningTypeCombining;
        }


        public SecStrAssignment GetSecStrAssignment()
        {
            //TODO connectivity
            throw new NotImplementedException();
            SecStrAssignment unjoined = InnerAssigner.GetSecStrAssignment();
            return new SecStrAssignment(LibAnnotation.JoinSSEs_GeomVersion(Protein, unjoined.SSEs, JoiningRmsdLimit, JoiningTypeCombining));
        }

        public String GetDescription()
        {
            return InnerAssigner.GetDescription() + " and then joining using geometric criterion with RMSD limit " + JoiningRmsdLimit.ToString("0.000") + " A";
        }
    }
}


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class OutputtingSecStrAssigner : ISecStrAssigner {
        public ISecStrAssigner InnerAssigner { get; private set; }
        public String OutputFile { get; private set; }

        public OutputtingSecStrAssigner(ISecStrAssigner innerAssigner, String outputFile){
            this.InnerAssigner = innerAssigner;
            this.OutputFile = outputFile;
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
            LibAnnotation.WriteAnnotationFile (OutputFile, result.SSEs, "Helix info obtained by " + InnerAssigner.GetDescription() + ".");
            //TODO write out connectivity and hbonds
            return result;
        }

        public String GetDescription(){
            return InnerAssigner.GetDescription ();
        }
    }	
}


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class OutputtingSecStrAssigner_Json : ISecStrAssigner {
        public ISecStrAssigner InnerAssigner { get; private set; }
        public String EntryName{ get; private set; }
        public String OutputFile { get; private set; }

        public OutputtingSecStrAssigner_Json(ISecStrAssigner innerAssigner, String outputFile, String entryName){
            if (entryName==null) throw new ArgumentNullException ();
            this.InnerAssigner = innerAssigner;
            this.EntryName=entryName;
            this.OutputFile = outputFile;
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
            LibAnnotation.WriteAnnotationFile_Json (OutputFile, 
                EntryName, 
                result.SSEs,
                null, 
                result.Connectivity,
                result.HBonds,
                "Helix info obtained by " + InnerAssigner.GetDescription() + "."
            );
            //TODO write out connectivity and hbonds
            return result;
        }

        public String GetDescription(){
            return InnerAssigner.GetDescription ();
        }
    }	
}


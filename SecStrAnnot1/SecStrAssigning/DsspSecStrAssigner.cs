using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class DsspSecStrAssigner : ISecStrAssigner {
        public String DSSPExecutable { get; private set; }
        public String PDBFile { get; private set; }
        public String DSSPFile { get; private set; }
        public string[] ChainIDs { get; private set; }
        public char[] AcceptedSSETypes { get; private set; }

        public DsspSecStrAssigner(String dsspExecutable, String PDBFile, String DSSPFile, IEnumerable<string> chainIDs, char[] acceptedSSETypes){
            this.DSSPExecutable = dsspExecutable;
            this.PDBFile = PDBFile;
            this.DSSPFile = DSSPFile;
            this.ChainIDs=chainIDs.ToArray ();
            this.AcceptedSSETypes = acceptedSSETypes;
        }

        public SecStrAssignment GetSecStrAssignment(){
            Lib.WriteInColor (ConsoleColor.Yellow, "Running DSSP:\n");
            if (!Lib.RunDSSP (DSSPExecutable, PDBFile, DSSPFile))
                throw new  SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.");
            
            try { 
                return new SecStrAssignment(LibAnnotation.ReadSSEsFromDSSP (DSSPFile, AcceptedSSETypes)
                    .Where (x => ChainIDs.Contains (x.ChainID))
                    .Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ()))
                    .ToList ()); 
            }
            catch (Exception e) { throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e); }
        }

        public String GetDescription(){
            return "DSSP method for " + Path.GetFileName (PDBFile) + " (accepted types: " + Lib.EnumerateWithCommas (AcceptedSSETypes) + ")";
        }
    }	
}


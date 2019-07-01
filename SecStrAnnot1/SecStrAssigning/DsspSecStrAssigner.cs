using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    // public class DsspSecStrAssigner : ISecStrAssigner {
    //     public String DSSPExecutable { get; private set; }
    //     public String PDBFile { get; private set; }
    //     public String DSSPFile { get; private set; }
    //     public string[] ChainIDs { get; private set; }
    //     public char[] AcceptedSSETypes { get; private set; }

    //     public DsspSecStrAssigner(String dsspExecutable, String PDBFile, String DSSPFile, IEnumerable<string> chainIDs, char[] acceptedSSETypes){
    //         this.DSSPExecutable = dsspExecutable;
    //         this.PDBFile = PDBFile;
    //         this.DSSPFile = DSSPFile;
    //         this.ChainIDs=chainIDs.ToArray ();
    //         this.AcceptedSSETypes = acceptedSSETypes;
    //     }

    //     public SecStrAssignment GetSecStrAssignment(){
    //         Lib.WriteInColor (ConsoleColor.Yellow, "Running DSSP:\n");
    //         if (!Lib.RunDSSP (DSSPExecutable, PDBFile, DSSPFile))
    //             throw new  SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.");
            
    //         try { 
    //             return new SecStrAssignment(LibAnnotation.ReadSSEsFromDSSP (DSSPFile, AcceptedSSETypes)
    //                 .Where (x => ChainIDs.Contains (x.ChainID))
    //                 .Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ()))
    //                 .ToList ()); 
    //         }
    //         catch (Exception e) { throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e); }
    //     }

    //     public String GetDescription(){
    //         return "DSSP method for " + Path.GetFileName (PDBFile) + " (accepted types: " + Lib.EnumerateWithCommas (AcceptedSSETypes) + ")";
    //     }
    // }	


    public class DsspSecStrAssigner : ISecStrAssigner {
        public Protein Protein { get; private set; }
        public String DSSPExecutable { get; private set; }
        public String RenumberedPDBFile { get; private set; }  // CIF file with label_ values in auth_ fields
        public String DSSPFile { get; private set; }
        public string[] ChainIDs { get; private set; }
        public char[] AcceptedSSETypes { get; private set; }

        public DsspSecStrAssigner(Protein protein, String dsspExecutable, String renumberedPDBFile, String dsspFile, IEnumerable<string> chainIDs, char[] acceptedSSETypes){
            this.Protein = protein;
            this.DSSPExecutable = dsspExecutable;
            this.RenumberedPDBFile = renumberedPDBFile;
            this.DSSPFile = dsspFile;
            this.ChainIDs = chainIDs.ToArray ();
            this.AcceptedSSETypes = acceptedSSETypes;
        }

        public SecStrAssignment GetSecStrAssignment(){
            SecStrAssignment result;
            this.Protein.SaveCif(RenumberedPDBFile, fillAuthFieldsWithLabelValues: true, fabulateOccupancyAndBFactor: true);
            Lib.WriteInColor (ConsoleColor.Yellow, "Running DSSP:\n");
            if (!Lib.RunDSSP (DSSPExecutable, RenumberedPDBFile, DSSPFile)){
                throw new  SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.");
            }            
            try { 
                IEnumerable<SSE> sses = LibAnnotation.ReadSSEsFromDSSP (DSSPFile, AcceptedSSETypes)
                    .Where (x => ChainIDs.Contains (x.ChainID))
                    .Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ()));
                result = new SecStrAssignment(sses); 
                // TODO add support for beta-connectivity
            }
            catch (Exception e) { throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e); }
            // File.Delete(RenumberedPDBFile);
            // File.Delete(DSSPFile);
            return result;
        }

        public String GetDescription(){
            return "DSSP method for " + Path.GetFileName (RenumberedPDBFile) + " (accepted types: " + Lib.EnumerateWithCommas (AcceptedSSETypes) + ")";
        }
    }	
}


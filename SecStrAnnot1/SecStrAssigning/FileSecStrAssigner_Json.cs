using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class FileSecStrAssigner_Json : ISecStrAssigner {
        public String SsaFile { get; private set; }
        public String EntryName{ get; private set; }
        public string[] ChainIDs { get; private set; }

        /* If entryName==null, it will take the first entry in the file. */
        public FileSecStrAssigner_Json(String ssaFile, String entryName, IEnumerable<string> chainIDs){
            this.SsaFile=ssaFile;
            this.EntryName=entryName;
            this.ChainIDs=chainIDs.ToArray ();
        }

        public SecStrAssignment GetSecStrAssignment(){
            try {
                String dump;
                List<(int, int, int)> connectivity;
                List<(String, int, int)> mergeable;
                Lib.Shuffler shuffler;
                SecStrAssignment result = new SecStrAssignment (LibAnnotation.ReadAnnotationFile_Json (SsaFile, EntryName, out dump, out connectivity, out mergeable, true).WhereAndGetShuffler (x => ChainIDs.Contains (x.ChainID), out shuffler).ToList ());
                connectivity=shuffler.UpdateIndices (connectivity).ToList ();
                result.Connectivity=connectivity;
                result.MergeableSSEs=mergeable;
                return result;
            } catch (FormatException e){
                Lib.WriteError ("Invalid format of the template annotation file:\n    " + e.Message );
                throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
            } catch (IOException e) {
                Lib.WriteError ("Could not open \"" + SsaFile + "\".");
                throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
            }
        }

        public String GetDescription(){
            return "reading from JSON file " + Path.GetFileName (SsaFile) + " (chainIDs: " + ChainIDs.EnumerateWithCommas () + ")";
        }
    }	
}


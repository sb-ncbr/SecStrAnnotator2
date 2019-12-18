using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using protein.Libraries;

namespace protein.SecStrAssigning
{
    public class FileSecStrAssigner : ISecStrAssigner {
        public String SsaFile { get; private set; }
        public string[] ChainIDs { get; private set; }

        public FileSecStrAssigner(String ssaFile,  IEnumerable<string> chainIDs){
            this.SsaFile=ssaFile;
            this.ChainIDs=chainIDs.ToArray ();
        }

        public SecStrAssignment GetSecStrAssignment(){
            try {
                return new SecStrAssignment(LibAnnotation.ReadAnnotationFile (SsaFile).Where (x => ChainIDs.Contains (x.ChainID)).ToList ());
            } catch (IOException e) {
                Lib.WriteError ("Could not open \"" + SsaFile + "\".");
                throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
            }
        }

        public String GetDescription(){
            return "reading from file " + Path.GetFileName (SsaFile) + " (chainIDs: " + ChainIDs.EnumerateWithCommas () + ")";
        }
    }	
}


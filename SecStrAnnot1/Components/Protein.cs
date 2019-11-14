using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace protein.Components
{
    public class Protein
    {
		public const bool IGNORE_ALTERNATIVE_LOCATIONS = true;
		private static bool alternativeLocationWarningPrinted = false;

        private Dictionary<string,Chain> chains;

		private void InitializeThis(IEnumerable<Atom> atoms){
			chains = new Dictionary<string, Chain>();

			if (atoms != null) {
				foreach (Atom a in atoms) {
					if (a.AltLoc != ' ' && a.AltLoc != 'A' && IGNORE_ALTERNATIVE_LOCATIONS) {
						if (!alternativeLocationWarningPrinted) {
							Lib.WriteWarning ("Alternative locations found. Ignoring alternative locations except for ' ' and 'A'.");
							alternativeLocationWarningPrinted = true;
						}				
					} else {
						this.AddAtom (a);
					}
				}
			}
		}

        private Protein()
        {
			InitializeThis (null);
        }
		
		public Protein(IEnumerable<Atom> atoms)
		{	
			InitializeThis (atoms);
		}

		public Protein(StreamReader reader)
		{
			List<Residue> modifiedResidues;
			IEnumerable<Atom> atoms = ReadAtomsFromPDB (reader);
			InitializeThis (atoms);
		}

		public Protein(StreamReader reader, string[] chainIDs, IEnumerable<(int, int)> resSeqRanges)
		{
			List<Residue> modifiedResidues;
			IEnumerable<Atom> atoms = ReadAtomsFromPDB (reader);
			if (chainIDs != null) {
				atoms = atoms.Where (a => chainIDs.Contains (a.ChainID));
			}
			if (resSeqRanges != null) {
				atoms = atoms.Where (a => resSeqRanges.InRanges (a.ResSeq));
			}
			InitializeThis (atoms);
		}

		private static List<Atom> ReadAtomsFromPDB(StreamReader reader){
			List<Atom> atoms = new List<Atom> ();
			String line;
			while ((line = reader.ReadLine ()) != null) {
				if (line.Length >= 6
				    && (line.Substring (0, 6) == "ATOM  " || line.Substring (0, 6) == "HETATM")) 
				{
					Atom atom = new Atom (line);
					atoms.Add (atom);
				} 
				else if (line.Length >= 6 && line.Substring (0, 6) == "ENDMDL") {
					break;
				}
			}
			return atoms;
		}

		/** Print protein to a file.*/
		public void Save(String outputFile)
		{
			try {
				StreamWriter w = new StreamWriter (outputFile);
				foreach (Atom a in this.GetAtoms()) {
					w.WriteLine (a.ToString ());
				}
				w.Close ();
			} catch (IOException e) {
				Console.Error.WriteLine ("Error: Could not open \"" + outputFile + "\" for writing.");
				throw new IOException ("Protein.Save: Could not open \"" + outputFile + "\" for writing.",e);
			}
		}

        /**
         * Returns a protein containing only those residues which have a C alpha atom.
         */
		public Protein KeepOnlyNormalResidues(bool doPrintWarningForHet)
		{
			Protein result = new Protein ();
			List<Residue> hetResWithCA = new List<Residue> ();
			foreach (Residue res in this.GetResidues()) {
				if (res.HasCAlpha ()) {
					result.AddResidue (res);
					if (res.IsHet ()) {
						hetResWithCA.Add (res);
					}
				}
			}
			if (doPrintWarningForHet && hetResWithCA.Count > 0) {
				Lib.WriteWarning ("Found hetero residues with C-alpha. They will be treated as normal residues: \n{0}", hetResWithCA.Select(r => r.ToString(false)).EnumerateWithCommas());
			}
			return result;
		}

        public SortedSet<Chain> GetChains()
        {
            return new SortedSet<Chain>(chains.Values);
        }
        
        public SortedSet<Residue> GetResidues()
        {
            SortedSet<Residue> result = new SortedSet<Residue>();
            foreach (Chain chain in this.GetChains())
            {
                result.UnionWith(chain.GetResidues());
            }
            return result;
        }
        
        public SortedSet<Atom> GetAtoms()
        {
            SortedSet<Atom> result = new SortedSet<Atom>();
            foreach (Residue residue in this.GetResidues())
            {
                result.UnionWith(residue.GetAtoms());
            }
            return result;
        }

		public Chain GetChain(string chainID)
		{
			return chains[chainID];
		}

		public bool HasChain(string chainID)
		{
			return chains.ContainsKey(chainID);
		}

        public void AddChain(string chainID, Chain newChain)
        {
            if (!chains.ContainsKey(chainID))
            {
                chains.Add(chainID, newChain);
            }
            else
            {
                throw new InvalidOperationException("Chain " + chainID + " already present.");
            }
        }

        public void AddResidue(Residue residue)
        {
            foreach (Atom atom in residue.GetAtoms())
            {
                AddAtom(atom);
            }
        }

		public void AddAtom(Atom atom)
		{
			if (atom.ICode != ' ') {
				if (!Setting.IgnoreInsertions) {
					Lib.WriteErrorAndExit ("Protein contains residues with insertion codes (not supported in current version). To ignore these residues, use option -i or --ignoreinsertions (loaded structure will not correspond fully to the input file!).");
				} else if (!Setting.IgnoreInsertionsWarningThrown) {
					Lib.WriteWarning ("Ignoring residues with insertion codes!");
					Setting.IgnoreInsertionsWarningThrown = true;
				}
			} else {
				if (!chains.ContainsKey (atom.ChainID)) {
					chains.Add (atom.ChainID, new Chain (atom.ChainID));
				}
				chains [atom.ChainID].AddAtom (atom);
			}
		}

		/*public static Protein FromCifModel(Model model){
			throw new NotImplementedException ();
		}*/
    }
}

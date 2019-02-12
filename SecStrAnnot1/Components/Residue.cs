using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace protein.Components
{
    public class Residue : IComparable<Residue>
    {
        private String name;
        public String Name
        {
            get
            {
                return name;
            }
            private set
            {
                name = value;
                // if (value.Length == 3)
                // {
                //     name = value;
                // }
                // else
                // {
                //     throw new ArgumentException("Residue name must consist of exactly 3 characters, not " + value);
                // }
            }
        }
		public char ShortName {
			get {
				char result;
				return namesLongToShort.TryGetValue (this.Name, out result) ? result : UNKNOWN_RESIDUE_1_LETTER;
			}
		}
        public string ChainID { get; private set; }
        public int ResSeq { get; private set; }
		//public bool IsModified { get; set; }
        private SortedSet<Atom> atoms;

		private static String[] normalNames = { "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PYL","PRO","GLN","ARG","SER","THR","SEC","VAL","TRP","TYR" };

		private static Dictionary<String,char> namesLongToShort = new Dictionary<String,char> {
			{ "ALA",'A' }, { "CYS",'C' }, { "ASP",'D' }, { "GLU",'E' }, { "PHE",'F' }, { "GLY",'G' }, { "HIS",'H' }, { "ILE",'I' },
			{ "LYS",'K' }, { "LEU",'L' }, { "MET",'M' }, { "ASN",'N' }, { "PYL",'O' }, { "PRO",'P' }, { "GLN",'Q' }, { "ARG",'R' },
			{ "SER",'S' }, { "THR",'T' }, { "SEC",'U' }, { "VAL",'V' }, { "TRP",'W' }, { "TYR",'Y' }
		};
		public const char UNKNOWN_RESIDUE_1_LETTER = '?';

		public Residue(String name, string chainID, int resSeq)
		{
			//Program.logger.Debug("new Residue(" + name + " " + chainID + " " + resSeq);
			Name = name;
			ChainID = chainID;
			ResSeq = resSeq;
			atoms = new SortedSet<Atom>();
		}

		public Residue(String name, string chainID, int resSeq, IEnumerable<Atom> atoms)
			: this(name,chainID,resSeq)
		{
			foreach (Atom a in atoms) {
				this.AddAtom (a);
			}
		}

		public Residue WithoutHydrogens(){
			var newAtoms = this.GetAtoms ().Where (a => !a.IsHydrogen);
			Residue result = new Residue (this.Name, this.ChainID, this.ResSeq, newAtoms);
			return result;
		}

        public override bool Equals(object obj)
        {
            if (obj is Residue)
            {
                Residue res = obj as Residue;
                return Name == res.Name && ChainID == res.ChainID && ResSeq == res.ResSeq;
            }
            else
            {
                return false;
            }
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode() + ChainID.GetHashCode() + ResSeq.GetHashCode();
        }

        public int CompareTo(Residue other)
        {
            int result = this.ChainID.CompareTo(other.ChainID);
            if (result != 0) return result;
            result = this.ResSeq - other.ResSeq;
            if (result != 0) return result;
            return this.Name.CompareTo(other.Name);
        }

        public SortedSet<Atom> GetAtoms() { return new SortedSet<Atom>(atoms); }

		public static SortedSet<Atom> GetAtomsOfAll(ICollection<Residue> residues)
		{
			SortedSet<Atom> result = new SortedSet<Atom>();
			foreach (Residue residue in residues)
			{
				result.UnionWith(residue.GetAtoms());
			}
			return result;
		}

        public void AddAtom(Atom atom)
        {
            if (atom.ResName != name)
				throw new ArgumentException("Residue.AddAtom(): residue name " + name + " expected " + atom.ResName + " found (in " + this.ToString (false) + ")");
            if (atom.ChainID != ChainID)
				throw new ArgumentException("Residue.AddAtom(): chain " + ChainID + " expected " + atom.ChainID + " found (in " + this.ToString (false) + ")");
            if (atom.ResSeq != ResSeq)
				throw new ArgumentException("Residue.AddAtom(): residue number " + ResSeq + " expected " + atom.ResSeq + " found (in " + this.ToString (false) + ")");
            atoms.Add(atom);
        }

		public bool IsProline{ get { return this.Name == "PRO"; } }

		public bool HasCAlpha(){
			return GetAtoms ().Any (a => a.IsCAlpha);
		}

		public bool IsHet(){
			try {
				return atoms.First().IsHetAtm;
			} catch (InvalidOperationException){ //no atoms
				return false;
			}
		}

		public override string ToString ()
			{
				return ToString (true);
			}

		public string ToString(bool shortString){
			if (shortString) {
				return this.ResSeq + this.Name;
			} else {
				return "Residue " + this.ChainID + " " + this.ResSeq + " " + this.Name;
			}
		}

        
    }
}

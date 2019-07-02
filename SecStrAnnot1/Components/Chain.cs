using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace protein.Components
{
    public class Chain : IComparable<Chain>
    {
        public string ID { get; private set; }
        private ChainBuilder chainBuilder;
        public int FirstResidueNumber { get { return chainBuilder.FirstResidueNumber; } }
        public int LastResidueNumber { get { return chainBuilder.LastResidueNumber; } }
        public int Count { get { return chainBuilder.Count; } }

		public List<Residue> GetResidues() { return chainBuilder.GetResidues().ToList (); }

		public List<List<Residue>> GetResidues(IEnumerable<(int, int)> fromTos) { 
			List<(int, int, int)> ranges = fromTos.Select((t, i) => (t.Item1, t.Item2, i)).OrderBy (t => t.Item1).ToList (); //Adding to each tuple Item3 with its index in the original fromTos, than ordering.
			int nextRange = 0;
			List<int> activeRanges = new List<int> ();
			List<Residue> residues = this.GetResidues ();
			List<List<Residue>> result = fromTos.Select (x=>new List<Residue>()).ToList ();
			foreach (Residue r in residues) {
				while (nextRange < ranges.Count && ranges [nextRange].Item1 <= r.ResSeq) {
					activeRanges.Add (nextRange);
					nextRange++;
				}
				activeRanges.RemoveAll (i => ranges [i].Item2 < r.ResSeq);
				foreach (int i in activeRanges) {
					result [ranges[i].Item3].Add (r); // ranges[i].Item3 = original index of ranges[i] in fromTos
				}
			}
			return result;
		}

		public SortedSet<Atom> GetAtoms() { return chainBuilder.GetAtoms(); }
		public static SortedSet<Residue> GetResiduesOfAll(ICollection<Chain> chains)
		{
			SortedSet<Residue> result = new SortedSet<Residue>();
			foreach (Chain chain in chains)
			{
				result.UnionWith(chain.GetResidues());
			}
			return result;
		}
		public static SortedSet<Atom> GetAtomsOfAll(ICollection<Chain> chains)
		{
			SortedSet<Atom> result = new SortedSet<Atom>();
			foreach (Chain chain in chains)
			{
				result.UnionWith(chain.GetAtoms());
			}
			return result;
		}

        public Residue GetResidue(int resSeq) {
            if (GetResidues().Count == 0 || resSeq < FirstResidueNumber || resSeq > LastResidueNumber)
            {
                throw new IndexOutOfRangeException();
            }
            return chainBuilder.GetResidue(resSeq);
        }
		//public List<Residue> ListResidues

        public void AddResidue(int resSeq, Residue newResidue)
        {
            if (Count == 0 || resSeq == FirstResidueNumber - 1 || resSeq == LastResidueNumber + 1)
            {
                chainBuilder.AddResidue(newResidue);
            }
            else
            {
                throw new InvalidOperationException();
            }
        }

		public Chain(string id)
        {
			this.ID = id;
            this.chainBuilder = new ChainBuilder();
        }
        
        public void AddAtom(Atom atom) { chainBuilder.AddAtom(atom); }

        public int CompareTo(Chain other) { 
			int res1 = this.ID.CompareTo (other.ID);
			return  res1 != 0 ? res1 : this.chainBuilder.CompareTo (other.chainBuilder);
		}

        public class LengthComparer : IComparer<Chain>
        {
            public int Compare(Chain c1, Chain c2) {
                int result = c1.Count - c2.Count;
                if (result !=0) return result;
                else return c1.CompareTo(c2);
            }
        }
    }
}

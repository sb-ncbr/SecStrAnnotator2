using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace protein
{
    class ChainBuilder:IComparable<ChainBuilder>
    {
        private Dictionary<int, Residue> residues;
        public int FirstResidueNumber { get; private set; }
        public int LastResidueNumber { get; private set; }
        public int Count { get { return residues.Count; } }
        
        public ChainBuilder()
        {
            residues = new Dictionary<int, Residue>();
            FirstResidueNumber = 0;
            LastResidueNumber = 0;
        }

        public ChainBuilder(ICollection<Atom> atoms)
            : this()
        {
            foreach (Atom atom in atoms)
            {
                AddAtom(atom);
            }
        }

        public SortedSet<Residue> GetResidues()
        {
            return new SortedSet<Residue>(residues.Values);
        }
        public SortedSet<Atom> GetAtoms()
        {
			return new SortedSet<Atom> (GetResidues ().SelectMany (r => r.GetAtoms ()));
        }

        public Residue GetResidue(int resSeq)
        {
            if (residues.ContainsKey(resSeq))
            {
                return residues[resSeq];
            }
            else
            {
                return null;
            }
        }
        
        public void AddResidue(Residue newResidue)
        {
            if (residues.Count == 0)
            {
                FirstResidueNumber = newResidue.ResSeq;
                LastResidueNumber = newResidue.ResSeq;
            }
            FirstResidueNumber = Math.Min(FirstResidueNumber, newResidue.ResSeq);
            LastResidueNumber = Math.Max(LastResidueNumber, newResidue.ResSeq);
            if (!residues.ContainsKey(newResidue.ResSeq))
            {
                residues.Add(newResidue.ResSeq, newResidue);
            }
        }

        public void AddAtom(Atom atom)
        {
            if (!residues.ContainsKey(atom.ResSeq))
            {
                AddResidue(new Residue(atom.ResName, atom.ChainID, atom.ResSeq));
            }
            residues[atom.ResSeq].AddAtom(atom);
        }

        public SortedSet<Chain> GetChains()
        {
            //Program.logger.OpenBlockD("ChainBuilder.GetChains()");
            SortedSet<Chain> result = new SortedSet<Chain>();
            Chain currentChain = null;

            if (Count == 0)
            {
                //Program.logger.Debug("Empty ChainBuilder");
                return result;
            }

            //Program.logger.Debug("Nonempty ChainBuilder, FirstResidueNumber = " + FirstResidueNumber + ", LastResidueNumber = " + LastResidueNumber);
            for (int i = FirstResidueNumber; i <= LastResidueNumber; i++)
            {
                Residue currentResidue = GetResidue(i);
                if (currentChain != null)
                {
                    if (currentResidue != null)
                    {
                        //Program.logger.Debug("Residue " + i + " present, added to chain " + result.Count);
                        currentChain.AddResidue(i, currentResidue);
                    }
                    else
                    {
                        //Program.logger.Debug("Residue " + i + " not present, finishing chain " + result.Count);
                        result.Add(currentChain);
                        currentChain = null;
                    }
                }
                else
                {
                    if (currentResidue != null)
                    {
                        //Program.logger.Debug("Residue " + i + " present, added to new chain " + result.Count);
                        currentChain = new Chain(' ');
                        currentChain.AddResidue(i, currentResidue);
                    }
                    else
                    {
                        //Program.logger.Debug("Residue " + i + " not present");
                    }
                }
            }
            if (currentChain != null)
            {
                //Program.logger.Debug("End of the ChainBuilder, finishing chain " + result.Count);
                result.Add(currentChain);
            }
            //Program.logger.CloseBlockD("ChainBuilder.GetChains()");
            return result;
        }

        public int CompareTo(ChainBuilder other)
        {
            if (this.Count == 0)
            {
                if (other.Count == 0) return 0;
                else return -1;
            }
            else
            {
                if (other.Count == 0) return 1;
                else return this.residues[this.FirstResidueNumber].CompareTo(other.residues[other.FirstResidueNumber]);
            }
        }
    }
}

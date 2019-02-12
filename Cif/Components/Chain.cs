using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Tables;

namespace Cif.Components
{
    public struct Chain
    {
        public Model Model { get; private set; }
        public int ChainIndex { get; private set; }

        public Chain(Model model, int chainIndex){
            this.Model = model;
            this.ChainIndex = chainIndex;
        }

        // own properties
        public string ID => Model.Chains.Id[ChainIndex];
                        
        
        public IEnumerable<int> GetResidueIndices(){
            int start = Model.Chains.ResidueStartIndex[ChainIndex];
            int end = Model.Chains.ResidueEndIndex[ChainIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Residue> GetResidues(){
            Model model = this.Model;
            return GetResidueIndices().Select(ri => new Residue(model, ri));
        }
		public List<List<Residue>> GetResidues(IEnumerable<Tuple<int,int>> fromTos) { 
			List<Tuple<int,int,int>> ranges = fromTos.Select ((t,i)=>new Tuple<int,int,int> (t.Item1,t.Item2,i)).OrderBy (t => t.Item1).ToList (); //Adding to each tuple Item3 with its index in the original fromTos, than ordering.
			int nextRange = 0;
			List<int> activeRanges = new List<int> ();
			List<List<Residue>> result = fromTos.Select (x => new List<Residue>()).ToList ();
            Model model = this.Model;
			foreach (int ri in GetResidueIndices()) {
				while (nextRange < ranges.Count && ranges [nextRange].Item1 <= Model.Residues.SeqNumber[ri]) {
					activeRanges.Add (nextRange);
					nextRange++;
				}
				activeRanges.RemoveAll (i => ranges [i].Item2 < model.Residues.SeqNumber[ri]);
				foreach (int i in activeRanges) {
					result [ranges[i].Item3].Add (new Residue(Model, ri)); // ranges[i].Item3 = original index of ranges[i] in fromTos
				}
			}
			return result;
		}
        
        public IEnumerable<int> GetAtomIndices(){
            int start = Model.Chains.AtomStartIndex[ChainIndex];
            int end = Model.Chains.AtomEndIndex[ChainIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Atom> GetAtoms(){
            Model model = this.Model;
            return GetAtomIndices().Select(ai => new Atom(model, ai));
        }

    }
}
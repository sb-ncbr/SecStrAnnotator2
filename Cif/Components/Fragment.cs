using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Tables;

namespace Cif.Components
{
    public struct Fragment
    {
        public Model Model { get; private set; }
        public int FragmentIndex { get; private set; }

        public Fragment(Model model, int chainIndex){
            this.Model = model;
            this.FragmentIndex = chainIndex;
        }

        // own properties
        //--

        // indirect properties
        public string ChainId => Model.Chains.Id[Model.Fragments.ChainIndex[FragmentIndex]];
        public string ChainAuthId => Model.Chains.AuthId[Model.Fragments.ChainIndex[FragmentIndex]];
        public string EntityId => Model.Entities.Id[Model.Chains.EntityIndex[FragmentIndex]];
                        
        
        public IEnumerable<int> GetResidueIndices(){
            int start = Model.Fragments.ResidueStartIndex[FragmentIndex];
            int end = Model.Fragments.ResidueEndIndex[FragmentIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Residue> GetResidues(){
            Model model = this.Model;
            return GetResidueIndices().Select(ri => new Residue(model, ri));
        }
        
        public IEnumerable<int> GetAtomIndices(){
            int start = Model.Fragments.AtomStartIndex[FragmentIndex];
            int end = Model.Fragments.AtomEndIndex[FragmentIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Atom> GetAtoms(){
            Model model = this.Model;
            return GetAtomIndices().Select(ai => new Atom(model, ai));
        }

        public override string ToString(){
            int firstResIndex = Model.Fragments.ResidueStartIndex[FragmentIndex];
            int lastResIndex = Model.Fragments.ResidueEndIndex[FragmentIndex]-1;
            return $"Fragment {ChainId} {Model.Residues.SeqNumber[firstResIndex]}:{Model.Residues.SeqNumber[lastResIndex]}";
        }
    }
}
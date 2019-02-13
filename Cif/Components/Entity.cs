using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Tables;

namespace Cif.Components
{
    public struct Entity
    {
        public Model Model { get; private set; }
        public int EntityIndex { get; private set; }

        public Entity(Model model, int entityIndex){
            this.Model = model;
            this.EntityIndex = entityIndex;
        }

        // own properties
        public string Id => Model.Entities.Id[EntityIndex];
                        
        
        public IEnumerable<int> GetChainIndices(){
            int start = Model.Entities.ChainStartIndex[EntityIndex];
            int end = Model.Entities.ChainEndIndex[EntityIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Chain> GetChains(){
            Model model = this.Model;
            return GetChainIndices().Select(ci => new Chain(model, ci));
        }
        public bool HasChain(string chainID){
            foreach (int ci in GetChainIndices())
            {
                if (Model.Chains.Id[ci] == chainID){
                    return true;
                }
            }
            return false;
        }
        public Chain GetChain(string chainID){
            foreach (int ci in GetChainIndices())
            {
                if (Model.Chains.Id[ci] == chainID){
                    return new Chain(Model, ci);
                }
            }
            throw new KeyNotFoundException($"Chain with ID {chainID} was not found.");
        }

        public IEnumerable<int> GetResidueIndices(){
            int start = Model.Entities.ResidueStartIndex[EntityIndex];
            int end = Model.Entities.ResidueEndIndex[EntityIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Residue> GetResidues(){
            Model model = this.Model;
            return GetResidueIndices().Select(ri => new Residue(model, ri));
        }
        
        public IEnumerable<int> GetAtomIndices(){
            int start = Model.Entities.AtomStartIndex[EntityIndex];
            int end = Model.Entities.AtomEndIndex[EntityIndex];
            return Enumerable.Range(start, end-start);
        }
        public IEnumerable<Atom> GetAtoms(){
            Model model = this.Model;
            return GetAtomIndices().Select(ai => new Atom(model, ai));
        }

        public override string ToString(){
            return $"Entity {Id}";
        }
    }
}
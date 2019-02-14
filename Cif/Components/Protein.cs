using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Tables;
using Cif.Libraries;
using System.Text;

namespace Cif.Components
{
    public struct Protein
    {
        public Model Model { get; private set; }

        public Protein (Model model){
            this.Model = model;
        }


        public IEnumerable<int> GetEntityIndices(){
            return Enumerable.Range(0, Model.Entities.Count);
        }
        public IEnumerable<Entity> GetEntities(){
            Model model = this.Model;
            return GetEntityIndices().Select(ei => new Entity(model, ei));
        }
        
        public IEnumerable<int> GetChainIndices(){
            return Enumerable.Range(0, Model.Chains.Count);
        }
        public IEnumerable<Chain> GetChains(){
            Model model = this.Model;
            return GetChainIndices().Select(ci => new Chain(model, ci));
        }
        public bool HasChain(string chainID){
            for (int ci = 0; ci < Model.Chains.Count; ci++)
            {
                if (Model.Chains.Id[ci] == chainID){
                    return true;
                }
            }
            return false;
        }
        public Chain GetChain(string chainID){
            for (int ci = 0; ci < Model.Chains.Count; ci++)
            {
                if (Model.Chains.Id[ci] == chainID){
                    return new Chain(Model, ci);
                }
            }
            throw new KeyNotFoundException($"Chain with ID {chainID} was not found.");
        }

        public IEnumerable<int> GetResidueIndices(){
            return Enumerable.Range(0, Model.Residues.Count);
        }
        public IEnumerable<Residue> GetResidues(){
            Model model = this.Model;
            return GetResidueIndices().Select(ri => new Residue(model, ri));
        }

        public IEnumerable<int> GetAtomIndices(){
            return Enumerable.Range(0, Model.Atoms.Count);
        }
        public IEnumerable<Atom> GetAtoms(){
            Model model = this.Model;
            return GetAtomIndices().Select(ai => new Atom(model, ai));
        }

        /**
         * Returns a protein containing only those residues which have a C alpha atom.
         */
		public Protein KeepOnlyNormalResidues(bool doPrintWarningForHet){
			List<Residue> hetResWithCA = new List<Residue> ();

            ModelBuilder builder = new ModelBuilder();
            foreach (Entity entity in this.GetEntities()){
                builder.StartEntity(entity.Id);
                foreach (Chain chain in entity.GetChains()) {
                    builder.StartChain(chain.Id, chain.AuthId);
                    foreach (Residue residue in chain.GetResidues()) {
                        if (residue.HasCAlpha()) {
                            builder.StartResidue(residue.SeqNumber, residue.Compound);
                            foreach (Atom atom in residue.GetAtoms()) {
                                builder.AddAtom(atom.Id, atom.AtomInfo());
                            }
                            if (residue.GetAtoms().First().IsHetatm){
                                hetResWithCA.Add(residue);
                            }
                        }
                    }
                }
            }
			if (doPrintWarningForHet && hetResWithCA.Count > 0) {
				Lib.WriteWarning ("Found hetero residues with C-alpha. They will be treated as normal residues: \n{0}", string.Join(", ", hetResWithCA));
			}
            Protein result = new Protein(builder.GetModel(this.Model.ModelNumber));
            Lib.WriteLineDebug($"KeepOnlyNormalResidues(): {result.Model.Residues.Count}");
            return result;
            //TODO implement this somehow!
        }

        public string ToLongString() {
            StringBuilder b = new StringBuilder();
            b.AppendLine("Protein");
            foreach (var chain in this.GetChains()) {
                b.AppendLine($"    {chain}");
                foreach (var residue in chain.GetResidues()) {
                    b.AppendLine($"        {residue}");
                    foreach (var atom in residue.GetAtoms()) {
                        b.AppendLine($"            {atom}");
                    }
                }
            }
            return b.ToString();
        }
    }
}
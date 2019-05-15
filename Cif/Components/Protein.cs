using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Tables;
using Cif.Libraries;
using System.Text;
using System.IO;

namespace Cif.Components
{
    public struct Protein
    {
        public Model Model { get; private set; }

        public Protein (ModelBuilder modelBuilder, int modelNumber = Model.DEFAULT_MODEL_NUM) : this(modelBuilder.GetModel(modelNumber)) {}
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
                            builder.StartResidue(residue.ResidueInfo());
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
            Protein result = new Protein(builder, this.Model.ModelNumber);
            // Lib.WriteLineDebug($"KeepOnlyNormalResidues(): {result.Model.Residues.Count}");
            return result;
            //TODO implement this somehow!
        }

        /**
         * Returns a protein containing only one alternative location (the first occurring in mmCIF file) of each atom.
         */
		public Protein KeepOnlyOneAlternativeLocation(){
            var seenAtoms = new HashSet<ValueTuple<string,int,string>>(); // each tuple is (chain ID, residue index, atom name)
            int removedAtoms = 0;
            var keptAltLocs = new HashSet<string>();
            var removedAltLocs = new HashSet<string>();
            ModelBuilder builder = new ModelBuilder();
            foreach (Entity entity in this.GetEntities()){
                builder.StartEntity(entity.Id);
                foreach (Chain chain in entity.GetChains()) {
                    builder.StartChain(chain.Id, chain.AuthId);
                    foreach (Residue residue in chain.GetResidues()) {
                        builder.StartResidue(residue.ResidueInfo());
                        foreach (Atom atom in residue.GetAtoms()) {
                            var atomIdentification = (chain.Id, residue.SeqNumber, atom.Name);
                            bool isUnseen = seenAtoms.Add(atomIdentification);
                            if (isUnseen){
                                builder.AddAtom(atom.Id, atom.AtomInfo());
                                keptAltLocs.Add(atom.AltLoc);
                            } else {
                                removedAtoms++; // atom with same chainId, resi and name has already been added
                                removedAltLocs.Add(atom.AltLoc);
                            }
                        }
                    }
                }
            }
			if (removedAtoms > 0) {
                string kept = string.Join(", ", keptAltLocs.OrderBy(a=>a).Select(a=> $"'{a}'"));
                string removed = string.Join(", ", removedAltLocs.OrderBy(a=>a).Select(a=> $"'{a}'"));
				Lib.WriteWarning ($"Found some atoms with identical chain ID, residue number, and atom name -> removing redundant atoms. Kept alternative locations: {kept}; removed: {removed}");
			}
            Protein result = new Protein(builder, this.Model.ModelNumber);
            return result;
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

		/** Print protein to a file.*/
		public void SaveCif(String outputFile, string dataName = "structure")
		{
			try {
                string categoryString = this.Model.ToCifCategoryString();
				StreamWriter w = new StreamWriter (outputFile);
                w.WriteLine("data_" + dataName);
                w.WriteLine("#");
                w.Write(categoryString);
				w.Close ();
			} catch (IOException e) {
				Console.Error.WriteLine ("Error: Could not open \"" + outputFile + "\" for writing.");
				throw new IOException ("Protein.SaveCif: Could not open \"" + outputFile + "\" for writing.", e);
			}
		}

		/** Print label_ to auth_numbering conversion table to a file.*/
		public void SaveLabel2AuthTable(String outputFile)
		{
            using(StreamWriter w = new StreamWriter(outputFile)){
                w.WriteLine($"#{ChainTable.ID_COLUMN}\t{ResidueTable.SEQ_NUMBER_COLUMN}\t{ChainTable.AUTH_ID_COLUMN}\t{ResidueTable.AUTH_SEQ_NUMBER_COLUMN}\t{ResidueTable.AUTH_INSERTION_CODE_COLUMN}");
                foreach (Residue residue in this.GetResidues())
                {
                    w.WriteLine($"{residue.ChainId}\t{residue.SeqNumber}\t{residue.ChainAuthId}\t{residue.AuthSeqNumber}\t{residue.AuthInsertionCode}");
                }           
            }
		}
    }
}
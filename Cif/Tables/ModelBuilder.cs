using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Cif.Libraries;

namespace Cif.Tables
{
    public class ModelBuilder
    {
        const string DEFAULT_ENTITY_ID = "1";
        const string DEFAULT_CHAIN_ID = "1";
        const int DEFAULT_RESIDUE_SEQ_ID = 1;
        const string DEFAULT_RESIDUE_COMPOUND = "XXX";
        const string DEFAULT_ATOM_ID = "1";

        private class EntityTab {
            public int count = 0;
            // down
            public List<int> atomStartIndex = new List<int>();
            public List<int> residueStartIndex = new List<int>();
            public List<int> chainStartIndex = new List<int>();
            // own properties
            public List<string> Id = new List<string>();
            
            public void AddRow(int chain, int residue, int atom, string id){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    Id[count-1] = id;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    residueStartIndex.Add(residue);
                    chainStartIndex.Add(chain);
                    Id.Add(id);
                }
            }
        }

        private class ChainTab {
            public int count = 0;
            // down
            public List<int> atomStartIndex = new List<int>();
            public List<int> residueStartIndex = new List<int>();
            // own properties
            public List<string> Id = new List<string>();
            public List<string> AuthId = new List<string>();

            public void AddRow(int entity, int residue, int atom, string id, string authId){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    Id[count-1] = id;
                    AuthId[count-1] = authId;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    residueStartIndex.Add(residue);
                    Id.Add(id);
                    AuthId.Add(authId);
                }
            }
        }

        private class ResidueTab {
            public int count = 0;
            // down
            public List<int> atomStartIndex = new List<int>();
            // own properties
            public List<int> SeqNumber = new List<int>();
            public List<string> Compound = new List<string>();

            public void AddRow(int entity, int chain, int atom, int seqNumber, string compound){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    SeqNumber[count-1] = seqNumber;
                    Compound[count-1] = compound;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    SeqNumber.Add(seqNumber);
                    Compound.Add(compound);
                }
            }
        }

        private class AtomTab {
            public int count = 0;
            // down
            // own properties
            public List<string> Id = new List<string>();
            public List<AtomInfo> Info = new List<AtomInfo>();

            public void AddRow(int entity, int chain, int residue, string id, AtomInfo info){
                //add row
                count++;
                Id.Add(id);
                Info.Add(info);
            }
        }

        private EntityTab entities;
        private ChainTab chains;
        private ResidueTab residues;
        private AtomTab atoms;

        public ModelBuilder() {
            this.entities = new EntityTab();
            this.chains = new ChainTab();
            this.residues = new ResidueTab();
            this.atoms = new AtomTab();
            StartEntity();
        }

        public void StartEntity(string id){
            entities.AddRow(chains.count, residues.count, atoms.count, id);
            StartChain();
        }
        public void StartEntity(){
            string id = entities.count == 0 ? DEFAULT_ENTITY_ID : SuccessorIdentifier(entities.Id[entities.count-1]);
            StartEntity(id);
        }

        public void StartChain(string id, string authId){
            chains.AddRow(entities.count-1, residues.count, atoms.count, id, authId);
            StartResidue();
        }
        public void StartChain(){
            string id = chains.count == 0 ? DEFAULT_CHAIN_ID : SuccessorIdentifier(chains.Id[chains.count-1]);
            StartChain(id, id);
        }

        public void StartResidue(int seqNumber, string compound){
            residues.AddRow(entities.count-1, chains.count-1, atoms.count, seqNumber, compound);
        }
        public void StartResidue(){
            int seqNumber = residues.count == 0 ? DEFAULT_RESIDUE_SEQ_ID : residues.SeqNumber[residues.count-1] + 1;
            StartResidue(seqNumber, DEFAULT_RESIDUE_COMPOUND);
        }

        public void AddAtom(string id, AtomInfo atomInfo){
            atoms.AddRow(entities.count-1, chains.count-1, residues.count-1, id, atomInfo);
        }
        public void AddAtom(AtomInfo atomInfo){
            string id = atoms.count == 0 ? DEFAULT_ATOM_ID : SuccessorIdentifier(atoms.Id[atoms.count-1]);
            AddAtom(id, atomInfo);
        }

        public Model GetModel(int modelNumber = Model.DEFAULT_MODEL_NUM){

            // count valid residues, chains..., i.e. those which contain at least one atom, residue...
            bool lastResidueEmpty = atoms.count == residues.atomStartIndex[residues.count-1];
            int validResidues = lastResidueEmpty ? residues.count-1 : residues.count;
            bool lastChainEmpty = validResidues == chains.residueStartIndex[chains.count-1];
            int validChains = lastChainEmpty ? chains.count-1 : chains.count;
            bool lastEntityEmpty = validChains == entities.chainStartIndex[entities.count-1];
            int validEntities = lastEntityEmpty ? entities.count-1 : entities.count;
            
            return new Model(modelNumber, 
                             residues.atomStartIndex.Take(validResidues).Append(atoms.count).ToArray(),
                             chains.residueStartIndex.Take(validChains).Append(validResidues).ToArray(), 
                             entities.chainStartIndex.Take(validEntities).Append(validChains).ToArray(),
                             entities.Id.Take(validEntities).ToArray(),
                             chains.Id.Take(validChains).ToArray(), chains.AuthId.Take(validChains).ToArray(),
                             residues.SeqNumber.Take(validResidues).ToArray(), residues.Compound.Take(validResidues).ToArray(),
                             atoms.Id.ToArray(), atoms.Info.ToArray()
                             );
        }


        /** Divides a string into a prefix and number, e.g. "ABC" => ("ABC",0), "OMG123" => ("OMG",123), "123" => ("",123) */
        private static (string, int) PrefixAndNumber(string str){
            string prefix = str;
            int number = 0;
            for (int i = str.Length - 1; i >= 0 ; i--)
            {
                string suffix = str.Substring(i);
                int o;
                if (int.TryParse(suffix, out o)){
                    prefix = str.Substring(0, i);
                    number = o;
                } else {
                    break;
                }
            }
            return (prefix, number);
        }

        private static string SuccessorIdentifier(string id){
            (string prefix, int number) = PrefixAndNumber(id);
            return prefix + (number + 1).ToString();
        }

        public static void Test(){
            ModelBuilder b = new ModelBuilder();
            b.StartChain("A", "A");
            b.StartResidue(1, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(2, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(3, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(5, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(6, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartChain("B", "B");
            b.StartResidue(7, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(8, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(3, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(5, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue(6, "ALA");
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
                    
            Cif.Components.Protein p = new Cif.Components.Protein(b.GetModel());

            Console.WriteLine(p.ToLongString());
        }
    }
}
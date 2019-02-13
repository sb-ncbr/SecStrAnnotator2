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
            public List<int> fragmentStartIndex = new List<int>();
            public List<int> chainStartIndex = new List<int>();
            // up
            // own properties
            public List<string> Id = new List<string>();
            
            public void AddRow(int chain, int fragment, int residue, int atom, string id){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    Id[count-1] = id;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    residueStartIndex.Add(residue);
                    fragmentStartIndex.Add(fragment);
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
            public List<int> fragmentStartIndex = new List<int>();
            // up
            public List<int> entityIndex = new List<int>();
            // own properties
            public List<string> Id = new List<string>();
            public List<string> AuthId = new List<string>();

            public void AddRow(int entity, int fragment, int residue, int atom, string id, string authId){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    Id[count-1] = id;
                    AuthId[count-1] = authId;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    residueStartIndex.Add(residue);
                    fragmentStartIndex.Add(fragment);
                    entityIndex.Add(entity);
                    Id.Add(id);
                    AuthId.Add(authId);
                }
            }
        }

        private class FragmentTab {
            public int count = 0;
            // down
            public List<int> atomStartIndex = new List<int>();
            public List<int> residueStartIndex = new List<int>();
            // up
            public List<int> chainIndex = new List<int>();
            public List<int> entityIndex = new List<int>();
            // own properties

            public void AddRow(int entity, int chain, int residue, int atom){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    // do nothing
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    residueStartIndex.Add(residue);
                    chainIndex.Add(chain);
                    entityIndex.Add(entity);
                }
            }
        }

        private class ResidueTab {
            public int count = 0;
            // down
            public List<int> atomStartIndex = new List<int>();
            // up
            public List<int> fragmentIndex = new List<int>();
            public List<int> chainIndex = new List<int>();
            public List<int> entityIndex = new List<int>();
            // own properties
            public List<int> SeqNumber = new List<int>();
            public List<string> Compound = new List<string>();

            public void AddRow(int entity, int chain, int fragment, int atom, int seqNumber, string compound){
                if (count > 0 && atom == atomStartIndex[count-1]){
                    //update row
                    SeqNumber[count-1] = seqNumber;
                    Compound[count-1] = compound;
                } else {
                    //add row
                    count++;
                    atomStartIndex.Add(atom);
                    fragmentIndex.Add(fragment);
                    chainIndex.Add(chain);
                    entityIndex.Add(entity);
                    SeqNumber.Add(seqNumber);
                    Compound.Add(compound);
                }
            }
        }

        private class AtomTab {
            public int count = 0;
            // down
            // up
            public List<int> residueIndex = new List<int>();
            public List<int> fragmentIndex = new List<int>();
            public List<int> chainIndex = new List<int>();
            public List<int> entityIndex = new List<int>();
            // own properties
            public List<string> Id = new List<string>();
            public List<AtomInfo> Info = new List<AtomInfo>();
            // public List<string> Name = new List<string>();
            // public List<string> Element = new List<string>();
            // public List<string> AltLoc = new List<string>();
            // public List<bool> IsHetatm = new List<bool>();
            // public List<double> X = new List<double>();
            // public List<double> Y = new List<double>();
            // public List<double> Z = new List<double>();

            public void AddRow(int entity, int chain, int fragment, int residue, string id, AtomInfo info){
                //add row
                count++;
                residueIndex.Add(residue);
                fragmentIndex.Add(fragment);
                chainIndex.Add(chain);
                entityIndex.Add(entity);
                Id.Add(id);
                Info.Add(info);
            }
            // public void AddRow(int entity, int chain, int fragment, int residue, string id, AtomInfo info){
            //     //add row
            //     count++;
            //     residueIndex.Add(residue);
            //     fragmentIndex.Add(fragment);
            //     chainIndex.Add(chain);
            //     entityIndex.Add(entity);
            //     Id.Add(id);
            //     Name.Add(info.Name);
            //     Element.Add(info.Element);
            //     AltLoc.Add(info.AltLoc);
            //     IsHetatm.Add(info.IsHetatm);
            //     X.Add(info.X);
            //     Y.Add(info.Y);
            //     Z.Add(info.Z);
            // }
        }

        private EntityTab entities;
        private ChainTab chains;
        private FragmentTab fragments;
        private ResidueTab residues;
        private AtomTab atoms;

        public ModelBuilder() {
            this.entities = new EntityTab();
            this.chains = new ChainTab();
            this.fragments = new FragmentTab();
            this.residues = new ResidueTab();
            this.atoms = new AtomTab();
            StartEntity();
        }

        public void StartEntity(string id){
            entities.AddRow(chains.count, fragments.count, residues.count, atoms.count, id);
            StartChain();
        }
        public void StartEntity(){
            string id = entities.count == 0 ? DEFAULT_ENTITY_ID : SuccessorIdentifier(entities.Id[entities.count-1]);
            StartEntity(id);
        }

        public void StartChain(string id, string authId){
            chains.AddRow(entities.count-1, fragments.count, residues.count, atoms.count, id, authId);
            StartFragment();
            // StartResidue();
        }
        public void StartChain(){
            string id = chains.count == 0 ? DEFAULT_CHAIN_ID : SuccessorIdentifier(chains.Id[chains.count-1]);
            StartChain(id, id);
        }

        public void StartFragment(){
            fragments.AddRow(entities.count-1, chains.count-1, residues.count, atoms.count);
            StartResidue();
        }

        public void StartResidue(int seqNumber, string compound){
            residues.AddRow(entities.count-1, chains.count-1, fragments.count-1, atoms.count, seqNumber, compound);
        }
        public void StartResidue(){
            int seqNumber = residues.count == 0 ? DEFAULT_RESIDUE_SEQ_ID : residues.SeqNumber[residues.count-1] + 1;
            StartResidue(seqNumber, DEFAULT_RESIDUE_COMPOUND);
        }

        public void AddAtom(string id, AtomInfo atomInfo){
            atoms.AddRow(entities.count-1, chains.count-1, fragments.count-1, residues.count-1, id, atomInfo);
        }
        public void AddAtom(AtomInfo atomInfo){
            string id = atoms.count == 0 ? DEFAULT_ATOM_ID : SuccessorIdentifier(atoms.Id[atoms.count-1]);
            AddAtom(id, atomInfo);
        }

        public Model GetModel() => GetModel(Model.DEFAULT_MODEL_NUM);
        public Model GetModel(int modelNumber){

            // count valid residues, fragments..., i.e. those which contain at least one atom, residue...
            bool lastResidueEmpty = atoms.count == residues.atomStartIndex[residues.count-1];
            int validResidues = lastResidueEmpty ? residues.count-1 : residues.count;
            bool lastFragmentEmpty = validResidues == fragments.residueStartIndex[fragments.count-1];
            int validFragments = lastFragmentEmpty ? fragments.count-1 : fragments.count;
            bool lastChainEmpty = validFragments == chains.fragmentStartIndex[chains.count-1];
            int validChains = lastChainEmpty ? chains.count-1 : chains.count;
            bool lastEntityEmpty = validChains == entities.chainStartIndex[entities.count-1];
            int validEntities = lastEntityEmpty ? entities.count-1 : entities.count;
            
            return new Model(modelNumber, 
                             residues.atomStartIndex.Take(validResidues).Append(atoms.count).ToArray(), 
                             fragments.residueStartIndex.Take(validFragments).Append(validResidues).ToArray(), 
                             chains.fragmentStartIndex.Take(validChains).Append(validFragments).ToArray(), 
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
            b.StartChain("C", "E");
            b.AddAtom(new AtomInfo("CA","C",".",false,0,1,2));
            b.StartChain();
            b.StartResidue();
            b.AddAtom(new AtomInfo("N","N",".",false,0,1,2));
            b.StartResidue();
            b.AddAtom(new AtomInfo("O","O",".",false,0,1,2));
            b.AddAtom(new AtomInfo("C","C",".",false,0,1,2));
                    
            Cif.Components.Protein p = new Cif.Components.Protein(b.GetModel());

            Console.WriteLine(p.ToLongString());
            // foreach (var chain in p.GetChains()) {
            //     Console.WriteLine($"{chain}");
            //     foreach (var residue in chain.GetResidues()) {
            //         Console.WriteLine($"    {residue}");
            //         foreach (var atom in residue.GetAtoms()) {
            //             Console.WriteLine($"        {atom}");
            //         }
            //     }
            // }
            // Console.WriteLine("Entities");
            // Console.WriteLine(string.Join("  ", b.entities.Id));
            // Console.WriteLine(string.Join("  ", b.entities.atomStartIndex));

            // Console.WriteLine("Chains");
            // Console.WriteLine(string.Join("  ", b.chains.Id));
            // Console.WriteLine(string.Join("  ", b.chains.atomStartIndex));

            // Console.WriteLine("Fragments");
            // Console.WriteLine(string.Join("  ", b.fragments.atomStartIndex));

            // Console.WriteLine("Residues");
            // Console.WriteLine(string.Join("  ", b.residues.SeqNumber));
            // Console.WriteLine(string.Join("  ", b.residues.Compound));
            // Console.WriteLine(string.Join("  ", b.residues.atomStartIndex));

            // Console.WriteLine("Atoms");
            // Console.WriteLine(string.Join("  ", b.atoms.Id));
            // Console.WriteLine(string.Join("  ", b.atoms.Info)); 
            // Console.WriteLine(string.Join("  ", b.atoms.residueIndex)); 
            // Console.WriteLine(string.Join("  ", b.atoms.fragmentIndex)); 
            // Console.WriteLine(string.Join("  ", b.atoms.chainIndex)); 
            // Console.WriteLine(string.Join("  ", b.atoms.entityIndex)); 
        }
    }
}
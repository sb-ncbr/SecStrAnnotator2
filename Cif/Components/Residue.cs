using System.Collections.Generic;
using System.Linq;
using Cif.Tables;

namespace Cif.Components {
    public struct Residue {
        public Model Model { get; private set; }
        public int ResidueIndex { get; private set; }

        public Residue(Model model, int residueIndex) {
            this.Model = model;
            this.ResidueIndex = residueIndex;
        }

        // own properties
        public int SeqNumber => Model.Residues.SeqNumber[ResidueIndex];
        public string Compound => Model.Residues.Compound[ResidueIndex];
        public int? AuthSeqNumber => Model.Residues.AuthSeqNumber[ResidueIndex];
        public string AuthInsertionCode => Model.Residues.AuthInsertionCode[ResidueIndex];
        public string AuthCompound => Model.Residues.AuthCompound[ResidueIndex];

        // indirect properties
        public string ChainId => Model.Chains.Id[Model.Residues.ChainIndex[ResidueIndex]];
        public string ChainAuthId => Model.Chains.AuthId[Model.Residues.ChainIndex[ResidueIndex]];
        public string EntityId => Model.Entities.Id[Model.Residues.EntityIndex[ResidueIndex]];


        public ResidueInfo ResidueInfo() => new ResidueInfo(SeqNumber, Compound, AuthSeqNumber, AuthInsertionCode, AuthCompound);

        public override string ToString() => ToString(true);
        public string ToString(bool shortString) {
            if (shortString) {
                return $"{this.SeqNumber} {this.Compound}";
            } else {
                return $"Residue {this.ChainId} {this.SeqNumber} {this.Compound}";
            }
        }


        public IEnumerable<int> GetAtomIndices() {
            int start = Model.Residues.AtomStartIndex[ResidueIndex];
            int end = Model.Residues.AtomEndIndex[ResidueIndex];
            return Enumerable.Range(start, end - start);
        }
        public IEnumerable<Atom> GetAtoms() {
            Model model = this.Model;
            return GetAtomIndices().Select(ai => new Atom(model, ai));
        }

        private Atom? GetAtom(string element, string name) {
            for (int ai = Model.Residues.AtomStartIndex[ResidueIndex]; ai < Model.Residues.AtomEndIndex[ResidueIndex]; ai++) {
                if (Model.Atoms.Element[ai] == element && Model.Atoms.Name[ai] == name) {
                    return new Atom(Model, ai);
                }
            }
            return null;
        }
        public bool HasCAlpha() => GetCAlpha() != null;
        public Atom? GetCAlpha() => GetAtom("C", "CA");
        public Atom? GetNAmide() => GetAtom("N", "N");
        public Atom? GetHAmide() => GetAtom("H", "H");
        public Atom? GetCCarb() => GetAtom("C", "C");
        public Atom? GetOCarb() => GetAtom("O", "O");

        private IEnumerable<Atom> GetAtoms(string element, string name) {
            for (int ai = Model.Residues.AtomStartIndex[ResidueIndex]; ai < Model.Residues.AtomEndIndex[ResidueIndex]; ai++) {
                if (Model.Atoms.Element[ai] == element && Model.Atoms.Name[ai] == name) {
                    yield return new Atom(Model, ai);
                }
            }
        }
        public IEnumerable<Atom> GetCAlphas() => GetAtoms("C", "CA");
        public IEnumerable<Atom> GetNAmides() => GetAtoms("N", "N");
        public IEnumerable<Atom> GetHAmides() => GetAtoms("H", "H");
        public IEnumerable<Atom> GetCCarbs() => GetAtoms("C", "C");
        public IEnumerable<Atom> GetOCarbs() => GetAtoms("O", "O");


        private static readonly Dictionary<string, char> namesLongToShort = new Dictionary<string, char> {
            { "ALA",'A' }, { "CYS",'C' }, { "ASP",'D' }, { "GLU",'E' }, { "PHE",'F' }, { "GLY",'G' }, { "HIS",'H' }, { "ILE",'I' },
            { "LYS",'K' }, { "LEU",'L' }, { "MET",'M' }, { "ASN",'N' }, { "PYL",'O' }, { "PRO",'P' }, { "GLN",'Q' }, { "ARG",'R' },
            { "SER",'S' }, { "THR",'T' }, { "SEC",'U' }, { "VAL",'V' }, { "TRP",'W' }, { "TYR",'Y' }
        };
        public static readonly string[] STANDARD_RESIDUE_3LETTER_NAMES = namesLongToShort.Keys.ToArray();

        public const char UNKNOWN_RESIDUE_1_LETTER = '?';
        public char ShortName => namesLongToShort.GetValueOrDefault(this.Compound, UNKNOWN_RESIDUE_1_LETTER);
        public string ShortOrParenthesizedName {
            get {
                char shortName = this.ShortName;
                if (shortName == UNKNOWN_RESIDUE_1_LETTER){
                    return $"({this.Compound})";
                } else {
                    return shortName.ToString();
                }
            }
        }

        public bool IsProline() {
            return this.Compound == "PRO";
        }

    }
}
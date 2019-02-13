using System.Globalization;
using Cif.Tables;

namespace Cif.Components
{
    public struct Atom
    {
        public Model Model { get; private set; }
        public int AtomIndex { get; private set; }

        public Atom(Model model, int atomIndex){
            this.Model = model;
            this.AtomIndex = atomIndex;
        }

        // own properties
        public string Serial => Model.Atoms.Id[AtomIndex];
        public string Name => Model.Atoms.Name[AtomIndex];
        public string Element => Model.Atoms.Element[AtomIndex];
        public bool IsHetatm => Model.Atoms.IsHetatm[AtomIndex];        
        public string AltLoc => Model.Atoms.AltLoc[AtomIndex];
        public double X => Model.Atoms.X[AtomIndex];
        public double Y => Model.Atoms.Y[AtomIndex];
        public double Z => Model.Atoms.Z[AtomIndex];

        // indirect properties
        public int ResSeq => Model.Residues.SeqNumber[Model.Atoms.ResidueIndex[AtomIndex]];
        public string ResName => Model.Residues.Compound[Model.Atoms.ResidueIndex[AtomIndex]];
        public string ChainID => Model.Chains.Id[Model.Atoms.ChainIndex[AtomIndex]];


        public protein.Vector Position() => new protein.Vector (X, Y, Z);
        
        static readonly CultureInfo CULTURE_INFO = new CultureInfo("en-US");
        public override string ToString()
        {
			return (IsHetatm?"HETATM":"ATOM  ")
                + Serial.ToString().PadLeft(5)
                + " "
                + Name.PadRight(4)
                + " "
                + AltLoc
                + " "
                + ResName
                + " "
                + ChainID
                + " "
                + ResSeq.ToString().PadLeft(4)
                + "  "
                + X.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + Y.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + Z.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + " "
                + Element.PadLeft(2);

        }

    }
}
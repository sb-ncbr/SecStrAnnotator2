using System.Globalization;
using Cif.Tables;

namespace Cif.Components
{
    public struct Atom
    {
        public Model Model { get; private set; }
        public int AtomIndex { get; private set; }

        public Atom(Model model, int atomIndex)
        {
            this.Model = model;
            this.AtomIndex = atomIndex;
        }

        // own properties
        public string Id => Model.Atoms.Id[AtomIndex];
        public string Name => Model.Atoms.Name[AtomIndex];
        public string Element => Model.Atoms.Element[AtomIndex];
        public bool IsHetatm => Model.Atoms.IsHetatm[AtomIndex];
        public string AltLoc => Model.Atoms.AltLoc[AtomIndex];
        public double X => Model.Atoms.X[AtomIndex];
        public double Y => Model.Atoms.Y[AtomIndex];
        public double Z => Model.Atoms.Z[AtomIndex];

        // indirect properties
        public int ResidueSeqNumber => Model.Residues.SeqNumber[Model.Atoms.ResidueIndex[AtomIndex]];
        public string ResidueCompound => Model.Residues.Compound[Model.Atoms.ResidueIndex[AtomIndex]];
        public string ChainId => Model.Chains.Id[Model.Atoms.ChainIndex[AtomIndex]];
        public string ChainAuthId => Model.Chains.Id[Model.Atoms.ChainIndex[AtomIndex]];
        public string EntityId => Model.Entities.Id[Model.Atoms.EntityIndex[AtomIndex]];


        public protein.Geometry.Point Position() => new protein.Geometry.Point(X, Y, Z);
        public AtomInfo AtomInfo() => new AtomInfo(Name, Element, AltLoc, IsHetatm, X, Y, Z);

        static readonly CultureInfo CULTURE_INFO = new CultureInfo("en-US");
        public override string ToString()
        {
            return (IsHetatm ? "HETATM" : "ATOM  ")
                + Id.ToString().PadLeft(5)
                + " "
                + Name.PadRight(4)
                + " "
                + AltLoc
                + " "
                + ResidueCompound
                + " "
                + ChainId
                + " "
                + ResidueSeqNumber.ToString().PadLeft(4)
                + "  "
                + X.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + Y.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + Z.ToString("0.000", CULTURE_INFO).PadLeft(8)
                + " "
                + Element.PadLeft(2);

        }

        public const string NAME_C_ALPHA = "CA";
        public const string ELEMENT_C = "C";
        public const string NAME_H_AMIDE = "H";
        public const string ELEMENT_H = "H";

    }
}
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
        public string AltLoc => Model.Atoms.AltLoc[AtomIndex];
        public double X => Model.Atoms.X[AtomIndex];
        public double Y => Model.Atoms.Y[AtomIndex];
        public double Z => Model.Atoms.Z[AtomIndex];

        // indirect properties
        public int ResSeq => Model.Residues.SeqNumber[Model.Atoms.ResidueIndex[AtomIndex]];


        public protein.Vector Position() => new protein.Vector (X, Y, Z);
        

    }
}
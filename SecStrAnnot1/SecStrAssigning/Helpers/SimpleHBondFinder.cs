using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

using protein.Geometry;

namespace protein.SecStrAssigning.Helpers
{
    class SimpleHBondFinder : IHBondFinder
    {
        protected Residue[] residues;
        private Point[] vecH;
        private Point[] vecN;
        protected Point[] vecCA;
        private Point[] vecC;
        private Point[] vecO;
        private bool[] canBeDonor;
        private bool[] canBeAcceptor;
        private double energyCutoff;

        public SimpleHBondFinder(IEnumerable<Residue> residues, double dsspEnergyCutoff)
        {
            this.residues = residues.ToArray();
            energyCutoff = dsspEnergyCutoff;

            canBeDonor = new bool[this.residues.Length];
            canBeAcceptor = new bool[this.residues.Length];
            vecCA = new Point[this.residues.Length];
            vecH = new Point[this.residues.Length];
            vecN = new Point[this.residues.Length];
            vecC = new Point[this.residues.Length];
            vecO = new Point[this.residues.Length];
            for (int i = 0; i < this.residues.Length; i++)
            {
                Residue r = this.residues[i];
                Atom? cAlpha = r.GetCAlpha();
                if (cAlpha == null)
                {
                    throw new ArgumentException($"Residue {r} has no C-alpha atom.");
                }
                vecCA[i] = ((Atom)r.GetCAlpha()).Position();
                Atom? hAmide = r.GetHAmide();
                Atom? nAmide = r.GetNAmide();
                Atom? cCarb = r.GetCCarb();
                Atom? oCarb = r.GetOCarb();
                if (hAmide != null && nAmide != null)
                {
                    canBeDonor[i] = true;
                    vecH[i] = ((Atom)hAmide).Position();
                    vecN[i] = ((Atom)nAmide).Position();
                }
                if (cCarb != null && oCarb != null)
                {
                    canBeAcceptor[i] = true;
                    vecC[i] = ((Atom)cCarb).Position();
                    vecO[i] = ((Atom)oCarb).Position();
                }
            }
        }

        /** Return DSSP energy between donor's NH and acceptor's CO. Return null, if some important atoms are missing and the energy cannot be calculated. */
        public double? DsspEnergy(int donor, int acceptor)
        {
            if (canBeDonor[donor] && canBeAcceptor[acceptor])
            {
                Point n = vecN[donor];
                Point h = vecH[donor];
                Point c = vecC[acceptor];
                Point o = vecO[acceptor];
                double energy = 0.084 * 332 * (1 / (o - n).Size + 1 / (c - h).Size - 1 / (o - h).Size - 1 / (c - n).Size);
                // Console.WriteLine($"{residues[donor]} -> {residues[acceptor]}:  {energy}");
                return energy;
            }
            else
            {
                return null;
            }
        }

        public bool IsHBond(int donor, int acceptor)
        {
            if (donor < 0 || acceptor < 0)
                return false;
            if (residues[donor].ChainId == residues[acceptor].ChainId && Math.Abs(residues[donor].SeqNumber - residues[acceptor].SeqNumber) <= 1)
            {
                return false;
            }
            double? dsspEnergy = DsspEnergy(donor, acceptor);
            if (dsspEnergy != null)
            {
                return dsspEnergy.Value <= energyCutoff;
            }
            else
            {
                //some residues do not contain needed atoms
                return false;
            }
        }

        public virtual List<int> FindHAcceptors(int res)
        {
            return Enumerable.Range(0, residues.Length).Where(s => IsHBond(res, s)).ToList();
        }

        public virtual List<int> FindHDonors(int res)
        {
            return Enumerable.Range(0, residues.Length).Where(s => IsHBond(s, res)).ToList();
        }
    }
}


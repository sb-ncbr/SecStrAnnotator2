using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace protein
{
	static class LibProtein{
		/**
         * Calculates the distance between two given atoms.
         */
		public static double Distance(Atom a1, Atom a2)
		{
			return Math.Sqrt (Math.Pow (a1.X - a2.X, 2) + Math.Pow (a1.Y - a2.Y, 2) + Math.Pow (a1.Z - a2.Z, 2));
		}

		/**
         * Selects all atoms which are closer to central atom than radius.
         */
		public static SortedSet<Atom> InRadius(SortedSet<Atom> atoms, Atom centralAtom, double radius)
		{
			SortedSet<Atom> inRadius = new SortedSet<Atom>();
			foreach (Atom atom in atoms)
			{
				if (Distance(atom, centralAtom) <= radius)
				{
					inRadius.Add(atom);
				}
			}
			return inRadius;
		}

		/**
         * Selects all residues that have any atom closer to central atom than radius.
         */
		public static SortedSet<Residue> ResiduesInRadiusByAny(ICollection<Residue> residues, Atom centralAtom, double radius)
		{
			return new SortedSet<Residue>(residues.Where(
				delegate(Residue res)
				{
					return res.GetAtoms().Any(
						delegate(Atom at)
						{
							return Distance(at, centralAtom) <= radius;
						});
				}));
		}

		/**
         * Selects all residues whose C alpha atom is closer to central atom than radius.
         */
		public static SortedSet<Residue> ResiduesInRadiusByCA(ICollection<Residue> residues, Atom centralAtom, double radius)
		{
			return new SortedSet<Residue>(residues.Where(
				delegate(Residue res)
				{
					return res.GetAtoms().Any(
						delegate(Atom at)
						{
							return at.Name.Equals(" CA ") && Distance(at, centralAtom) <= radius;
						});
				}));
		}

		/** Gets atoms of the protein backbone of residues residueFrom-residueTo in chain chainID in protein p. 
			The atoms are ordered in direction from N-terminus to C-terminus.
		*/
		public static List<Atom> Backbone(Protein p, char chainID, int residueFrom, int residueTo){
			List<Atom> result = new List<Atom> ();
			for (int i=residueFrom; i<=residueTo; i++) {
				Residue r;
				try {
					r = p.GetChain(chainID).GetResidues ().First (delegate (Residue res) {
						return res.ResSeq == i;
					});
				} catch (InvalidOperationException) {
					throw new ArgumentException ("Missing residue " + i + " in chain " + chainID + ".");
				}

				try {
					result.Add (r.GetAtoms ().First (delegate (Atom a) {
						return a.Name == Atom.NAME_N_AMIDE;
					}));
				} catch (InvalidOperationException) {
					throw new ArgumentException ("Missing N atom in residue " + i + ".");
				}

				try {
					result.Add (r.GetAtoms ().First (delegate (Atom a) {
						return a.Name == Atom.NAME_C_ALPHA;
					}));
				} catch (InvalidOperationException) {
					throw new ArgumentException ("Missing CA atom in residue " + i + ".");
				}

				try {
					result.Add (r.GetAtoms ().First (delegate (Atom a) {
						return a.Name == Atom.NAME_C_CARB;
					}));
				} catch (InvalidOperationException) {
					throw new ArgumentException ("Missing C atom in residue " + i + ".");
				}
			}
			return result;
		}

	}
}


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

	}
}


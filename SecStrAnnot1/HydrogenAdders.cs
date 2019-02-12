using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein
{
	public static class HydrogenAdders
	{
		public interface IHydrogenAdder {
			List<Residue> AddHydrogens (IEnumerable<Residue> residues);
		}

		public class DsspLikeAmideHydrogenAdder : IHydrogenAdder {
			const double NH_BOND_LENGTH = 1.0;

			public DsspLikeAmideHydrogenAdder(){
			}

			/**
			 * First removes all hydrogens and then adds amide hydrogens using simple approach described in DSSP.
			 * (Position of H is 1.0 Angstrom from N, in the opposite direction than carboxylic O is from C.)
			 */
			public List<Residue> AddHydrogens(IEnumerable<Residue> residues){
				throw new NotImplementedException();
				//TODO implement this somehow!
				/*int nextAtomID = residues.SelectMany(r => r.GetAtoms()).Max (a => a.Serial) + 1;
				List<Residue> newResidues = residues.Select(r => r.WithoutHydrogens()).ToList();
				List<Residue> relevantResidues = newResidues.Where (r => r.HasCAlpha ()).ToList();

				for (int i = 1; i < relevantResidues.Count; i++) {
					Residue r0 = relevantResidues [i - 1];
					Residue r1 = relevantResidues [i];
					if (r1.ChainID == r0.ChainID && r1.ResSeq == r0.ResSeq + 1 && !r1.IsProline) {
						try {
							Atom nAtom = r1.GetAtoms ().First (a => a.IsNAmide);
							Vector n = nAtom.Position ();
							Vector c = r0.GetAtoms ().First (a => a.IsCCarb).Position ();
							Vector o = r0.GetAtoms ().First (a => a.IsOCarb).Position ();
							Vector h = n + NH_BOND_LENGTH * (c - o).Normalize ();
							Atom newHydrogen = new Atom (nextAtomID++, Atom.NAME_H_AMIDE, nAtom.AltLoc, r1.Name, r1.ChainID, r1.ResSeq, nAtom.ICode, 
								h.X, h.Y, h.Z, nAtom.Occupancy, nAtom.TempFactor, Atom.ELEMENT_H, Atom.CHARGE_ZERO, false);
							r1.AddAtom(newHydrogen);
						} catch (InvalidOperationException){ // Some required atoms are missing
							//do nothing
						}
					}
				}
				return newResidues;*/
			}
		}

	}
}

